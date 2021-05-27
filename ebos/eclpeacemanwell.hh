// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/**
 * \file
 *
 * \copydoc Opm::EclPeacemanWell
 */
#ifndef EWOMS_ECL_PEACEMAN_WELL_HH
#define EWOMS_ECL_PEACEMAN_WELL_HH

#include <ebos/eclgenericpeacemanwell.hh>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/alignedallocator.hh>

#include <opm/simulators/wells/WGState.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

#include <map>
#include <unordered_set>

namespace Opm {

template <class TypeTag>
class EcfvDiscretization;

class WellState;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief The well model of Peaceman.
 *
 * This class is tailored for the element centered finite volume
 * discretization, assumes a vertical borehole and is intended to be
 * used by the EclWellManager.
 *
 * See:
 *
 * Z. Chen, G. Huan, Y. Ma: Computational Methods for Multiphase
 * Flows in Porous Media, 1st edition, SIAM, 2006, pp. 445-446
 *
 * and
 *
 * D. W. Peaceman: Interpretation of well-block pressures in numerical
 * reservoir simulation, The 52nd Annual SPE Fall Technical Conference
 * and Exhibition, Denver, CO., 1977
 */
template <class TypeTag>
class EclPeacemanWell : public BaseAuxiliaryModule<TypeTag>
                      , public EclGenericPeacemanWell<GetPropType<TypeTag, Properties::FluidSystem>,
                                                      GetPropType<TypeTag, Properties::Scalar>,
                                                      getPropValue<TypeTag, Properties::NumPhases>()>
{
    using AuxModule = BaseAuxiliaryModule<TypeTag>;

    using NeighborSet = typename AuxModule::NeighborSet;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Discretization = GetPropType<TypeTag, Properties::Discretization>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    typedef MathToolbox<Evaluation> Toolbox;

    typedef typename GridView::template Codim<0>::Entity        Element;
    typedef Element  ElementStorage;

    // the dimension of the simulator's world
    static const int dimWorld = GridView::dimensionworld;

    // convenient access to the number of phases and the number of
    // components
    static const unsigned numComponents = getPropValue<TypeTag, Properties::NumComponents>();
    static const unsigned numPhases = getPropValue<TypeTag, Properties::NumPhases>();

    using WellBase = EclGenericPeacemanWell<FluidSystem,Scalar,numPhases>;

    // convenient access to the phase and component indices. If the compiler bails out
    // here, you're probably using an incompatible fluid system. This class has only been
    // tested with Opm::FluidSystems::BlackOil...
    static const unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static const unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static const unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static const unsigned oilCompIdx = FluidSystem::oilCompIdx;
    static const unsigned waterCompIdx = FluidSystem::waterCompIdx;
    static const unsigned gasCompIdx = FluidSystem::gasCompIdx;

    static const unsigned numModelEq = getPropValue<TypeTag, Properties::NumEq>();
    static const unsigned conti0EqIdx = GetPropType<TypeTag, Properties::Indices>::conti0EqIdx;
    static const unsigned contiEnergyEqIdx = GetPropType<TypeTag, Properties::Indices>::contiEnergyEqIdx;

    static constexpr unsigned historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>();

    static constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();

    typedef CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/true> FluidState;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    // all quantities that need to be stored per degree of freedom that intersects the
    // well.
    struct DofVariables {
        DofVariables() = default;
        DofVariables(const DofVariables&) = default;

        // retrieve the solution dependent quantities that are only updated at the
        // beginning of a time step from the IntensiveQuantities of the model
        void updateBeginTimestep(const IntensiveQuantities&)
        {}

        // retrieve the solution dependent quantities from the IntensiveQuantities of the
        // model
        void update(const IntensiveQuantities& intQuants)
        {
            const auto& fs = intQuants.fluidState();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                pressure[phaseIdx] = fs.pressure(phaseIdx);
                density[phaseIdx] = fs.density(phaseIdx);
                mobility[phaseIdx] = intQuants.mobility(phaseIdx);
            }

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                oilMassFraction[compIdx] = fs.massFraction(oilPhaseIdx, compIdx);
                gasMassFraction[compIdx] = fs.massFraction(gasPhaseIdx, compIdx);
            }
        }

        // the depth of the centroid of the DOF
        Scalar depth;

        // the volume in m^3 of the DOF
        Scalar totalVolume;

        // the effective size of an element in each direction. This is defined as the
        // distance of the face centers along the respective axis.
        std::array<Scalar, dimWorld> effectiveSize;

        // the intrinsic permeability matrix for the degree of freedom
        DimMatrix permeability;

        // the effective permeability of the connection. usually that's the geometric
        // mean of the X and Y permeabilities of the DOF times the DOF's height
        Scalar effectivePermeability;

        // The connection transmissibility factor to be used for a given DOF. this is
        // usually computed from the values above but it can be explicitly specified by
        // the user...
        Scalar connectionTransmissibilityFactor;

        // the radius of the well for the given degree of freedom
        Scalar boreholeRadius;

        // The skin factor of the well at the given degree of freedom
        Scalar skinFactor;

        //////////////
        // the following quantities depend on the considered solution and are thus updated
        // at the beginning of each Newton-Raphson iteration.
        //////////////

        // the phase pressures inside a DOF
        std::array<Evaluation, numPhases> pressure;

        // the phase densities at the DOF
        std::array<Evaluation, numPhases> density;

        // the phase mobilities of the DOF
        std::array<Evaluation, numPhases> mobility;

        // the composition of the oil phase at the DOF
        std::array<Evaluation, numComponents> oilMassFraction;

        // the composition of the gas phase at the DOF
        std::array<Evaluation, numComponents> gasMassFraction;

        ElementStorage element;
        unsigned pvtRegionIdx;
        unsigned localDofIdx;
    };

    // some safety checks/caveats
    static_assert(std::is_same<Discretization, EcfvDiscretization<TypeTag> >::value,
                  "The Peaceman well model is only implemented for the "
                  "element-centered finite volume discretization!");
    static_assert(dimWorld == 3,
                  "The Peaceman well model is only implemented for 3D grids!");

public:
    EclPeacemanWell(const Simulator& simulator)
        : WellBase(simulator.vanguard().gridView().comm())
        , simulator_(simulator)
    {
        // set the composition of the injected fluids based. If
        // somebody is stupid enough to inject oil, we assume he wants
        // to loose his fortune on dry oil...
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
                injectionFluidState_.setMoleFraction(phaseIdx, compIdx, 0.0);
        injectionFluidState_.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);
        injectionFluidState_.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);
        injectionFluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, 1.0);

        // set the temperature to 25 deg C, just so that it is set
        injectionFluidState_.setTemperature(273.15 + 25);
    }

    /*!
     * \copydoc Opm::BaseAuxiliaryModule::numDofs()
     */
    unsigned numDofs() const override
    { return 1; }

    /*!
     * \copydoc Opm::BaseAuxiliaryModule::addNeighbors()
     */
    void addNeighbors(std::vector<NeighborSet>& neighbors) const override
    {
        int wellGlobalDof = AuxModule::localToGlobalDof(/*localDofIdx=*/0);

        // the well's bottom hole pressure always affects itself...
        neighbors[wellGlobalDof].insert(wellGlobalDof);

        // add the grid DOFs which are influenced by the well, and add the well dof to
        // the ones neighboring the grid ones
        auto wellDofIt = dofVariables_.begin();
        const auto& wellDofEndIt = dofVariables_.end();
        for (; wellDofIt != wellDofEndIt; ++ wellDofIt) {
            neighbors[wellGlobalDof].insert(wellDofIt->first);
            neighbors[wellDofIt->first].insert(wellGlobalDof);
        }
    }

    /*!
     * \copydoc Opm::BaseAuxiliaryModule::addNeighbors()
     */
    void applyInitial() override
    {
        auto& sol = const_cast<SolutionVector&>(simulator_.model().solution(/*timeIdx=*/0));

        int wellGlobalDof = AuxModule::localToGlobalDof(/*localDofIdx=*/0);
        sol[wellGlobalDof] = 0.0;

        // make valgrind shut up about the DOFs for the well even if the PrimaryVariables
        // class contains some "holes" due to alignment
        Valgrind::SetDefined(sol[wellGlobalDof]);

        // also apply the initial solution of the well to the "old" time steps
        for (unsigned timeIdx = 1; timeIdx < historySize; ++timeIdx) {
            auto& oldSol = const_cast<SolutionVector&>(simulator_.model().solution(timeIdx));

            oldSol[wellGlobalDof] = sol[wellGlobalDof];
        }
    }

    /*!
     * \copydoc Opm::BaseAuxiliaryModule::linearize()
     */
    void linearize(SparseMatrixAdapter& matrix, GlobalEqVector& residual) override
    {
        const SolutionVector& curSol = simulator_.model().solution(/*timeIdx=*/0);

        typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlock;

        unsigned wellGlobalDofIdx = AuxModule::localToGlobalDof(/*localDofIdx=*/0);
        residual[wellGlobalDofIdx] = 0.0;

        MatrixBlock diagBlock(0.0);
        for (unsigned i = 0; i < numModelEq; ++ i)
            diagBlock[i][i] = 1.0;

        MatrixBlock block(0.0);

        if (this->wellStatus() == WellBase::Shut) {
            // if the well is shut, make the auxiliary DOFs a trivial equation in the
            // matrix: the main diagonal is already set to the identity matrix, the
            // off-diagonal matrix entries must be set to 0.
            auto wellDofIt = dofVariables_.begin();
            const auto& wellDofEndIt = dofVariables_.end();
            for (; wellDofIt != wellDofEndIt; ++ wellDofIt) {
                matrix.setBlock(wellGlobalDofIdx, wellDofIt->first, block);
                matrix.setBlock(wellDofIt->first, wellGlobalDofIdx, block);
            }
            matrix.setBlock(wellGlobalDofIdx, wellGlobalDofIdx, diagBlock);
            residual[wellGlobalDofIdx] = 0.0;
            return;
        }
        else if (dofVariables_.empty()) {
            // the well does not feature any perforations on the local process
            matrix.setBlock(wellGlobalDofIdx, wellGlobalDofIdx, diagBlock);
            residual[wellGlobalDofIdx] = 0.0;
            return;
        }

        Scalar wellResid = wellResidual_(this->actualBottomHolePressure_);
        residual[wellGlobalDofIdx][0] = wellResid;

        // account for the effect of the grid DOFs which are influenced by the well on
        // the well equation and the effect of the well on the grid DOFs
        auto wellDofIt = dofVariables_.begin();
        const auto& wellDofEndIt = dofVariables_.end();

        ElementContext elemCtx(simulator_);
        for (; wellDofIt != wellDofEndIt; ++ wellDofIt) {
            unsigned gridDofIdx = wellDofIt->first;
            const auto& dofVars = *dofVariables_[gridDofIdx];
            DofVariables tmpDofVars(dofVars);
            auto priVars(curSol[gridDofIdx]);

            /////////////
            // influence of grid on well
            elemCtx.updateStencil( dofVars.element );
            // reset block from previous values
            block = 0.0;
            for (unsigned priVarIdx = 0; priVarIdx < numModelEq; ++priVarIdx) {
                // calculate the derivative of the well equation w.r.t. the current
                // primary variable using forward differences
                Scalar eps =
                    1e3
                    *std::numeric_limits<Scalar>::epsilon()
                    *std::max<Scalar>(1.0, priVars[priVarIdx]);
                priVars[priVarIdx] += eps;

                elemCtx.updateIntensiveQuantities(priVars, dofVars.localDofIdx, /*timeIdx=*/0);
                tmpDofVars.update(elemCtx.intensiveQuantities(dofVars.localDofIdx, /*timeIdx=*/0));

                Scalar dWellEq_dPV =
                    (wellResidual_(this->actualBottomHolePressure_, &tmpDofVars, gridDofIdx) - wellResid)
                    / eps;
                block[0][priVarIdx] = dWellEq_dPV;

                // go back to the original primary variables
                priVars[priVarIdx] -= eps;
            }
            matrix.setBlock(wellGlobalDofIdx, gridDofIdx, block);

            //
            /////////////

            /////////////
            // influence of well on grid:
            RateVector q(0.0);
            RateVector modelRate;
            std::array<Scalar, numPhases> resvRates;

            elemCtx.updateIntensiveQuantities(priVars, dofVars.localDofIdx, /*timeIdx=*/0);

            const auto& fluidState = elemCtx.intensiveQuantities(dofVars.localDofIdx, /*timeIdx=*/0).fluidState();

            // first, we need the source term of the grid for the slightly disturbed well.
            Scalar eps =
                1e3
                *std::numeric_limits<Scalar>::epsilon()
                *std::max<Scalar>(1e5, this->actualBottomHolePressure_);
            computeVolumetricDofRates_(resvRates, this->actualBottomHolePressure_ + eps, *dofVariables_[gridDofIdx]);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                modelRate.setVolumetricRate(fluidState, phaseIdx, resvRates[phaseIdx]);
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    q[compIdx] += modelRate[compIdx];
            }

            // then, we subtract the source rates for a undisturbed well.
            computeVolumetricDofRates_(resvRates, this->actualBottomHolePressure_, *dofVariables_[gridDofIdx]);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                modelRate.setVolumetricRate(fluidState, phaseIdx, resvRates[phaseIdx]);
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    q[compIdx] -= modelRate[compIdx];
            }

            // and finally, we divide by the epsilon to get the derivative
            for (unsigned eqIdx = 0; eqIdx < numModelEq; ++eqIdx)
                q[eqIdx] /= eps;

            // now we put this derivative into the right place in the Jacobian
            // matrix. This is a bit hacky because it assumes that the model uses a mass
            // rate for each component as its first conservation equation, but we require
            // the black-oil model for now anyway, so this should not be too much of a
            // problem...
            Valgrind::CheckDefined(q);
            block = 0.0;
            for (unsigned eqIdx = 0; eqIdx < numModelEq; ++ eqIdx)
                block[eqIdx][0] = - getValue(q[eqIdx])/dofVars.totalVolume;

            matrix.setBlock(gridDofIdx, wellGlobalDofIdx, block);

            //
            /////////////
        }

        // effect of changing the well's bottom hole pressure on the well equation
        Scalar eps =
            1e3
            *std::numeric_limits<Scalar>::epsilon()
            *std::max<Scalar>(1e7, this->targetBottomHolePressure_);
        Scalar wellResidStar = wellResidual_(this->actualBottomHolePressure_ + eps);
        diagBlock[0][0] = (wellResidStar - wellResid)/eps;

        matrix.setBlock(wellGlobalDofIdx, wellGlobalDofIdx, diagBlock);
    }

    Scalar volumetricSurfaceRateForConnection(int globalDofIdx, int phaseIdx) const
    {
        const DofVariables& dofVars = *dofVariables_.at(globalDofIdx);
        std::array<Scalar, numPhases> volumetricReservoirRates;
        computeVolumetricDofRates_(volumetricReservoirRates, this->actualBottomHolePressure_, dofVars);
        std::array<Scalar, numPhases> volumetricSurfaceRates;
        computeSurfaceRates_(volumetricSurfaceRates, volumetricReservoirRates, dofVars);
        return volumetricSurfaceRates[phaseIdx];
    }


    // reset the well to the initial state, i.e. remove all degrees of freedom...
    void clear()
    {
        dofVarsStore_.clear();
        dofVariables_.clear();
    }

    /*!
     * \brief Add a degree of freedom to the well.
     */
    template <class Context>
    void addDof(const Context& context, unsigned dofIdx)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        if (applies(globalDofIdx))
            // we already have this DOF in the well!
            return;

        const auto& dofPos = context.pos(dofIdx, /*timeIdx=*/0);

        dofVarsStore_.push_back(DofVariables());
        dofVariables_[globalDofIdx] = &dofVarsStore_.back();
        DofVariables& dofVars = *dofVariables_[globalDofIdx];
        this->wellTotalVolume_ += context.model().dofTotalVolume(globalDofIdx);

        dofVars.element = context.element();

        dofVars.localDofIdx = dofIdx;
        dofVars.pvtRegionIdx = context.problem().pvtRegionIndex(context, dofIdx, /*timeIdx=*/0);
        assert(dofVars.pvtRegionIdx == 0);

        // determine the size of the element
        dofVars.effectiveSize.fill(0.0);

        // we assume all elements to be hexahedrons!
        assert(context.element().subEntities(/*codim=*/dimWorld) == 8);

        const auto& refElem = Dune::ReferenceElements<Scalar, /*dim=*/3>::cube();

        // determine the current element's effective size
        const auto& elem = context.element();
        unsigned faceIdx = 0;
        unsigned numFaces = refElem.size(/*codim=*/1);
        for (; faceIdx < numFaces; ++faceIdx) {
            const auto& faceCenterLocal = refElem.position(faceIdx, /*codim=*/1);
            const auto& faceCenter = elem.geometry().global(faceCenterLocal);

            switch (faceIdx) {
            case 0:
                dofVars.effectiveSize[0] -= faceCenter[0];
                break;
            case 1:
                dofVars.effectiveSize[0] += faceCenter[0];
                break;
            case 2:
                dofVars.effectiveSize[1] -= faceCenter[1];
                break;
            case 3:
                dofVars.effectiveSize[1] += faceCenter[1];
                break;
            case 4:
                dofVars.depth += faceCenter[2];
                dofVars.effectiveSize[2] -= faceCenter[2];
                break;
            case 5:
                dofVars.depth += faceCenter[2];
                dofVars.effectiveSize[2] += faceCenter[2];
                break;
            }
        }

        // the volume associated with the DOF
        dofVars.totalVolume = context.model().dofTotalVolume(globalDofIdx);

        // the depth of the degree of freedom
        dofVars.depth /= 2;

        // default borehole radius: 1/2 foot
        dofVars.boreholeRadius = 0.3048/2;

        // default skin factor: 0
        dofVars.skinFactor = 0;

        // the intrinsic permeability tensor of the DOF
        const auto& K = context.problem().intrinsicPermeability(context, dofIdx, /*timeIdx=*/0);
        dofVars.permeability = K;

        // default the effective permeability: Geometric mean of the x and y components
        // of the intrinsic permeability of DOF times the DOF's height.
        assert(K[0][0] > 0);
        assert(K[1][1] > 0);
        dofVars.effectivePermeability =
            std::sqrt(K[0][0]*K[1][1])*dofVars.effectiveSize[2];

        // from that, compute the default connection transmissibility factor
        computeConnectionTransmissibilityFactor_(globalDofIdx);

        // we assume that the z-coordinate represents depth (and not
        // height) here...
        if (dofPos[2] < this->refDepth_)
            this->refDepth_ = dofPos[2];
    }

    int numConnections() const
    { return dofVariables_.size(); }

    /*!
     * \brief Set the connection transmissibility factor for a given degree of freedom.
     */
    template <class Context>
    void setConnectionTransmissibilityFactor(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx]->connectionTransmissibilityFactor = value;
    }

    /*!
     * \brief Set the effective permeability Kh to be used for a given degree of freedom.
     *
     * By default, Kh is sqrt(K_xx * K_yy) * h, where K_xx and K_yy is the permeability
     * for the DOF in X and Y directions and h is the height associated with the degree
     * of freedom.
     *
     * Note: The connection transmissibility factor is updated after calling this method,
     *       so if setConnectionTransmissibilityFactor() is to have any effect, it should
     *       be called after setEffectivePermeability()!
     */
    template <class Context>
    void setEffectivePermeability(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx].effectivePermeability = value;

        computeConnectionTransmissibilityFactor_(globalDofIdx);
    }

    /*!
     * \brief Return true iff a degree of freedom is directly affected
     *        by the well
     */
    bool applies(unsigned globalDofIdx) const
    { return dofVariables_.count(globalDofIdx) > 0; }

    /*!
     * \brief Set the skin factor of the well
     *
     * Note: The connection transmissibility factor is updated after calling this method,
     *       so if setConnectionTransmissibilityFactor() is to have any effect, it should
     *       be called after setSkinFactor()!
     */
    template <class Context>
    void setSkinFactor(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx].skinFactor = value;

        computeConnectionTransmissibilityFactor_(globalDofIdx);
    }

    /*!
     * \brief Return the well's skin factor at a DOF [-].
     */
    Scalar skinFactor(unsigned gridDofIdx) const
    { return dofVariables_.at(gridDofIdx).skinFactor_; }

    /*!
     * \brief Set the borehole radius of the well
     *
     * Note: The connection transmissibility factor is updated after calling this method,
     *       so if setConnectionTransmissibilityFactor() is to have any effect, it should
     *       be called after setRadius()!
     */
    template <class Context>
    void setRadius(const Context& context, unsigned dofIdx, Scalar value)
    {
        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx]->boreholeRadius = value;

        computeConnectionTransmissibilityFactor_(globalDofIdx);
    }

    /*!
     * \brief Return the well's radius at a cell [m].
     */
    Scalar radius(unsigned gridDofIdx) const
    { return dofVariables_.at(gridDofIdx)->radius_; }

    /*!
     * \brief Informs the well that an iteration has just begun.
     *
     * The beginIteration*() methods, the well calculates the bottom
     * and tubing head pressures, the actual unconstraint production and
     * injection rates, etc. The callback is split into three parts as
     * this arrangement avoids iterating over the whole grid and to
     * re-calculate the volume variables for each well.
     *
     * This is supposed to prepare the well object to do the
     * computations which are required to do the DOF specific
     * things.
     */
    void beginIterationPreProcess()
    { }

    /*!
     * \brief Do the DOF specific part at the beginning of each iteration
     */
    template <class Context>
    void beginIterationAccumulate(Context& context, unsigned timeIdx)
    {
        if (this->wellStatus() == WellBase::Shut)
            return;

        for (unsigned dofIdx = 0; dofIdx < context.numPrimaryDof(timeIdx); ++dofIdx) {
            unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
            if (!applies(globalDofIdx))
                continue;

            DofVariables& dofVars = *dofVariables_.at(globalDofIdx);
            const auto& intQuants = context.intensiveQuantities(dofIdx, timeIdx);

            if (this->iterationIdx_ == 0)
                dofVars.updateBeginTimestep(intQuants);

            dofVars.update(intQuants);
        }
    }

    /*!
     * \brief Informs the well that an iteration has just begun.
     *
     * This is the post-processing part which uses the results of the
     * accumulation callback.
     */
    void beginIterationPostProcess()
    {
        if (this->wellStatus() == WellBase::Shut)
            return;

        auto& sol = const_cast<SolutionVector&>(simulator_.model().solution(/*timeIdx=*/0));
        int wellGlobalDof = AuxModule::localToGlobalDof(/*localDofIdx=*/0);

        if (!dofVariables_.empty()) {
            // retrieve the bottom hole pressure from the global system of equations
            this->actualBottomHolePressure_ = Toolbox::value(dofVariables_.begin()->second->pressure[0]);
            this->actualBottomHolePressure_ = computeRateEquivalentBhp_();
        }
        else
            // start with 300 bars if we don't have anything better
            this->actualBottomHolePressure_ = 300 * 1e5;

        sol[wellGlobalDof][0] = this->actualBottomHolePressure_;

        computeOverallRates_(this->actualBottomHolePressure_,
                             this->actualResvRates_,
                             this->actualSurfaceRates_);

        this->actualWeightedResvRate_ = computeWeightedRate_(this->actualResvRates_);
        this->actualWeightedSurfaceRate_ = computeWeightedRate_(this->actualSurfaceRates_);
    }

    /*!
     * \brief Computes the source term for a degree of freedom.
     */
    template <class Context>
    void computeTotalRatesForDof(RateVector& q,
                                 const Context& context,
                                 unsigned dofIdx,
                                 unsigned timeIdx) const
    {
        q = 0.0;

        unsigned globalDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        if (this->wellStatus() == WellBase::Shut || !applies(globalDofIdx))
            return;

        // create a DofVariables object for the current evaluation point
        DofVariables tmp(*dofVariables_.at(globalDofIdx));

        tmp.update(context.intensiveQuantities(dofIdx, timeIdx));

        std::array<Evaluation, numPhases> volumetricRates;
        computeVolumetricDofRates_(volumetricRates, this->actualBottomHolePressure_, tmp);

        // convert to mass rates
        RateVector modelRate(0.0);
        const auto& intQuants = context.intensiveQuantities(dofIdx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            // energy is disabled or we have production for the given phase, i.e., we
            // can use the intensive quantities' fluid state
            modelRate.setVolumetricRate(intQuants.fluidState(), phaseIdx, volumetricRates[phaseIdx]);

            if (enableEnergy) {
                if (volumetricRates[phaseIdx] < 0.0) {
                    // producer
                    const auto& fs = intQuants.fluidState();
                    modelRate[contiEnergyEqIdx] += volumetricRates[phaseIdx]*fs.density(phaseIdx)*fs.enthalpy(phaseIdx);
                }
                else if (volumetricRates[phaseIdx] > 0.0
                         && this->injectedPhaseIdx_ == phaseIdx)
                {
                    // injector for the right phase. we need to use the thermodynamic
                    // quantities from the borehole as upstream
                    //
                    // TODO: This is not implemented in a very efficient way, the
                    // required quantities could be precomputed at initialization!
                    auto fs = injectionFluidState_;

                    // TODO: maybe we need to use a depth dependent pressure here. the
                    // difference is probably not very large, and for wells that span
                    // multiple perforations it is unclear what "well temperature" means
                    // anyway.
                    fs.setPressure(phaseIdx, this->actualBottomHolePressure_);

                    fs.setTemperature(this->wellTemperature_);

                    typename FluidSystem::template ParameterCache<Evaluation> paramCache;
                    unsigned globalSpaceIdx = context.globalSpaceIndex(dofIdx, timeIdx);
                    unsigned pvtRegionIdx = context.primaryVars(dofIdx, timeIdx).pvtRegionIndex();
                    paramCache.setRegionIndex(pvtRegionIdx);
                    paramCache.setMaxOilSat(context.problem().maxOilSaturation(globalSpaceIdx));
                    paramCache.updatePhase(fs, phaseIdx);

                    const auto& rho = FluidSystem::density(fs, paramCache, phaseIdx);
                    fs.setDensity(phaseIdx, rho);

                    const auto& h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
                    fs.setEnthalpy(phaseIdx, h);

                    modelRate[contiEnergyEqIdx] += volumetricRates[phaseIdx]*fs.density(phaseIdx)*fs.enthalpy(phaseIdx);
                }
            }

            for (unsigned eqIdx = 0; eqIdx < modelRate.size(); ++eqIdx)
                q[conti0EqIdx + eqIdx] += modelRate[conti0EqIdx + eqIdx];
        }

        Valgrind::CheckDefined(q);
    }

protected:
    // compute the connection transmissibility factor based on the effective permeability
    // of a connection, the radius of the borehole and the skin factor.
    void computeConnectionTransmissibilityFactor_(unsigned globalDofIdx)
    {
        auto& dofVars = *dofVariables_[globalDofIdx];

        const auto& D = dofVars.effectiveSize;
        const auto& K = dofVars.permeability;
        Scalar Kh = dofVars.effectivePermeability;
        Scalar S = dofVars.skinFactor;
        Scalar rWell = dofVars.boreholeRadius;

        // compute the "equivalence radius" r_0 of the connection
        assert(K[0][0] > 0.0);
        assert(K[1][1] > 0.0);
        Scalar tmp1 = std::sqrt(K[1][1]/K[0][0]);
        Scalar tmp2 = 1.0 / tmp1;
        Scalar r0 = std::sqrt(D[0]*D[0]*tmp1 + D[1]*D[1]*tmp2);
        r0 /= std::sqrt(tmp1) + std::sqrt(tmp2);
        r0 *= 0.28;

        // we assume the well borehole in the center of the dof and that it is vertical,
        // i.e., the area which is exposed to the flow is 2*pi*r0*h. (for non-vertical
        // wells this would need to be multiplied with the cosine of the angle and the
        // height must be adapted...)
        const Scalar exposureFactor = 2*M_PI;

        dofVars.connectionTransmissibilityFactor = exposureFactor*Kh/(std::log(r0 / rWell) + S);
    }

    template <class ResultEval, class BhpEval>
    void computeVolumetricDofRates_(std::array<ResultEval, numPhases>& volRates,
                                    const BhpEval& bottomHolePressure,
                                    const DofVariables& dofVars) const
    {
        typedef MathToolbox<Evaluation> DofVarsToolbox;
        typedef typename std::conditional<std::is_same<BhpEval, Scalar>::value,
                                          ResultEval,
                                          Scalar>::type DofEval;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            volRates[phaseIdx] = 0.0;

        // connection transmissibility factor for the current DOF.
        Scalar Twj = dofVars.connectionTransmissibilityFactor;

        // bottom hole pressure and depth of the degree of freedom
        ResultEval pbh = bottomHolePressure;
        Scalar depth = dofVars.depth;

        // gravity constant
        Scalar g = simulator_.problem().gravity()[dimWorld - 1];

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            // well model due to Peaceman; see Chen et al., p. 449

            // phase pressure in grid cell
            const DofEval& p = DofVarsToolbox::template decay<DofEval>(dofVars.pressure[phaseIdx]);

            // density and mobility of fluid phase
            const DofEval& rho = DofVarsToolbox::template decay<DofEval>(dofVars.density[phaseIdx]);
            DofEval lambda;
            if (this->wellType_ == WellBase::Producer) {
                //assert(p < pbh);
                lambda = DofVarsToolbox::template decay<DofEval>(dofVars.mobility[phaseIdx]);
            }
            else if (this->wellType_ == WellBase::Injector) {
                //assert(p > pbh);
                if (phaseIdx != this->injectedPhaseIdx_)
                    continue;

                // use the total mobility, i.e. the sum of all phase mobilities at the
                // injector cell. this seems a bit weird: at the wall of the borehole,
                // there should only be injected phase present, so its mobility should be
                // 1/viscosity...
                lambda = 0.0;
                for (unsigned phase2Idx = 0; phase2Idx < numPhases; ++phase2Idx) {
                    if (!FluidSystem::phaseIsActive(phase2Idx))
                        continue;

                    lambda += DofVarsToolbox::template decay<DofEval>(dofVars.mobility[phase2Idx]);
                }
            }
            else
                throw std::logic_error("Type of well \""+this->name()+"\" is undefined");

            Valgrind::CheckDefined(pbh);
            Valgrind::CheckDefined(p);
            Valgrind::CheckDefined(g);
            Valgrind::CheckDefined(rho);
            Valgrind::CheckDefined(lambda);
            Valgrind::CheckDefined(depth);
            Valgrind::CheckDefined(this->refDepth_);

            // pressure in the borehole ("hole pressure") at the given location
            ResultEval ph = pbh + rho*g*(depth - this->refDepth_);

            // volumetric reservoir rate for the phase
            volRates[phaseIdx] = Twj*lambda*(ph - p);

            Valgrind::CheckDefined(g);
            Valgrind::CheckDefined(ph);
            Valgrind::CheckDefined(volRates[phaseIdx]);
        }
    }

    /*!
     * \brief Given the volumetric rates for all phases, return the
     *        corresponding weighted rate
     *
     * The weights are user-specified and can be set using
     * setVolumetricPhaseWeights()
     */
    template <class Eval>
    Eval computeWeightedRate_(const std::array<Eval, numPhases>& volRates) const
    {
        Eval result = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            result += volRates[phaseIdx]*this->volumetricWeight_[phaseIdx];
        }
        return result;
    }

    /*!
     * \brief Convert volumetric reservoir rates into volumetric volume rates.
     *
     * This requires the density and composition of the phases and
     * thus the applicable fluid state.
     */
    template <class Eval>
    void computeSurfaceRates_(std::array<Eval, numPhases>& surfaceRates,
                              const std::array<Eval, numPhases>& reservoirRate,
                              const DofVariables& dofVars) const
    {
        // the array for the surface rates and the one for the reservoir rates must not
        // be the same!
        assert(&surfaceRates != &reservoirRate);

        int regionIdx = dofVars.pvtRegionIdx;

        // If your compiler bails out here, you have not chosen the correct fluid
        // system. Currently, only Opm::FluidSystems::BlackOil is supported, sorry...
        Scalar rhoOilSurface = FluidSystem::referenceDensity(oilPhaseIdx, regionIdx);
        Scalar rhoGasSurface = FluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rhoWaterSurface = FluidSystem::referenceDensity(waterPhaseIdx, regionIdx);

        // oil
        if (FluidSystem::phaseIsActive(oilPhaseIdx))
            surfaceRates[oilPhaseIdx] =
                // oil in gas phase
                reservoirRate[gasPhaseIdx]
                * Toolbox::value(dofVars.density[gasPhaseIdx])
                * Toolbox::value(dofVars.gasMassFraction[oilCompIdx])
                / rhoOilSurface
                +
                // oil in oil phase
                reservoirRate[oilPhaseIdx]
                * Toolbox::value(dofVars.density[oilPhaseIdx])
                * Toolbox::value(dofVars.oilMassFraction[oilCompIdx])
                / rhoOilSurface;

        // gas
        if (FluidSystem::phaseIsActive(gasPhaseIdx))
            surfaceRates[gasPhaseIdx] =
                // gas in gas phase
                reservoirRate[gasPhaseIdx]
                * Toolbox::value(dofVars.density[gasPhaseIdx])
                * Toolbox::value(dofVars.gasMassFraction[gasCompIdx])
                / rhoGasSurface
                +
                // gas in oil phase
                reservoirRate[oilPhaseIdx]
                * Toolbox::value(dofVars.density[oilPhaseIdx])
                * Toolbox::value(dofVars.oilMassFraction[gasCompIdx])
                / rhoGasSurface;

        // water
        if (FluidSystem::phaseIsActive(waterPhaseIdx))
            surfaceRates[waterPhaseIdx] =
                reservoirRate[waterPhaseIdx]
                * Toolbox::value(dofVars.density[waterPhaseIdx])
                / rhoWaterSurface;
    }

    const WellState& wellState() const
    {
        throw std::logic_error("wellState() method not implemented for class eclpeacemanwell");
    }

    WellState& wellState()
    {
        throw std::logic_error("wellState() method not implemented for class eclpeacemanwell");
    }

    void commitWGState()
    {
        throw std::logic_error("commitWellState() method not implemented for class eclpeacemanwell");
    }

    void commitWGState(WGState)
    {
        throw std::logic_error("commitWellState() method not implemented for class eclpeacemanwell");
    }

    void resetWGState()
    {
        throw std::logic_error("resetWellState() method not implemented for class eclpeacemanwell");
    }

    void updateNupcolWGState()
    {
        throw std::logic_error("updateNupcolWellState() method not implemented for class eclpeacemanwell");
    }

    void
    updateEclWell(int, int)
    {
        throw std::logic_error("updateEclWell() method not implemented for class eclpeacemanwell");
    }


    void
    updateEclWells(int, const std::unordered_set<std::string>&) {
        throw std::logic_error("updateEclWells() method not implemented for class eclpeacemanwell");
    }


    double
    wellPI(int) const
    {
        throw std::logic_error("wellPI() method not implemented for class eclpeacemanwell");
    }

    double
    wellPI(const std::string& ) const
    {
        throw std::logic_error("wellPI() method not implemented for class eclpeacemanwell");
    }



    /*!
     * \brief Compute the volumetric phase rate of the complete well given a bottom hole
     *        pressure.
     *
     * A single degree of freedom may be different from the evaluation point.
     */
    void computeOverallRates_(Scalar bottomHolePressure,
                              std::array<Scalar, numPhases>& overallResvRates,
                              std::array<Scalar, numPhases>& overallSurfaceRates,
                              const DofVariables *evalDofVars = 0,
                              int globalEvalDofIdx = -1) const

    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            overallResvRates[phaseIdx] = 0.0;
            overallSurfaceRates[phaseIdx] = 0.0;
        }

        auto dofVarsIt = dofVariables_.begin();
        const auto& dofVarsEndIt = dofVariables_.end();
        for (; dofVarsIt != dofVarsEndIt; ++ dofVarsIt) {
            std::array<Scalar, numPhases> volumetricReservoirRates;
            const DofVariables *tmp;
            if (dofVarsIt->first == globalEvalDofIdx)
                tmp = evalDofVars;
            else
                tmp = dofVarsIt->second;

            computeVolumetricDofRates_<Scalar, Scalar>(volumetricReservoirRates, bottomHolePressure, *tmp);

            std::array<Scalar, numPhases> volumetricSurfaceRates;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                volumetricSurfaceRates[ phaseIdx ] = 0;
            }
            computeSurfaceRates_(volumetricSurfaceRates, volumetricReservoirRates, *tmp);

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                overallResvRates[phaseIdx] += volumetricReservoirRates[phaseIdx];
                overallSurfaceRates[phaseIdx] += volumetricSurfaceRates[phaseIdx];
            }
        }
    }

    /*!
     * \brief Compute the weighted volumetric rate of the complete well given a bottom
     *        hole pressure.
     *
     * A single degree of freedom may be different from the evaluation point.
     */
    Scalar computeOverallWeightedSurfaceRate_(Scalar bottomHolePressure,
                                              std::array<Scalar, numPhases>& overallSurfaceRates,
                                              const DofVariables& evalDofVars,
                                              int globalEvalDofIdx) const

    {
        static std::array<Scalar, numPhases> resvRatesDummy;
        computeOverallRates_(bottomHolePressure,
                             overallSurfaceRates,
                             resvRatesDummy,
                             evalDofVars,
                             globalEvalDofIdx);
        return computeWeightedRate_(overallSurfaceRates);
    }

    // this is a more convenient version of the method above if all degrees of freedom
    // are supposed to be at their evaluation points.
    Scalar computeOverallWeightedSurfaceRate_(Scalar bottomHolePressure,
                                              std::array<Scalar, numPhases>& overallSurfaceRates) const
    {
        // create a dummy DofVariables object and call the method above using an index
        // that is guaranteed to never be part of a well...
        static DofVariables dummyDofVars;
        return computeOverallWeightedSurfaceRate_(bottomHolePressure,
                                                  overallSurfaceRates,
                                                  dummyDofVars,
                                                  /*globalEvalDofIdx=*/-1);
    }

    /*!
     * \brief Compute the "rate-equivalent bottom hole pressure"
     *
     * I.e. The bottom hole pressure where the well rate is exactly the one which is
     * targeted. This is zero of the "rate-equivalent bottom hole pressure" would be
     * smaller than 1 bar.
     */
    Scalar computeRateEquivalentBhp_() const
    {
        if (this->wellStatus() == WellBase::Shut)
            // there is no flow happening in the well, so we return 0...
            return 0.0;

        // initialize the bottom hole pressure which we would like to calculate
        Scalar bhpScalar = this->actualBottomHolePressure_;
        if (bhpScalar > 1e8)
            bhpScalar = 1e8;
        if (bhpScalar < 1e5)
            bhpScalar = 1e5;

        // if the BHP goes below 1 bar for the first time, we reset it to 10 bars and
        // are "on bail", i.e. if it goes below 1 bar again, we give up because the
        // final pressure would be below 1 bar...
        bool onBail = false;

        // Newton-Raphson method
        typedef DenseAd::Evaluation<Scalar, 1> BhpEval;

        BhpEval bhpEval(bhpScalar);
        bhpEval.setDerivative(0, 1.0);
        const Scalar tolerance = 1e3*std::numeric_limits<Scalar>::epsilon();
        const int maxIter = 20;
        for (int iterNum = 0; iterNum < maxIter; ++iterNum) {
            const auto& f = wellResidual_<BhpEval>(bhpEval);

            if (std::abs(f.derivative(0)) < 1e-20)
                throw NumericalIssue("Cannot determine the bottom hole pressure for well "+this->name()
                                     +": Derivative of the well residual is too small");
            Scalar delta = f.value()/f.derivative(0);

            bhpEval.setValue(bhpEval.value() - delta);
            if (bhpEval < 1e5) {
                bhpEval.setValue(1e5);
                if (onBail)
                    return bhpEval.value();
                else
                    onBail = true;
            }
            else
                onBail = false;

            if (std::abs(delta/bhpEval.value()) < tolerance)
                return bhpEval.value();
        }

        throw NumericalIssue("Could not determine the bottom hole pressure of well '"+this->name()
                              +"' within " + std::to_string(maxIter) + " iterations.");
    }

    template <class BhpEval>
    BhpEval wellResidual_(const BhpEval& bhp,
                          const DofVariables *replacementDofVars = 0,
                          int replacedGridIdx = -1) const
    {
        typedef MathToolbox<BhpEval> BhpEvalToolbox;

        // compute the volumetric reservoir and surface rates for the complete well
        BhpEval resvRate = 0.0;

        std::array<BhpEval, numPhases> totalSurfaceRates;
        std::fill(totalSurfaceRates.begin(), totalSurfaceRates.end(), 0.0);

        auto dofVarsIt = dofVariables_.begin();
        const auto& dofVarsEndIt = dofVariables_.end();
        for (; dofVarsIt != dofVarsEndIt; ++ dofVarsIt) {
            std::array<BhpEval, numPhases> resvRates;
            const DofVariables *dofVars = dofVarsIt->second;
            if (replacedGridIdx == dofVarsIt->first)
                dofVars = replacementDofVars;
            computeVolumetricDofRates_(resvRates, bhp, *dofVars);

            std::array<BhpEval, numPhases> surfaceRates;
            computeSurfaceRates_(surfaceRates, resvRates, *dofVars);

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                totalSurfaceRates[phaseIdx] += surfaceRates[phaseIdx];
            }

            resvRate += computeWeightedRate_(resvRates);
        }

        BhpEval surfaceRate = computeWeightedRate_(totalSurfaceRates);

        // compute the residual of well equation. we currently use max(rateMax - rate,
        // bhp - targetBhp) for producers and max(rateMax - rate, bhp - targetBhp) for
        // injectors. (i.e., the target bottom hole pressure is an upper limit for
        // injectors and a lower limit for producers.) Note that with this approach, one
        // of the limits must always be reached to get the well equation to zero...
        Valgrind::CheckDefined(this->maximumSurfaceRate_);
        Valgrind::CheckDefined(this->maximumReservoirRate_);
        Valgrind::CheckDefined(surfaceRate);
        Valgrind::CheckDefined(resvRate);

        BhpEval result = 1e30;

        BhpEval maxSurfaceRate = this->maximumSurfaceRate_;
        BhpEval maxResvRate = this->maximumReservoirRate_;
        if (this->wellStatus() == WellBase::Closed) {
            // make the weight of the fluids on the surface equal and require that no
            // fluids are produced on the surface...
            maxSurfaceRate = 0.0;
            surfaceRate = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                surfaceRate += totalSurfaceRates[phaseIdx];
            }

            // don't care about the reservoir rate...
            maxResvRate = 1e30;
        }

        if (this->wellType_ == WellBase::Injector) {
            // for injectors the computed rates are positive and the target BHP is the
            // maximum allowed pressure ...
            result = BhpEvalToolbox::min(maxSurfaceRate - surfaceRate, result);
            result = BhpEvalToolbox::min(maxResvRate - resvRate, result);
            result = BhpEvalToolbox::min(1e-7*(this->targetBottomHolePressure_ - bhp), result);
        }
        else {
            assert(this->wellType_ == WellBase::Producer);
            // ... for producers the rates are negative and the bottom hole pressure is
            // is the minimum
            result = BhpEvalToolbox::min(maxSurfaceRate + surfaceRate, result);
            result = BhpEvalToolbox::min(maxResvRate + resvRate, result);
            result = BhpEvalToolbox::min(1e-7*(bhp - this->targetBottomHolePressure_), result);
        }

        const Scalar scalingFactor = 1e-3;
        return scalingFactor*result;
    }

    const Simulator& simulator_;

    std::vector<DofVariables, aligned_allocator<DofVariables, alignof(DofVariables)> > dofVarsStore_;
    std::map<int, DofVariables*> dofVariables_;

    // The thermodynamic state of the fluid which gets injected
    //
    // The fact that this attribute is mutable is kind of an hack
    // which can be avoided using a PressureOverlayFluidState, but
    // then performance would be slightly worse...
    mutable FluidState injectionFluidState_;
};

} // namespace Opm

#endif
