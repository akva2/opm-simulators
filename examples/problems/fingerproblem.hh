/*
  Copyright (C) 2008-2013 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \copydoc Ewoms::FingerProblem
 */
#ifndef EWOMS_FINGER_PROBLEM_HH
#define EWOMS_FINGER_PROBLEM_HH

#include "fingergridmanager.hh"

#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/ParkerLenhard.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <opm/material/fluidsystems/2pImmiscibleFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Air.hpp>

#include <ewoms/models/immiscible/immiscibleproperties.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <vector>
#include <string>

namespace Ewoms {
template <class TypeTag>
class FingerProblem;
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(FingerBaseProblem, INHERITS_FROM(FingerGridManager));

// declare the properties used by the finger problem
NEW_PROP_TAG(InitialWaterSaturation);

// Set the problem property
SET_TYPE_PROP(FingerBaseProblem, Problem, Ewoms::FingerProblem<TypeTag>);

// Set the wetting phase
SET_PROP(FingerBaseProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(FingerBaseProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::GasPhase<Scalar, Opm::Air<Scalar> > type;
};

// Set the material Law
SET_PROP(FingerBaseProblem, MaterialLaw)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx> Traits;

    // use the parker-lenhard hysteresis law
    typedef Opm::ParkerLenhard<Traits> ParkerLenhard;
    typedef ParkerLenhard type;
};

// Write the solutions of individual newton iterations?
SET_BOOL_PROP(FingerBaseProblem, NewtonWriteConvergence, false);

// Use forward differences instead of central differences
SET_INT_PROP(FingerBaseProblem, NumericDifferenceMethod, +1);

// Enable constraints
SET_INT_PROP(FingerBaseProblem, EnableConstraints, true);

// Enable gravity
SET_BOOL_PROP(FingerBaseProblem, EnableGravity, true);

// define the properties specific for the finger problem
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeX, 0.1);
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeY, 0.3);
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeZ, 0.1);

SET_SCALAR_PROP(FingerBaseProblem, InitialWaterSaturation, 0.01);

SET_INT_PROP(FingerBaseProblem, CellsX, 20);
SET_INT_PROP(FingerBaseProblem, CellsY, 70);
SET_INT_PROP(FingerBaseProblem, CellsZ, 1);

// The default for the end time of the simulation
SET_SCALAR_PROP(FingerBaseProblem, EndTime, 1e3);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(FingerBaseProblem, InitialTimeStepSize, 10);
}} // namespace Opm, Properties

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 * \brief Two-phase problem featuring some gravity-driven saturation
 *        fingers.
 *
 * The domain of this problem is sized 10cm times 1m and is initially
 * dry. Water is then injected at three locations on the top of the
 * domain which leads to gravity fingering. The boundary conditions
 * used are no-flow for the left and right and top of the domain and
 * free-flow at the bottom. This problem uses the Parker-Lenhard
 * hystersis model which might lead to non-monotonic saturation in the
 * fingers if the right material parameters is chosen and the spatial
 * discretization is fine enough.
 */
template <class TypeTag>
class FingerProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    //!\cond SKIP_THIS
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum {
        // number of phases

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiWettingEqIdx = Indices::conti0EqIdx + wettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;

    typedef typename GET_PROP(TypeTag, MaterialLaw)::ParkerLenhard ParkerLenhard;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    //!\endcond

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    FingerProblem(Simulator &simulator)
        : ParentType(simulator)
    {
        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 20; // -> 20°C
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::name
     */
    std::string name() const
    { return std::string("finger_") + Model::name(); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitialWaterSaturation,
                             "The initial saturation in the domain [] of the "
                             "wetting phase");
    }

    /*!
     * \copydoc VcfvProblem::init
     */
    void init()
    {

        // parameters for the Van Genuchten law of the main imbibition
        // and the main drainage curves.
        micParams_.setVgAlpha(0.0037);
        micParams_.setVgN(4.7);
        micParams_.finalize();

        mdcParams_.setVgAlpha(0.0037);
        mdcParams_.setVgN(4.7);
        mdcParams_.finalize();

        // initialize the material parameter objects of the individual
        // finite volumes
        int n = this->model().numDof();
        materialParams_.resize(n);
        for (int i = 0; i < n; ++i) {
            materialParams_[i].setMicParams(&micParams_);
            materialParams_[i].setMdcParams(&mdcParams_);
            materialParams_[i].setSwr(0.0);
            materialParams_[i].setSnr(0.1);
            materialParams_[i].finalize();
            ParkerLenhard::reset(materialParams_[i]);
        }

        K_ = this->toDimMatrix_(4.6e-10);

        setupInitialFluidState_();

        ParentType::init();
    }

    /*!
     * \copydoc VcfvProblem::postTimeStep
     */
    void postTimeStep()
    {
        // update the history of the hysteresis law
        ElementContext elemCtx(this->simulator());

        auto elemIt = this->gridView().template begin<0>();
        const auto &elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateAll(*elemIt);
            for (int scvIdx = 0; scvIdx < elemCtx.numDof(/*timeIdx=*/0); ++scvIdx) {
                int globalIdx = elemCtx.globalSpaceIndex(scvIdx, /*timeIdx=*/0);
                const auto &fs
                    = elemCtx.volVars(scvIdx, /*timeIdx=*/0).fluidState();
                ParkerLenhard::update(materialParams_[globalIdx], fs);
            }
        }
    }

    //! \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx,
                                           int timeIdx) const
    { return K_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    { return 0.4; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    {
        int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return materialParams_[globalSpaceIdx];
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.cvCenter(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)
            || onLowerBoundary_(pos)) {
            values.setNoFlow();
        }
        else {
            assert(onUpperBoundary_(pos));

            values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidState_);
        };

        // override the value for the liquid phase by forced
        // imbibition of water on inlet boundary segments
        if (onInlet_(pos)) {
            values[contiWettingEqIdx] = -0.001; // [kg/(m^2 s)]
        }
    }

    //! \}

    /*!
     * \name Volume terms
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx,
                 int timeIdx) const
    {
        // assign the primary variables
        values.assignNaive(initialFluidState_);
    }

    /*!
     * \copydoc VcfvProblem::constraints
     */
    template <class Context>
    void constraints(Constraints &constraints, const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        if (onUpperBoundary_(pos) && !onInlet_(pos)) {
            constraints.setAllConstraint();
            constraints.assignNaive(initialFluidState_);
        }
        else if (onLowerBoundary_(pos)) {
            constraints.setAllConstraint();
            constraints.assignNaive(initialFluidState_);
        }
    }

    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    { rate = Scalar(0.0); }
    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    {
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        Scalar lambda = (this->boundingBoxMax()[0] - pos[0]) / width;

        if (!onUpperBoundary_(pos))
            return false;

        Scalar xInject[] = { 0.25, 0.75 };
        Scalar injectLen[] = { 0.1, 0.1 };
        for (unsigned i = 0; i < sizeof(xInject) / sizeof(Scalar); ++i) {
            if (xInject[i] - injectLen[i] / 2 < lambda
                && lambda < xInject[i] + injectLen[i] / 2)
                return true;
        }
        return false;
    }

    void setupInitialFluidState_()
    {
        auto &fs = initialFluidState_;
        fs.setPressure(wettingPhaseIdx, /*pressure=*/1e5);

        Scalar Sw = EWOMS_GET_PARAM(TypeTag, Scalar, InitialWaterSaturation);
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

        // set the absolute pressures
        Scalar pn = 1e5;
        fs.setPressure(nonWettingPhaseIdx, pn);
        fs.setPressure(wettingPhaseIdx, pn);
    }

    DimMatrix K_;

    typename MaterialLawParams::VanGenuchtenParams micParams_;
    typename MaterialLawParams::VanGenuchtenParams mdcParams_;

    std::vector<MaterialLawParams> materialParams_;

    Opm::ImmiscibleFluidState<Scalar, FluidSystem> initialFluidState_;

    Scalar temperature_;
    Scalar eps_;
};

} // namespace Ewoms

#endif
