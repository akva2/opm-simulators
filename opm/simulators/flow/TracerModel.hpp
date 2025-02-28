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
 * \copydoc Opm::TracerModel
 */
#ifndef OPM_TRACER_MODEL_HPP
#define OPM_TRACER_MODEL_HPP

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/GenericTracerModel.hpp>
#include <opm/simulators/utils/VectorVectorDataHandle.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableTracerModel {
    using type = UndefinedProperty;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief A class which handles tracers as specified in by ECL
 */
template <class TypeTag>
class TracerModel : public GenericTracerModel<GetPropType<TypeTag, Properties::Grid>,
                                              GetPropType<TypeTag, Properties::GridView>,
                                              GetPropType<TypeTag, Properties::DofMapper>,
                                              GetPropType<TypeTag, Properties::Stencil>,
                                              GetPropType<TypeTag, Properties::FluidSystem>,
                                              GetPropType<TypeTag, Properties::Scalar>>
{
    using BaseType = GenericTracerModel<GetPropType<TypeTag, Properties::Grid>,
                                        GetPropType<TypeTag, Properties::GridView>,
                                        GetPropType<TypeTag, Properties::DofMapper>,
                                        GetPropType<TypeTag, Properties::Stencil>,
                                        GetPropType<TypeTag, Properties::FluidSystem>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using TracerEvaluation = DenseAd::Evaluation<Scalar,1>;

    using TracerMatrix = typename BaseType::TracerMatrix;
    using TracerVector = typename BaseType::TracerVector;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

public:
    explicit TracerModel(Simulator& simulator)
        : BaseType(simulator.vanguard().gridView(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().cartesianIndexMapper(),
                   simulator.model().dofMapper(),
                   simulator.vanguard().cellCentroids())
        , simulator_(simulator)
        , tbatch({waterPhaseIdx, oilPhaseIdx, gasPhaseIdx})
        , wat_(tbatch[0])
        , oil_(tbatch[1])
        , gas_(tbatch[2])
    { }


    /*
      The initialization of the tracer model is a three step process:

      1. The init() method is called. This will allocate buffers and initialize
         some phase index stuff. If this is a normal run the initial tracer
         concentrations will be assigned from the TBLK or TVDPF keywords.

      2. [Restart only:] The tracer concentration are read from the restart
         file and the concentrations are applied with repeated calls to the
         setTracerConcentration() method. This is currently done in the
         eclwriter::beginRestart() method.

      3. Internally the tracer model manages the concentrations in "batches" for
         the oil, water and gas tracers respectively. The batches should be
         initialized with the initial concentration, that must be performed
         after the concentration values have been assigned. This is done in
         method prepareTracerBatches() called from eclproblem::finishInit().
    */
    void init(bool rst)
    {
        this->doInit(rst, simulator_.model().numGridDof(),
                     gasPhaseIdx, oilPhaseIdx, waterPhaseIdx);
    }

    void prepareTracerBatches()
    {
        for (std::size_t tracerIdx = 0; tracerIdx < this->tracerPhaseIdx_.size(); ++tracerIdx) {
            if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::waterPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)){
                    throw std::runtime_error("Water tracer specified for non-water fluid system:" + this->name(tracerIdx));
                }

                wat_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }
            else if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::oilPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)){
                    throw std::runtime_error("Oil tracer specified for non-oil fluid system:" + this->name(tracerIdx));
                }

                oil_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }
            else if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::gasPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)){
                    throw std::runtime_error("Gas tracer specified for non-gas fluid system:" + this->name(tracerIdx));
                }

                gas_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }

            // resize free and solution volume storages
            vol1_[this->tracerPhaseIdx_[tracerIdx]][Free].
                resize(this->splitTracerConcentration_[Free][tracerIdx].size());
            vol1_[this->tracerPhaseIdx_[tracerIdx]][Solution].
                resize(this->splitTracerConcentration_[Solution][tracerIdx].size());
            dVol_[this->tracerPhaseIdx_[tracerIdx]][Free].
                resize(this->splitTracerConcentration_[Free][tracerIdx].size());
            dVol_[this->tracerPhaseIdx_[tracerIdx]][Solution].
                resize(this->splitTracerConcentration_[Solution][tracerIdx].size());
        }

        // will be valid after we move out of tracerMatrix_
        TracerMatrix* base = this->tracerMatrix_.get();
        for (auto& tr : this->tbatch) {
            if (tr.numTracer() != 0) {
                if (this->tracerMatrix_) {
                    tr.mat = std::move(this->tracerMatrix_);
                }
                else {
                    tr.mat = std::make_unique<TracerMatrix>(*base);
                }
            }
        }
    }

    void beginTimeStep()
    {
        if (this->numTracers() == 0) {
            return;
        }

        OPM_TIMEBLOCK(tracerUpdateCache);
        updateStorageCache();
    }

    /*!
     * \brief Informs the tracer model that a time step has just been finished.
     */
    void endTimeStep()
    {
        if (this->numTracers() == 0) {
            return;
        }

        OPM_TIMEBLOCK(tracerAdvance);
        advanceTracerFields();
    }

    /*!
     * \brief This method writes the complete state of all tracer
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter&)
    { /* not implemented */ }

    /*!
     * \brief This method restores the complete state of the tracer
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter&)
    { /* not implemented */ }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<BaseType&>(*this));
        serializer(tbatch);
    }

protected:
    using TracerTypeIdx = typename BaseType::TracerTypeIdx;
    using BaseType::Free;
    using BaseType::Solution;

    template<TracerTypeIdx Index>
    Scalar computeVolume_(const int tracerPhaseIdx,
                          const unsigned globalDofIdx,
                          const unsigned timeIdx) const
    {
        const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        Scalar phaseVolume;
        if constexpr (Index == Free) {
            phaseVolume = decay<Scalar>(fs.saturation(tracerPhaseIdx)) *
                          decay<Scalar>(fs.invB(tracerPhaseIdx)) *
                          decay<Scalar>(intQuants.porosity());
        } else {
            // vaporized oil
            if (tracerPhaseIdx == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                phaseVolume =
                    decay<Scalar>(fs.saturation(FluidSystem::gasPhaseIdx)) *
                    decay<Scalar>(fs.invB(FluidSystem::gasPhaseIdx)) *
                    decay<Scalar>(fs.Rv()) *
                    decay<Scalar>(intQuants.porosity());
            }

            // dissolved gas
            else if (tracerPhaseIdx == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                phaseVolume =
                    decay<Scalar>(fs.saturation(FluidSystem::oilPhaseIdx)) *
                    decay<Scalar>(fs.invB(FluidSystem::oilPhaseIdx)) *
                    decay<Scalar>(fs.Rs()) *
                    decay<Scalar>(intQuants.porosity());
            }
            else {
                phaseVolume = 0.0;
            }
        }

        return max(phaseVolume, 1e-10);
    }

    template<TracerTypeIdx Index>
    void computeFlux_(TracerEvaluation& flux,
                      bool& isUp,
                      const int tracerPhaseIdx,
                      const ElementContext& elemCtx,
                      const unsigned scvfIdx,
                      const unsigned timeIdx) const
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const unsigned inIdx = extQuants.interiorIndex();

        Scalar v;
        unsigned upIdx;
        if constexpr (Index == Free) {
            upIdx = extQuants.upstreamIndex(tracerPhaseIdx);
            const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
            const auto& fs = intQuants.fluidState();
            v = decay<Scalar>(extQuants.volumeFlux(tracerPhaseIdx)) *
                decay<Scalar>(fs.invB(tracerPhaseIdx));
        } else {
            if (tracerPhaseIdx == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                upIdx = extQuants.upstreamIndex(FluidSystem::gasPhaseIdx);

                const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
                const auto& fs = intQuants.fluidState();
                v = decay<Scalar>(fs.invB(FluidSystem::gasPhaseIdx)) *
                    decay<Scalar>(extQuants.volumeFlux(FluidSystem::gasPhaseIdx)) *
                    decay<Scalar>(fs.Rv());
            }
            // dissolved gas
            else if (tracerPhaseIdx == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                upIdx = extQuants.upstreamIndex(FluidSystem::oilPhaseIdx);
                const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
                const auto& fs = intQuants.fluidState();
                v = decay<Scalar>(fs.invB(FluidSystem::oilPhaseIdx)) *
                    decay<Scalar>(extQuants.volumeFlux(FluidSystem::oilPhaseIdx)) *
                    decay<Scalar>(fs.Rs());
            }
            else {
                upIdx = 0;
                v = 0.0;
            }
        }

        const Scalar A = scvf.area();
        if (inIdx == upIdx) {
            flux = A*v*variable<TracerEvaluation>(1.0, 0);
            isUp = true;
        }
        else {
            flux = A*v;
            isUp = false;
        }
    }

    template<TracerTypeIdx Index, class TrRe>
    void assembleTracerEquationVolume_(TrRe& tr,
                                       const ElementContext& elemCtx,
                                       const Scalar scvVolume,
                                       const Scalar dt,
                                       unsigned I,
                                       unsigned I1)
    {
        // Storage terms at previous time step (timeIdx = 1)
        auto storage1 = [&tr, this, &I, &I1,
                         cache = elemCtx.enableStorageCache()](const unsigned tIdx)
        {
            if (cache) {
                return tr.storageOfTimeIndex1_[tIdx][I][Index];
            }  else {
                const Scalar volume = computeVolume_<Index>(tr.phaseIdx_, I1, 1);
                return volume * tr.concentration_[tIdx][I][Index];
            }
        };

        const Scalar scdt = scvVolume / dt;

        const TracerEvaluation vol = computeVolume_<Index>(tr.phaseIdx_, I, 0) * variable<TracerEvaluation>(1.0, 0);
        dVol_[tr.phaseIdx_][I][Index] += vol.value() * scvVolume - vol1_[tr.phaseIdx_][I][Index];
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            const Scalar storageOfTimeIndex0 = vol.value() * tr.concentration_[tIdx][I][Index];
            const Scalar localStorage = (storageOfTimeIndex0 - storage1(tIdx)) * scdt;
            tr.residual_[tIdx][I][Index] += localStorage; // residual + flux
        }

        // Derivative matrix
        (*tr.mat)[I][I][Index][Index] += vol.derivative(0) * scdt;
    }

    template<class TrRe>
    void assembleTracerEquationVolume(TrRe& tr,
                                      const ElementContext& elemCtx,
                                      const Scalar scvVolume,
                                      const Scalar dt,
                                      unsigned I,
                                      unsigned I1)

    {
        if (tr.numTracer() == 0) {
            return;
        }

        assembleTracerEquationVolume_<Free>(tr, elemCtx, scvVolume, dt, I, I1);
        assembleTracerEquationVolume_<Solution>(tr, elemCtx, scvVolume, dt, I, I1);
    }

    template<TracerTypeIdx Index, class TrRe>
    void assembleTracerEquationFlux_(TrRe& tr,
                                     const ElementContext& elemCtx,
                                     unsigned scvfIdx,
                                     unsigned I,
                                     unsigned J,
                                     const Scalar dt)
    {
        TracerEvaluation flux;
        bool isUp;
        computeFlux_<Index>(flux, isUp, tr.phaseIdx_, elemCtx, scvfIdx, 0);
        dVol_[tr.phaseIdx_][I][Index] += flux.value() * dt;
        const int globalUpIdx = isUp ? I : J;
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            // Fluxes
            tr.residual_[tIdx][I][Index] += flux.value() *
                                            tr.concentration_[tIdx][globalUpIdx][Index]; // residual + flux
        }

        // Derivative matrix
        if (isUp) {
            (*tr.mat)[J][I][Index][Index] = -flux.derivative(0);
            (*tr.mat)[I][I][Index][Index] += flux.derivative(0);
        }
    }

    template<class TrRe>
    void assembleTracerEquationFlux(TrRe& tr,
                                    const ElementContext& elemCtx,
                                    unsigned scvfIdx,
                                    unsigned I,
                                    unsigned J,
                                    const Scalar dt)
    {
        if (tr.numTracer() == 0) {
            return;
        }

        assembleTracerEquationFlux_<Free>(tr, elemCtx, scvfIdx, I, J, dt);
        assembleTracerEquationFlux_<Solution>(tr, elemCtx, scvfIdx, I, J, dt);
    }

    template<class TrRe, class Well>
    void assembleTracerEquationWell(TrRe& tr,
                                    const Well& well)
    {
        if (tr.numTracer() == 0) {
            return;
        }

        const auto& eclWell = well.wellEcl();

        // Init. well output to zero
        auto& tracerRate = this->wellTracerRate_[eclWell.seqIndex()];
        auto& solTracerRate = this->wellTracerRate_[eclWell.seqIndex()];
        auto& freeTracerRate = this->wellFreeTracerRate_[eclWell.seqIndex()];
        auto* mswTracerRate = eclWell.isMultiSegment() ? &this->mSwTracerRate_[eclWell.seqIndex()] : nullptr;
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            tracerRate.emplace_back(this->name(tr.idx_[tIdx]), 0.0);
            freeTracerRate.emplace_back(this->wellfname(tr.idx_[tIdx]), 0.0);
            solTracerRate.emplace_back(this->wellsname(tr.idx_[tIdx]), 0.0);
            if (eclWell.isMultiSegment()) {
                auto& wtr = mswTracerRate->emplace_back(this->name(tr.idx_[tIdx]));
                for (std::size_t i = 0; i < eclWell.getConnections().size(); ++i) {
                    wtr.rate[eclWell.getConnections().get(i).segment()] = 0.0;
                }
            }
        }

        std::vector<Scalar> wtracer(tr.numTracer());
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            wtracer[tIdx] = this->currentConcentration_(eclWell, this->name(tr.idx_[tIdx]));
        }

        const Scalar dt = simulator_.timeStepSize();
        const std::size_t well_index = simulator_.problem().wellModel().wellState().index(well.name()).value();
        const auto& ws = simulator_.problem().wellModel().wellState().well(well_index);
        for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
            const auto I = ws.perf_data.cell_index[i];
            const Scalar rate = well.volumetricSurfaceRateForConnection(I, tr.phaseIdx_);
            Scalar rate_s;
            if (tr.phaseIdx_ == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                rate_s = ws.perf_data.phase_mixing_rates[i][ws.vaporized_oil];
            }
            else if (tr.phaseIdx_ == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                rate_s = ws.perf_data.phase_mixing_rates[i][ws.dissolved_gas];
            }
            else {
                rate_s = 0.0;
            }

            const Scalar rate_f = rate - rate_s;
            if (rate_f > 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    // Injection of free tracer only
                    tr.residual_[tIdx][I][Free] -= rate_f*wtracer[tIdx];

                    // Store _injector_ tracer rate for reporting
                    // (can be done here since WTRACER is constant)
                    tracerRate[tIdx].rate += rate_f*wtracer[tIdx];
                    freeTracerRate[tIdx].rate += rate_f*wtracer[tIdx];
                    if (eclWell.isMultiSegment()) {
                        (*mswTracerRate)[tIdx].rate[eclWell.getConnections().get(i).segment()] += rate_f*wtracer[tIdx];
                    }
                }
                dVol_[tr.phaseIdx_][I][Free] -= rate_f * dt;
            }
            else if (rate_f < 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    // Store _injector_ tracer rate for cross-flowing well connections
                    // (can be done here since WTRACER is constant)
                    tracerRate[tIdx].rate += rate_f*wtracer[tIdx];
                    freeTracerRate[tIdx].rate += rate_f*wtracer[tIdx];

                    // Production of free tracer
                    tr.residual_[tIdx][I][Free] -= rate_f * tr.concentration_[tIdx][I][Free];
                }
                dVol_[tr.phaseIdx_][I][Free] -= rate_f * dt;

                // Derivative matrix for free tracer producer
                (*tr.mat)[I][I][Free][Free] -= rate_f * variable<TracerEvaluation>(1.0, 0).derivative(0);
            }
            if (rate_s < 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    // Production of solution tracer
                    tr.residual_[tIdx][I][Solution] -= rate_s * tr.concentration_[tIdx][I][Solution];
                }
                dVol_[tr.phaseIdx_][I][Solution] -= rate_s * dt;

                // Derivative matrix for solution tracer producer
                (*tr.mat)[I][I][Solution][Solution] -= rate_s * variable<TracerEvaluation>(1.0, 0).derivative(0);
            }
        }
    }

    template<class TrRe>
    void assembleTracerEquationSource(TrRe& tr,
                                      const Scalar dt,
                                      unsigned I)
    {
        if (tr.numTracer() == 0) {
            return;
        }

        // Skip if solution tracers do not exist
        if (tr.phaseIdx_ ==  FluidSystem::waterPhaseIdx ||
            (tr.phaseIdx_ ==  FluidSystem::gasPhaseIdx && !FluidSystem::enableDissolvedGas()) ||
            (tr.phaseIdx_ ==  FluidSystem::oilPhaseIdx && !FluidSystem::enableVaporizedOil()))
        {
            return;
        }

        const Scalar& dsVol = dVol_[tr.phaseIdx_][I][Solution];
        const Scalar& dfVol = dVol_[tr.phaseIdx_][I][Free];

        // Source term determined by sign of dsVol: if dsVol > 0 then ms -> mf, else mf -> ms
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            if (dsVol >= 0) {
                const auto delta = (dfVol / dt) * tr.concentration_[tIdx][I][Free];
                tr.residual_[tIdx][I][Free] -= delta;
                tr.residual_[tIdx][I][Solution] += delta;
            }
            else {
                const auto delta = (dsVol / dt) * tr.concentration_[tIdx][I][Solution];
                tr.residual_[tIdx][I][Free] += delta;
                tr.residual_[tIdx][I][Solution] -= delta;
            }
        }

        // Derivative matrix
        if (dsVol >= 0) {
            const auto delta = (dfVol / dt) * variable<TracerEvaluation>(1.0, 0).derivative(0);
            (*tr.mat)[I][I][Free][Free] -= delta;
            (*tr.mat)[I][I][Solution][Free] += delta;
        }
        else {
            const auto delta = (dsVol / dt) * variable<TracerEvaluation>(1.0, 0).derivative(0);
            (*tr.mat)[I][I][Free][Solution] += delta;
            (*tr.mat)[I][I][Solution][Solution] -= delta;
        }
    }

    void assembleTracerEquations_()
    {
        // Note that we formulate the equations in terms of a concentration update
        // (compared to previous time step) and not absolute concentration.
        // This implies that current concentration (tr.concentration_[][]) contributes
        // to the rhs both through storage and flux terms.
        // Compare also advanceTracerFields(...) below.

        OPM_TIMEBLOCK(tracerAssemble);
        for (auto& tr : tbatch) {
            if (tr.numTracer() != 0) {
                (*tr.mat) = 0.0;
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    tr.residual_[tIdx] = 0.0;
                }
            }
        }

        this->wellTracerRate_.clear();
        this->wellFreeTracerRate_.clear();
        this->wellSolTracerRate_.clear();

        // Well terms
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        for (const auto& wellPtr : wellPtrs) {
            for (auto& tr : tbatch) {
                this->assembleTracerEquationWell(tr, *wellPtr);
            }
        }

        ElementContext elemCtx(simulator_);
        const Scalar dt = elemCtx.simulator().timeStepSize();
        for (const auto& elem : elements(simulator_.gridView())) {
            elemCtx.updateStencil(elem);

            const std::size_t I = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/0);

            if (elem.partitionType() != Dune::InteriorEntity) {
                // Dirichlet boundary conditions needed for the parallel matrix
                for (const auto& tr : tbatch) {
                    if (tr.numTracer() != 0) {
                        (*tr.mat)[I][I][Free][Free] = 1.;
                        (*tr.mat)[I][I][Solution][Solution] = 1.;
                    }
                }
                continue;
            }
            elemCtx.updateAllIntensiveQuantities();
            elemCtx.updateAllExtensiveQuantities();

            const Scalar extrusionFactor =
                    elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
            Valgrind::CheckDefined(extrusionFactor);
            assert(isfinite(extrusionFactor));
            assert(extrusionFactor > 0.0);
            const Scalar scvVolume =
                    elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume()
                    * extrusionFactor;

            const std::size_t I1 = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/1);

            for (auto& tr : tbatch) {
                this->assembleTracerEquationVolume(tr, elemCtx, scvVolume, dt, I, I1);
            }

            const std::size_t numInteriorFaces = elemCtx.numInteriorFaces(/*timeIdx=*/0);
            for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
                const auto& face = elemCtx.stencil(0).interiorFace(scvfIdx);
                const unsigned j = face.exteriorIndex();
                const unsigned J = elemCtx.globalSpaceIndex(/*dofIdx=*/ j, /*timeIdx=*/0);
                for (auto& tr : tbatch) {
                    this->assembleTracerEquationFlux(tr, elemCtx, scvfIdx, I, J, dt);
                }
            }

            // Source terms (mass transfer between free and solution tracer)
            for (auto& tr : tbatch) {
                this->assembleTracerEquationSource(tr, dt, I);
            }
        }

        // Communicate overlap using grid Communication
        for (auto& tr : tbatch) {
            if (tr.numTracer() == 0) {
                continue;
            }
            auto handle = VectorVectorDataHandle<GridView, std::vector<TracerVector>>(tr.residual_,
                                                                                      simulator_.gridView());
            simulator_.gridView().communicate(handle, Dune::InteriorBorder_All_Interface,
                                              Dune::ForwardCommunication);
        }
    }

    template<TracerTypeIdx Index, class TrRe>
    void updateElem(TrRe& tr,
                    const Scalar scvVolume,
                    const unsigned globalDofIdx)
    {
        const Scalar vol1 = computeVolume_<Index>(tr.phaseIdx_, globalDofIdx, 0);
        vol1_[tr.phaseIdx_][globalDofIdx][Index] = vol1 * scvVolume;
        dVol_[tr.phaseIdx_][globalDofIdx][Index] = 0.0;
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            tr.storageOfTimeIndex1_[tIdx][globalDofIdx][Index] =
                vol1 * tr.concentrationInitial_[tIdx][globalDofIdx][Index];
        }
    }

    void updateStorageCache()
    {
        for (auto& tr : tbatch) {
            if (tr.numTracer() != 0) {
                tr.concentrationInitial_ = tr.concentration_;
            }
        }

        ElementContext elemCtx(simulator_);
        for (const auto& elem : elements(simulator_.gridView())) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const Scalar extrusionFactor = elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
            const Scalar scvVolume = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume() * extrusionFactor;
            const unsigned globalDofIdx = elemCtx.globalSpaceIndex(0, /*timeIdx=*/0);
            for (auto& tr : tbatch) {
                if (tr.numTracer() == 0) {
                    continue;
                }
                updateElem<Free>(tr, scvVolume, globalDofIdx);
                updateElem<Solution>(tr, scvVolume, globalDofIdx);
            }
        }
    }

    void advanceTracerFields()
    {
        assembleTracerEquations_();

        for (auto& tr : tbatch) {
            if (tr.numTracer() == 0) {
                continue;
            }

            // Note that we solve for a concentration update (compared to previous time step)
            // Confer also assembleTracerEquations_(...) above.
            std::vector<TracerVector> dx(tr.concentration_);
            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                dx[tIdx] = 0.0;
            }

            const bool converged = this->linearSolveBatchwise_(*tr.mat, dx, tr.residual_);
            if (!converged) {
                OpmLog::warning("### Tracer model: Linear solver did not converge. ###");
            }

            OPM_TIMEBLOCK(tracerPost);

            auto limit = [&tr, &dx,
                          &splitConcentration = this->splitTracerConcentration_]
                          (TracerTypeIdx index,
                           const Scalar S,
                           const unsigned tIdx,
                           const unsigned globalDofIdx)
            {
                constexpr Scalar tol_gas_sat = 1e-6;
                if (tr.concentration_[tIdx][globalDofIdx][index] - dx[tIdx][globalDofIdx][index] < 0.0 ||
                    S < tol_gas_sat)
                {
                    tr.concentration_[tIdx][globalDofIdx][index] = 0.0;
                }
                else {
                    tr.concentration_[tIdx][globalDofIdx][index] -= dx[tIdx][globalDofIdx][index];
                }
                splitConcentration[index][tr.idx_[tIdx]][globalDofIdx] =
                    tr.concentration_[tIdx][globalDofIdx][index];
            };

            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                for (std::size_t globalDofIdx = 0; globalDofIdx < tr.concentration_[tIdx].size(); ++globalDofIdx) {
                    // New concetration. Concentrations that are negative or where free/solution phase is not
                    // present are set to zero
                    const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, 0);
                    const auto& fs = intQuants.fluidState();
                    const Scalar Sf = decay<Scalar>(fs.saturation(tr.phaseIdx_));
                    Scalar Ss = 0.0;

                    if (tr.phaseIdx_ == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                        Ss = decay<Scalar>(fs.saturation(FluidSystem::oilPhaseIdx));
                    }
                    else if (tr.phaseIdx_ == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                        Ss = decay<Scalar>(fs.saturation(FluidSystem::gasPhaseIdx));
                    }

                    limit(Free, Sf, tIdx, globalDofIdx);
                    limit(Solution, Ss, tIdx, globalDofIdx);
                }
            }

            // Store _producer_ tracer rate for reporting
            const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
            for (const auto& wellPtr : wellPtrs) {
                const auto& eclWell = wellPtr->wellEcl();

                // Injection rates already reported during assembly
                if (!eclWell.isProducer()) {
                    continue;
                }

                Scalar rateWellPos = 0.0;
                Scalar rateWellNeg = 0.0;
                const std::size_t well_index = simulator_.problem().wellModel().wellState().index(eclWell.name()).value();
                const auto& ws = simulator_.problem().wellModel().wellState().well(well_index);
                auto& tracerRate = this->wellTracerRate_[well_index];
                auto& freeTracerRate = this->wellFreeTracerRate_[well_index];
                auto& solTracerRate = this->wellSolTracerRate_[well_index];
                auto* mswTracerRate = eclWell.isMultiSegment() ? &this->mSwTracerRate_[well_index] : nullptr;

                auto assign = [&tr, &eclWell,
                               &tracerRate, &mswTracerRate](const TracerTypeIdx index,
                                                            const std::size_t i,
                                                            const unsigned I,
                                                            const Scalar rate,
                                                            std::vector<TracerRate<Scalar>>& splitRate)
                {
                    if (rate < 0) {
                        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                            // Store _producer_ free tracer rate for reporting
                            tracerRate[tIdx].rate += rate * tr.concentration_[tIdx][I][index];
                            splitRate[tIdx].rate += rate * tr.concentration_[tIdx][I][index];
                            if (eclWell.isMultiSegment()) {
                                (*mswTracerRate)[tIdx].rate[eclWell.getConnections().get(i).segment()] +=
                                    rate * tr.concentration_[tIdx][I][index];
                            }
                        }
                    }
                };

                for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
                    const auto I = ws.perf_data.cell_index[i];
                    const Scalar rate = wellPtr->volumetricSurfaceRateForConnection(I, tr.phaseIdx_);

                    Scalar rate_s;
                    if (tr.phaseIdx_ == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                        rate_s = ws.perf_data.phase_mixing_rates[i][ws.vaporized_oil];
                    }
                    else if (tr.phaseIdx_ == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                        rate_s = ws.perf_data.phase_mixing_rates[i][ws.dissolved_gas];
                    }
                    else {
                        rate_s = 0.0;
                    }

                    const Scalar rate_f = rate - rate_s;
                    assign(Free, i, I, rate_f, freeTracerRate);
                    assign(Solution, i, I, rate_s, solTracerRate);

                    if (rate < 0) {
                        rateWellNeg += rate;
                    }
                    else {
                        rateWellPos += rate;
                    }
                }

                //Scalar rateWellTotal = rateWellNeg + rateWellPos;

                // TODO: Some inconsistencies here that perhaps should be clarified.
                // The "offical" rate as reported below is occasionally significant
                // different from the sum over connections (as calculated above). Only observed
                // for small values, neglible for the rate itself, but matters when used to
                // calculate tracer concentrations.
                const Scalar official_well_rate_total =
                    simulator_.problem().wellModel().wellState().well(well_index).surface_rates[tr.phaseIdx_];

                const Scalar rateWellTotal = official_well_rate_total;

                if (rateWellTotal > rateWellNeg) { // Cross flow
                    const Scalar bucketPrDay = 10.0/(1000.*3600.*24.); // ... keeps (some) trouble away
                    const Scalar factor = (rateWellTotal < -bucketPrDay) ? rateWellTotal/rateWellNeg : 0.0;
                    for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                        tracerRate[tIdx].rate *= factor;
                    }
                }
            }
        }
    }

    Simulator& simulator_;

    // This struct collects tracers of the same type (i.e, transported in same phase).
    // The idea being that, under the assumption of linearity, tracers of same type can
    // be solved in concert, having a common system matrix but separate right-hand-sides.

    // Since oil or gas tracers appears in dual compositions when VAPOIL respectively DISGAS
    // is active, the template argument is intended to support future extension to these
    // scenarios by supplying an extended vector type.

    template <typename TV>
    struct TracerBatch
    {
        std::vector<int> idx_;
        const int phaseIdx_;
        std::vector<TV> concentrationInitial_;
        std::vector<TV> concentration_;
        std::vector<TV> storageOfTimeIndex1_;
        std::vector<TV> residual_;
        std::unique_ptr<TracerMatrix> mat;

        bool operator==(const TracerBatch& rhs) const
        {
            return this->concentrationInitial_ == rhs.concentrationInitial_ &&
                   this->concentration_ == rhs.concentration_;
        }

        static TracerBatch serializationTestObject()
        {
            TracerBatch<TV> result(4);
            result.idx_ = {1,2,3};
            result.concentrationInitial_ = {5.0, 6.0};
            result.concentration_ = {7.0, 8.0};
            result.storageOfTimeIndex1_ = {9.0, 10.0, 11.0};
            result.residual_ = {12.0, 13.0};

            return result;
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(concentrationInitial_);
            serializer(concentration_);
        }

        TracerBatch(int phaseIdx = 0) : phaseIdx_(phaseIdx) {}

        int numTracer() const
        { return idx_.size(); }

        void addTracer(const int idx, const TV & concentration)
        {
            const int numGridDof = concentration.size();
            idx_.emplace_back(idx);
            concentrationInitial_.emplace_back(concentration);
            concentration_.emplace_back(concentration);
            residual_.emplace_back(numGridDof);
            storageOfTimeIndex1_.emplace_back(numGridDof);
        }
    };

    std::array<TracerBatch<TracerVector>,3> tbatch;
    TracerBatch<TracerVector>& wat_;
    TracerBatch<TracerVector>& oil_;
    TracerBatch<TracerVector>& gas_;
    std::array<std::array<std::vector<Scalar>, 3>, 2> vol1_;
    std::array<std::array<std::vector<Scalar>, 3>, 2> dVol_;
};

} // namespace Opm

#endif // OPM_TRACER_MODEL_HPP
