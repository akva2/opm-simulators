/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/wells/MultisegmentWellAssemble.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/material/densead/DynamicEvaluation.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/linalg/matrixblock.hh>

#include <opm/simulators/wells/MultisegmentWellEquations.hpp>
#include <opm/simulators/wells/WellAssemble.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

template<typename FluidSystem, typename Indices, typename Scalar>
template<class EvalWell>
void
MultisegmentWellAssemble<FluidSystem,Indices,Scalar>::
assembleControlEq(const WellState& well_state,
                  const GroupState& group_state,
                  const Schedule& schedule,
                  const SummaryState& summaryState,
                  const Well::InjectionControls& inj_controls,
                  const Well::ProductionControls& prod_controls,
                  const EvalWell& wqTotal,
                  const EvalWell& bhp,
                  const double rho,
                  const int SPres,
                  const std::function<EvalWell(int)>& getQs,
                  MultisegmentWellEquations<Indices,Scalar>& eqns,
                  DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    EvalWell control_eq(0.0);

    const auto& well = well_.wellEcl();

    auto getRates = [&]() {
        std::vector<EvalWell> rates(3, EvalWell(eqns.numWellEq + Indices::numEq, 0.0));
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            rates[Water] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            rates[Oil] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            rates[Gas] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx));
        }
        return rates;
    };

    if (well_.wellIsStopped()) {
        control_eq = wqTotal;
    } else if (well_.isInjector() ) {
        // Find scaling factor to get injection rate,
        const InjectorType injectorType = inj_controls.injector_type;
        double scaling = 1.0;
        const auto& pu = well_.phaseUsage();
        switch (injectorType) {
        case InjectorType::WATER:
        {
            scaling = well_.scalingFactor(pu.phase_pos[BlackoilPhases::Aqua]);
            break;
        }
        case InjectorType::OIL:
        {
            scaling = well_.scalingFactor(pu.phase_pos[BlackoilPhases::Liquid]);
            break;
        }
        case InjectorType::GAS:
        {
            scaling = well_.scalingFactor(pu.phase_pos[BlackoilPhases::Vapour]);
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injectors " + well.name());
        }
        const EvalWell injection_rate = wqTotal / scaling;
        // Setup function for evaluation of BHP from THP (used only if needed).
        std::function<EvalWell()> bhp_from_thp = [&]() {
            const auto rates = getRates();
            return WellBhpThpCalculator(well_).calculateBhpFromThp(well_state,
                                                                   rates,
                                                                   well,
                                                                   summaryState,
                                                                   rho,
                                                                   deferred_logger);
        };
        // Call generic implementation.
        WellAssemble(well_).assembleControlEqInj(well_state,
                                                 group_state,
                                                 schedule,
                                                 summaryState,
                                                 inj_controls,
                                                 bhp,
                                                 injection_rate,
                                                 bhp_from_thp,
                                                 control_eq,
                                                 deferred_logger);
    } else {
        // Find rates.
        const auto rates = getRates();
        // Setup function for evaluation of BHP from THP (used only if needed).
        std::function<EvalWell()> bhp_from_thp = [&]() {
            return WellBhpThpCalculator(well_).calculateBhpFromThp(well_state,
                                                                   rates,
                                                                   well,
                                                                   summaryState,
                                                                   rho,
                                                                   deferred_logger);
        };
        // Call generic implementation.
        WellAssemble(well_).assembleControlEqProd(well_state,
                                                  group_state,
                                                  schedule,
                                                  summaryState,
                                                  prod_controls,
                                                  bhp,
                                                  rates,
                                                  bhp_from_thp,
                                                  control_eq,
                                                  deferred_logger);
    }

    // using control_eq to update the matrix and residuals
    eqns.resWell_[0][SPres] = control_eq.value();
    for (int pv_idx = 0; pv_idx < eqns.numWellEq; ++pv_idx) {
        eqns.duneD_[0][0][SPres][pv_idx] = control_eq.derivative(pv_idx + Indices::numEq);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
template<class EvalWell>
void MultisegmentWellAssemble<FluidSystem,Indices,Scalar>::
assemblePressureEq(const int seg,
                   const int seg_upwind,
                   const EvalWell& pressure_equation,
                   const EvalWell& outlet_pressure,
                   const int WFrac,
                   const int GFrac,
                   const int SPres,
                   const int WQTotal,
                   const int outlet_segment_index,
                   MultisegmentWellEquations<Indices,Scalar>& eqns) const
{
    eqns.resWell_[seg][SPres] = pressure_equation.value();
    eqns.duneD_[seg][seg][SPres][SPres] += pressure_equation.derivative(SPres + Indices::numEq);
    eqns.duneD_[seg][seg][SPres][WQTotal] += pressure_equation.derivative(WQTotal + Indices::numEq);
    if (WFrac > -1) {
        eqns.duneD_[seg][seg_upwind][SPres][WFrac] += pressure_equation.derivative(WFrac + Indices::numEq);
    }
    if (GFrac > -1) {
        eqns.duneD_[seg][seg_upwind][SPres][GFrac] += pressure_equation.derivative(GFrac + Indices::numEq);
    }

    // contribution from the outlet segment
    eqns.resWell_[seg][SPres] -= outlet_pressure.value();
    for (int pv_idx = 0; pv_idx < eqns.numWellEq; ++pv_idx) {
        eqns.duneD_[seg][outlet_segment_index][SPres][pv_idx] = -outlet_pressure.derivative(pv_idx + Indices::numEq);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void MultisegmentWellAssemble<FluidSystem,Indices,Scalar>::
assembleTrivialEq(const int seg,
                  Scalar value,
                  const int SPres,
                  const int WQTotal,
                  MultisegmentWellEquations<Indices,Scalar>& eqns) const
{
    eqns.resWell_[seg][SPres] = value;
    eqns.duneD_[seg][seg][SPres][WQTotal] = 1.;
}

template<typename FluidSystem, typename Indices, typename Scalar>
template<class EvalWell>
void MultisegmentWellAssemble<FluidSystem,Indices,Scalar>::
assembleAccelerationTerm(const int seg,
                         const int comp_idx,
                         const EvalWell& accumulation_term,
                         MultisegmentWellEquations<Indices,Scalar>& eqns) const
{
    eqns.resWell_[seg][comp_idx] += accumulation_term.value();
    for (int pv_idx = 0; pv_idx < eqns.numWellEq; ++pv_idx) {
        eqns.duneD_[seg][seg][comp_idx][pv_idx] += accumulation_term.derivative(pv_idx + Indices::numEq);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
template<class EvalWell>
void MultisegmentWellAssemble<FluidSystem,Indices,Scalar>::
assembleFlowTerm(const int seg,
                 const int seg_upwind,
                 const int comp_idx,
                 const int WFrac,
                 const int GFrac,
                 const int WQTotal,
                 const EvalWell& segment_rate,
                 MultisegmentWellEquations<Indices,Scalar>& eqns) const
{
    // pressure derivative should be zero

    eqns.resWell_[seg][comp_idx] -= segment_rate.value();
    eqns.duneD_[seg][seg][comp_idx][WQTotal] -= segment_rate.derivative(WQTotal + Indices::numEq);
    if (WFrac > -1) {
        eqns.duneD_[seg][seg_upwind][comp_idx][WFrac] -= segment_rate.derivative(WFrac + Indices::numEq);
    }
    if (GFrac > -1) {
        eqns.duneD_[seg][seg_upwind][comp_idx][GFrac] -= segment_rate.derivative(GFrac + Indices::numEq);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
template<class EvalWell>
void MultisegmentWellAssemble<FluidSystem,Indices,Scalar>::
assemblePerfTerm(const int seg,
                 const int cell_idx,
                 const int comp_idx,
                 const EvalWell& cq_s_effective,
                 MultisegmentWellEquations<Indices,Scalar>& eqns) const
{
    // subtract sum of phase fluxes in the well equations.
    eqns.resWell_[seg][comp_idx] += cq_s_effective.value();

         // assemble the jacobians
    for (int pv_idx = 0; pv_idx < eqns.numWellEq; ++pv_idx) {

             // also need to consider the efficiency factor when manipulating the jacobians.
        eqns.duneC_[seg][cell_idx][pv_idx][comp_idx] -= cq_s_effective.derivative(pv_idx + Indices::numEq); // intput in transformed matrix

                  // the index name for the D should be eq_idx / pv_idx
        eqns.duneD_[seg][seg][comp_idx][pv_idx] += cq_s_effective.derivative(pv_idx + Indices::numEq);
    }

    for (int pv_idx = 0; pv_idx < Indices::numEq; ++pv_idx) {
        // also need to consider the efficiency factor when manipulating the jacobians.
        eqns.duneB_[seg][cell_idx][comp_idx][pv_idx] += cq_s_effective.derivative(pv_idx);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
template<class PressureMatrix, class BVector>
void MultisegmentWellAssemble<FluidSystem,Indices,Scalar>::
addWellPressureEqs(const BVector& weights,
                   const int pressureVarIndex,
                   const int seg_pressure_var_ind,
                   const WellState& well_state,
                   const MultisegmentWellEquations<Indices,Scalar>& eqns,
                   PressureMatrix& jacobian) const
{
    // Add the pressure contribution to the cpr system for the well

    // Add for coupling from well to reservoir
    const int welldof_ind = eqns.duneC_.M() + well_.indexOfWell();
    if (!well_.isPressureControlled(well_state)) {
        for (size_t rowC = 0; rowC < eqns.duneC_.N(); ++rowC) {
            for (auto colC = eqns.duneC_[rowC].begin(),
                      endC = eqns.duneC_[rowC].end(); colC != endC; ++colC) {
                const auto row_index = colC.index();
                const auto& bw = weights[row_index];
                double matel = 0.0;

                for(size_t i = 0; i< bw.size(); ++i){
                    matel += bw[i]*(*colC)[seg_pressure_var_ind][i];
                }
                jacobian[row_index][welldof_ind] += matel;
            }
        }
    }
    // make cpr weights for well by pure avarage of reservoir weights of the perforations
    if (!well_.isPressureControlled(well_state)) {
        auto well_weight = weights[0];
        well_weight = 0.0;
        int num_perfs = 0;
        for (size_t rowB = 0; rowB < eqns.duneB_.N(); ++rowB) {
            for (auto colB = eqns.duneB_[rowB].begin(),
                      endB = eqns.duneB_[rowB].end(); colB != endB; ++colB) {
                const auto col_index = colB.index();
                const auto& bw = weights[col_index];
                well_weight += bw;
                num_perfs += 1;
            }
        }

        well_weight /= num_perfs;
        assert(num_perfs>0);

        // Add for coupling from reservoir to well and caclulate diag elelement corresping to incompressible standard well
        double diag_ell = 0.0;
        for (size_t rowB = 0; rowB < eqns.duneB_.N(); ++rowB) {
            const auto& bw = well_weight;
            for (auto colB = eqns.duneB_[rowB].begin(),
                      endB = eqns.duneB_[rowB].end(); colB != endB; ++colB) {
                const auto col_index = colB.index();
                double matel = 0.0;
                for(size_t i = 0; i< bw.size(); ++i){
                    matel += bw[i] *(*colB)[i][pressureVarIndex];
                }
                jacobian[welldof_ind][col_index] += matel;
                diag_ell -= matel;
            }
        }

        assert(diag_ell > 0.0);
        jacobian[welldof_ind][welldof_ind] = diag_ell;
    } else {
        jacobian[welldof_ind][welldof_ind] = 1.0; // maybe we could have used diag_ell if calculated
    }
}

#define INSTANCE(Dim,Block,...) \
template class MultisegmentWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>; \
template void \
MultisegmentWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assembleControlEq(const WellState&, \
                  const GroupState&, \
                  const Schedule&, \
                  const SummaryState&, \
                  const Well::InjectionControls&, \
                  const Well::ProductionControls&, \
                  const DenseAd::Evaluation<double,Dim,0u>&, \
                  const DenseAd::Evaluation<double,Dim,0u>&, \
                  const double, \
                  const int, \
                  const std::function<DenseAd::Evaluation<double,Dim,0u>(int)>&, \
                  MultisegmentWellEquations<__VA_ARGS__,double>&, \
                  DeferredLogger&) const; \
template void \
MultisegmentWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assemblePressureEq(const int, \
                   const int, \
                   const DenseAd::Evaluation<double,Dim,0u>&, \
                   const DenseAd::Evaluation<double,Dim,0u>&, \
                   const int, \
                   const int, \
                   const int, \
                   const int, \
                   const int, \
                   MultisegmentWellEquations<__VA_ARGS__,double>&) const; \
template void \
MultisegmentWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assembleAccelerationTerm(const int, \
                         const int, \
                         const DenseAd::Evaluation<double,Dim,0u>&, \
                         MultisegmentWellEquations<__VA_ARGS__,double>&) const; \
template void \
MultisegmentWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assembleFlowTerm(const int, \
                 const int, \
                 const int, \
                 const int, \
                 const int, \
                 const int, \
                 const DenseAd::Evaluation<double,Dim,0u>&, \
                 MultisegmentWellEquations<__VA_ARGS__,double>&) const; \
template void \
MultisegmentWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assemblePerfTerm(const int, \
                 const int, \
                 const int, \
                 const DenseAd::Evaluation<double,Dim,0u>&, \
                 MultisegmentWellEquations<__VA_ARGS__,double>&) const; \
template void \
MultisegmentWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
addWellPressureEqs(const Dune::BlockVector<Dune::FieldVector<double,Block>>&, \
                   const int, \
                   const int, \
                   const WellState&, \
                   const MultisegmentWellEquations<__VA_ARGS__,double>&, \
                   Dune::BCRSMatrix<MatrixBlock<double,1,1>>& jacobian) const;

// One phase
INSTANCE(3, 1, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(4, 2, BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(8, 6, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(5, 2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(5, 2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(5, 2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(6, 3, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(6, 3, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(6, 3, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(6, 3, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(7, 4, BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)

// Blackoil
INSTANCE(7, 3, BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(7, 3, BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(8, 4, BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(8, 4, BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(8, 4, BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(8, 4, BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(8, 4, BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(8, 4, BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(8, 4, BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)
INSTANCE(9, 5, BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

}
