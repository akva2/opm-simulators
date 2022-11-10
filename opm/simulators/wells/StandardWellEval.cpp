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
#include <opm/simulators/wells/StandardWellEval.hpp>

#include <opm/material/densead/DynamicEvaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <cassert>
#include <cmath>



namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
StandardWellEval<FluidSystem,Indices,Scalar>::
StandardWellEval(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif,
                 const bool has_polymermw)
    : StandardWellConnections<Scalar>(baseif)
    , baseif_(baseif)
    , primary_variables_(baseif_, has_polymermw)
    , F0_(StandardWellEquations<Indices,Scalar>::numWellConservationEq)
    , linSys_(baseif_.parallelWellInfo())
{
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellEval<FluidSystem,Indices,Scalar>::EvalWell
StandardWellEval<FluidSystem,Indices,Scalar>::
extendEval(const Eval& in) const
{
    EvalWell out(numWellEq_ + Indices::numEq, in.value());
    for(int eqIdx = 0; eqIdx < Indices::numEq;++eqIdx) {
        out.setDerivative(eqIdx, in.derivative(eqIdx));
    }
    return out;
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updateWellStateFromPrimaryVariables(WellState& well_state,
                                    DeferredLogger& deferred_logger) const
{
    this->primary_variables_.copyToWellState(well_state, deferred_logger);

    WellBhpThpCalculator(baseif_).
            updateThp(this->getRho(),
                      [this,&well_state]() { return this->baseif_.getALQ(well_state); },
                      {FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                       FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx),
                       FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)},
                      well_state, deferred_logger);
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
computeAccumWell()
{
    for (size_t eq_idx = 0; eq_idx < F0_.size(); ++eq_idx) {
        F0_[eq_idx] = this->primary_variables_.wellSurfaceVolumeFraction(eq_idx, numWellEq_).value();
    }
}

template<class FluidSystem, class Indices, class Scalar>
ConvergenceReport
StandardWellEval<FluidSystem,Indices,Scalar>::
getWellConvergence(const WellState& well_state,
                   const std::vector<double>& B_avg,
                   const double maxResidualAllowed,
                   const double tol_wells,
                   const double relaxed_tolerance_flow,
                   const bool relax_tolerance,
                   std::vector<double>& res,
                   DeferredLogger& deferred_logger) const
{
    res.resize(numWellEq_);
    for (int eq_idx = 0; eq_idx < numWellEq_; ++eq_idx) {
        // magnitude of the residual matters
        res[eq_idx] = std::abs(linSys_.getResidual()[0][eq_idx]);
    }

    std::vector<double> well_flux_residual(baseif_.numComponents());

    // Finish computation
    for (int compIdx = 0; compIdx < baseif_.numComponents(); ++compIdx )
    {
        well_flux_residual[compIdx] = B_avg[compIdx] * res[compIdx];
    }

    ConvergenceReport report;
    using CR = ConvergenceReport;
    CR::WellFailure::Type type = CR::WellFailure::Type::MassBalance;
    // checking if any NaN or too large residuals found
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }

        const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
        const int compIdx = Indices::canonicalToActiveComponentIndex(canonicalCompIdx);

        if (std::isnan(well_flux_residual[compIdx])) {
            report.setWellFailed({type, CR::Severity::NotANumber, compIdx, baseif_.name()});
        } else if (well_flux_residual[compIdx] > maxResidualAllowed) {
            report.setWellFailed({type, CR::Severity::TooLarge, compIdx, baseif_.name()});
        } else if (!relax_tolerance && well_flux_residual[compIdx] > tol_wells) {
            report.setWellFailed({type, CR::Severity::Normal, compIdx, baseif_.name()});
        } else if (well_flux_residual[compIdx] > relaxed_tolerance_flow) {
            report.setWellFailed({type, CR::Severity::Normal, compIdx, baseif_.name()});
        }
    }

    WellConvergence(baseif_).
        checkConvergenceControlEq(well_state,
                                  {1.e3, 1.e4, 1.e-4, 1.e-6, maxResidualAllowed},
                                   std::abs(linSys_.getResidual()[0][Bhp]),
                                  report,
                                  deferred_logger);

    return report;
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
computeConnectionDensities(const std::vector<double>& perfComponentRates,
                           const std::vector<double>& b_perf,
                           const std::vector<double>& rsmax_perf,
                           const std::vector<double>& rvmax_perf,
                           const std::vector<double>& rvwmax_perf,
                           const std::vector<double>& surf_dens_perf,
                           DeferredLogger& deferred_logger)
{
    // Verify that we have consistent input.
    const int nperf = baseif_.numPerfs();
    const int num_comp = baseif_.numComponents();

    // 1. Compute the flow (in surface volume units for each
    //    component) exiting up the wellbore from each perforation,
    //    taking into account flow from lower in the well, and
    //    in/out-flow at each perforation.
    std::vector<double> q_out_perf((nperf)*num_comp, 0.0);

    // Step 1 depends on the order of the perforations. Hence we need to
    // do the modifications globally.
    // Create and get the global perforation information and do this sequentially
    // on each process

    const auto& factory = baseif_.parallelWellInfo().getGlobalPerfContainerFactory();
    auto global_q_out_perf = factory.createGlobal(q_out_perf, num_comp);
    auto global_perf_comp_rates = factory.createGlobal(perfComponentRates, num_comp);

    // TODO: investigate whether we should use the following techniques to calcuate the composition of flows in the wellbore
    // Iterate over well perforations from bottom to top.
    for (int perf = factory.numGlobalPerfs() - 1; perf >= 0; --perf) {
        for (int component = 0; component < num_comp; ++component) {
            auto index = perf * num_comp + component;
            if (perf == factory.numGlobalPerfs() - 1) {
                // This is the bottom perforation. No flow from below.
                global_q_out_perf[index] = 0.0;
            } else {
                // Set equal to flow from below.
                global_q_out_perf[index] = global_q_out_perf[index + num_comp];
            }
            // Subtract outflow through perforation.
            global_q_out_perf[index] -= global_perf_comp_rates[index];
        }
    }

    // Copy the data back to local view
    factory.copyGlobalToLocal(global_q_out_perf, q_out_perf, num_comp);

    // 2. Compute the component mix at each perforation as the
    //    absolute values of the surface rates divided by their sum.
    //    Then compute volume ratios (formation factors) for each perforation.
    //    Finally compute densities for the segments associated with each perforation.
    std::vector<double> mix(num_comp,0.0);
    std::vector<double> x(num_comp);
    std::vector<double> surf_dens(num_comp);

    for (int perf = 0; perf < nperf; ++perf) {
        // Find component mix.
        const double tot_surf_rate = std::accumulate(q_out_perf.begin() + num_comp*perf,
                                                     q_out_perf.begin() + num_comp*(perf+1), 0.0);
        if (tot_surf_rate != 0.0) {
            for (int component = 0; component < num_comp; ++component) {
                mix[component] = std::fabs(q_out_perf[perf*num_comp + component]/tot_surf_rate);
            }
        } else if (num_comp == 1) {
            mix[num_comp-1] = 1.0;
        } else {
            std::fill(mix.begin(), mix.end(), 0.0);
            // No flow => use well specified fractions for mix.
            if (baseif_.isInjector()) {
                switch (baseif_.wellEcl().injectorType()) {
                case InjectorType::WATER:
                    mix[FluidSystem::waterCompIdx] = 1.0;
                    break;
                case InjectorType::GAS:
                    mix[FluidSystem::gasCompIdx] = 1.0;
                    break;
                case InjectorType::OIL:
                    mix[FluidSystem::oilCompIdx] = 1.0;
                    break;
                case InjectorType::MULTI:
                    // Not supported.
                    // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                    //                         "Multi phase injectors are not supported, requested for well " + name());
                    break;
                }
            } else {
                assert(baseif_.isProducer());
                // For the frist perforation without flow we use the preferred phase to decide the mix initialization.
                if (perf == 0) { //
                    switch (baseif_.wellEcl().getPreferredPhase()) {
                    case Phase::OIL:
                        mix[FluidSystem::oilCompIdx] = 1.0;
                        break;
                    case Phase::GAS:
                        mix[FluidSystem::gasCompIdx] = 1.0;
                        break;
                    case Phase::WATER:
                        mix[FluidSystem::waterCompIdx] = 1.0;
                        break;
                    default:
                        // No others supported.
                        break;
                    }
                // For the rest of the perforation without flow we use mix from the above perforation.
                } else {
                    mix = x;
                }

            }
        }
        // Compute volume ratio.
        x = mix;

        // Subtract dissolved gas from oil phase and vapporized oil from gas phase and vaporized water from gas phase
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned oilpos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            double rs = 0.0;
            double rv = 0.0;
            if (!rsmax_perf.empty() && mix[oilpos] > 1e-12) {
                rs = std::min(mix[gaspos]/mix[oilpos], rsmax_perf[perf]);
            }
            if (!rvmax_perf.empty() && mix[gaspos] > 1e-12) {
                rv = std::min(mix[oilpos]/mix[gaspos], rvmax_perf[perf]);
            }
            const double d = 1.0 - rs*rv;
            if (d <= 0.0) {
                std::ostringstream sstr;
                sstr << "Problematic d value " << d << " obtained for well " << baseif_.name()
                     << " during ccomputeConnectionDensities with rs " << rs
                     << ", rv " << rv
                     << " obtaining d " << d
                     << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                     << " for this connection.";
                deferred_logger.debug(sstr.str());
            } else {
                if (rs > 0.0) {
                    // Subtract gas in oil from gas mixture
                    x[gaspos] = (mix[gaspos] - mix[oilpos]*rs)/d;
                }
                if (rv > 0.0) {
                    // Subtract oil in gas from oil mixture
                    x[oilpos] = (mix[oilpos] - mix[gaspos]*rv)/d;
                }
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                //matrix system: (mix[oilpos] = q_os, x[oilpos] = bo*q_or, etc...)
                //┌             ┐   ┌                ┐  ┌           ┐
                //│mix[oilpos]  │   | 1     Rv     0 |  |x[oilpos]  |
                //│mix[gaspos]  │ = │ Rs    1      0 │  │x[gaspos]  │
                //│mix[waterpos]│   │ 0     Rvw    1 │  │x[waterpos │
                //└             ┘   └                ┘  └           ┘
                const unsigned waterpos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                double rvw = 0.0;
                if (!rvwmax_perf.empty() && mix[gaspos] > 1e-12) {
                    rvw = std::min(mix[waterpos]/mix[gaspos], rvwmax_perf[perf]);
                }
                if (rvw > 0.0) {
                    // Subtract water in gas from water mixture
                    x[waterpos] = mix[waterpos] - x[gaspos] * rvw;
                }
            }
        } else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            //no oil
            const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned waterpos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            double rvw = 0.0;
            if (!rvwmax_perf.empty() && mix[gaspos] > 1e-12) {
                rvw = std::min(mix[waterpos]/mix[gaspos], rvwmax_perf[perf]);
            }
            if (rvw > 0.0) {
               // Subtract water in gas from water mixture
               x[waterpos] = mix[waterpos] - mix[gaspos] * rvw;
            }
        }

        double volrat = 0.0;
        for (int component = 0; component < num_comp; ++component) {
            volrat += x[component] / b_perf[perf*num_comp+ component];
        }
        for (int component = 0; component < num_comp; ++component) {
            surf_dens[component] = surf_dens_perf[perf*num_comp+ component];
        }

        // Compute segment density.
        this->perf_densities_[perf] = std::inner_product(surf_dens.begin(), surf_dens.end(), mix.begin(), 0.0) / volrat;
    }
}


template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
init(std::vector<double>& perf_depth,
     const std::vector<double>& depth_arg,
     const int num_cells,
     const bool has_polymermw)
{
    perf_depth.resize(baseif_.numPerfs(), 0.);
    for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
        const int cell_idx = baseif_.cells()[perf];
        perf_depth[perf] = depth_arg[cell_idx];
    }

    // counting/updating primary variable numbers
    if (has_polymermw) {
        if (baseif_.isInjector()) {
            // adding a primary variable for water perforation rate per connection
            numWellEq_ += baseif_.numPerfs();
            // adding a primary variable for skin pressure per connection
            numWellEq_ += baseif_.numPerfs();
        }
    }

    // with the updated numWellEq_, we can initialize the primary variables and matrices now
    primary_variables_.resize(numWellEq_);

    // setup sparsity pattern for the matrices
    this->linSys_.init(num_cells, this->numWellEq_, baseif_.numPerfs(), baseif_.cells());
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellEval<FluidSystem,Indices,Scalar>::
addWellContribution(WellContributions& wellContribs) const
{
    linSys_.addWellContribution(wellContribs);
}

#define INSTANCE(...) \
template class StandardWellEval<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)

// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

}
