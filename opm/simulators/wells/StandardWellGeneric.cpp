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
#include <opm/simulators/wells/StandardWellGeneric.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>
#include <stdexcept>

namespace Opm
{

template<class Scalar>
StandardWellGeneric<Scalar>::
StandardWellGeneric(int Bhp,
                    const WellInterfaceGeneric& baseif)
    : baseif_(baseif)
    , perf_densities_(baseif_.numPerfs())
    , perf_pressure_diffs_(baseif_.numPerfs())
    , parallelB_(duneB_, baseif_.parallelWellInfo())
    , Bhp_(Bhp)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise);
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}


template<class Scalar>
double
StandardWellGeneric<Scalar>::
relaxationFactorRate(const std::vector<double>& primary_variables,
                     const BVectorWell& dwells)
{
    double relaxation_factor = 1.0;
    static constexpr int WQTotal = 0;

    // For injector, we only check the total rates to avoid sign change of rates
    const double original_total_rate = primary_variables[WQTotal];
    const double newton_update = dwells[0][WQTotal];
    const double possible_update_total_rate = primary_variables[WQTotal] - newton_update;

    // 0.8 here is a experimental value, which remains to be optimized
    // if the original rate is zero or possible_update_total_rate is zero, relaxation_factor will
    // always be 1.0, more thoughts might be needed.
    if (original_total_rate * possible_update_total_rate < 0.) { // sign changed
        relaxation_factor = std::abs(original_total_rate / newton_update) * 0.8;
    }

    assert(relaxation_factor >= 0.0 && relaxation_factor <= 1.0);

    return relaxation_factor;
}

template<class Scalar>
double
StandardWellGeneric<Scalar>::
relaxationFactorFraction(const double old_value,
                         const double dx)
{
    assert(old_value >= 0. && old_value <= 1.0);

    double relaxation_factor = 1.;

    // updated values without relaxation factor
    const double possible_updated_value = old_value - dx;

    // 0.95 is an experimental value remains to be optimized
    if (possible_updated_value < 0.0) {
        relaxation_factor = std::abs(old_value / dx) * 0.95;
    } else if (possible_updated_value > 1.0) {
        relaxation_factor = std::abs((1. - old_value) / dx) * 0.95;
    }
    // if possible_updated_value is between 0. and 1.0, then relaxation_factor
    // remains to be one

    assert(relaxation_factor >= 0. && relaxation_factor <= 1.);

    return relaxation_factor;
}

template<class Scalar>
void
StandardWellGeneric<Scalar>::
computeConnectionPressureDelta()
{
    // Algorithm:

    // We'll assume the perforations are given in order from top to
    // bottom for each well.  By top and bottom we do not necessarily
    // mean in a geometric sense (depth), but in a topological sense:
    // the 'top' perforation is nearest to the surface topologically.
    // Our goal is to compute a pressure delta for each perforation.

    // 1. Compute pressure differences between perforations.
    //    dp_perf will contain the pressure difference between a
    //    perforation and the one above it, except for the first
    //    perforation for each well, for which it will be the
    //    difference to the reference (bhp) depth.

    const int nperf = baseif_.numPerfs();
    perf_pressure_diffs_.resize(nperf, 0.0);
    auto z_above = baseif_.parallelWellInfo().communicateAboveValues(baseif_.refDepth(), baseif_.perfDepth());

    for (int perf = 0; perf < nperf; ++perf) {
        const double dz = baseif_.perfDepth()[perf] - z_above[perf];
        perf_pressure_diffs_[perf] = dz * perf_densities_[perf] * baseif_.gravity();
    }

    // 2. Compute pressure differences to the reference point (bhp) by
    //    accumulating the already computed adjacent pressure
    //    differences, storing the result in dp_perf.
    //    This accumulation must be done per well.
    const auto beg = perf_pressure_diffs_.begin();
    const auto end = perf_pressure_diffs_.end();
    baseif_.parallelWellInfo().partialSumPerfValues(beg, end);
}

template<class Scalar>
std::optional<double>
StandardWellGeneric<Scalar>::
computeBhpAtThpLimitProdWithAlq(const std::function<std::vector<double>(const double)>& frates,
                                const SummaryState& summary_state,
                                DeferredLogger& deferred_logger,
                                double maxPerfPress,
                                double alq_value) const
{
    return WellBhpThpCalculator(baseif_).computeBhpAtThpLimitProd(frates,
                                                                  summary_state,
                                                                  maxPerfPress,
                                                                  this->getRho(),
                                                                  alq_value,
                                                                  deferred_logger);
}

template<class Scalar>
void
StandardWellGeneric<Scalar>::
checkConvergenceControlEq(const WellState& well_state,
                          ConvergenceReport& report,
                          DeferredLogger& deferred_logger,
                          const double max_residual_allowed) const
{
    double control_tolerance = 0.;
    using CR = ConvergenceReport;
    CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;

    const int well_index = baseif_.indexOfWell();
    const auto& ws = well_state.well(well_index);
    if (baseif_.wellIsStopped()) {
        ctrltype = CR::WellFailure::Type::ControlRate;
        control_tolerance = 1.e-6; // use smaller tolerance for zero control?
    }
    else if (baseif_.isInjector() )
    {
        auto current = ws.injection_cmode;
        switch(current) {
        case Well::InjectorCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = 1.e4; // 0.1 bar
            break;
        case Well::InjectorCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = 1.e3; // 0.01 bar
            break;
        case Well::InjectorCMode::RATE:
        case Well::InjectorCMode::RESV:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-4; //
            break;
        case Well::InjectorCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-6; //
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
        }
    }
    else if (baseif_.isProducer() )
    {
        auto current = ws.production_cmode;
        switch(current) {
        case Well::ProducerCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = 1.e4; // 0.1 bar
            break;
        case Well::ProducerCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = 1.e3; // 0.01 bar
            break;
        case Well::ProducerCMode::ORAT:
        case Well::ProducerCMode::WRAT:
        case Well::ProducerCMode::GRAT:
        case Well::ProducerCMode::LRAT:
        case Well::ProducerCMode::RESV:
        case Well::ProducerCMode::CRAT:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-4; // smaller tolerance for rate control
            break;
        case Well::ProducerCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-6; // smaller tolerance for rate control
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
        }
    }

    const double well_control_residual = std::abs(this->resWell_[0][Bhp_]);
    const int dummy_component = -1;
    if (std::isnan(well_control_residual)) {
        report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, baseif_.name()});
    } else if (well_control_residual > max_residual_allowed * 10.) {
        report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, baseif_.name()});
    } else if ( well_control_residual > control_tolerance) {
        report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, baseif_.name()});
    }
}

template<class Scalar>
void
StandardWellGeneric<Scalar>::
checkConvergencePolyMW(const std::vector<double>& res,
                       ConvergenceReport& report,
                       const double maxResidualAllowed) const
{
  if (baseif_.isInjector()) {
      //  checking the convergence of the perforation rates
      const double wat_vel_tol = 1.e-8;
      const int dummy_component = -1;
      using CR = ConvergenceReport;
      const auto wat_vel_failure_type = CR::WellFailure::Type::MassBalance;
      for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
          const double wat_vel_residual = res[Bhp_ + 1 + perf];
          if (std::isnan(wat_vel_residual)) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::NotANumber, dummy_component, baseif_.name()});
          } else if (wat_vel_residual > maxResidualAllowed * 10.) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::TooLarge, dummy_component, baseif_.name()});
          } else if (wat_vel_residual > wat_vel_tol) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::Normal, dummy_component, baseif_.name()});
          }
      }

      // checking the convergence of the skin pressure
      const double pskin_tol = 1000.; // 1000 pascal
      const auto pskin_failure_type = CR::WellFailure::Type::Pressure;
      for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
          const double pskin_residual = res[Bhp_ + 1 + perf + baseif_.numPerfs()];
          if (std::isnan(pskin_residual)) {
              report.setWellFailed({pskin_failure_type, CR::Severity::NotANumber, dummy_component, baseif_.name()});
          } else if (pskin_residual > maxResidualAllowed * 10.) {
              report.setWellFailed({pskin_failure_type, CR::Severity::TooLarge, dummy_component, baseif_.name()});
          } else if (pskin_residual > pskin_tol) {
              report.setWellFailed({pskin_failure_type, CR::Severity::Normal, dummy_component, baseif_.name()});
          }
      }
  }
}


template<class Scalar>
void
StandardWellGeneric<Scalar>::
getNumBlocks(unsigned int& numBlocks) const
{
    numBlocks = duneB_.nonzeroes();
}

template class StandardWellGeneric<double>;

}
