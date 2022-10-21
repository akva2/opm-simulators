/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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
#include <opm/simulators/wells/WellConstraints.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

namespace Opm
{

Well::ProducerCMode WellConstraints::
activeProductionConstraint(const SingleWellState& ws,
                           const SummaryState& summaryState,
                           const RateConvFunc& calcReservoirVoidageRates,
                           bool& thp_limit_violated_but_not_switched,
                           DeferredLogger& deferred_logger) const
{
    const PhaseUsage& pu = well_.phaseUsage();
    const auto controls = well_.wellEcl().productionControls(summaryState);
    const auto currentControl = ws.production_cmode;

    if (controls.hasControl(Well::ProducerCMode::BHP) && currentControl != Well::ProducerCMode::BHP) {
        const double bhp_limit = controls.bhp_limit;
        double current_bhp = ws.bhp;
        if (bhp_limit > current_bhp)
            return Well::ProducerCMode::BHP;
    }

    if (controls.hasControl(Well::ProducerCMode::ORAT) && currentControl != Well::ProducerCMode::ORAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]];
        if (controls.oil_rate < current_rate)
            return Well::ProducerCMode::ORAT;
    }

    if (controls.hasControl(Well::ProducerCMode::WRAT) && currentControl != Well::ProducerCMode::WRAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];
        if (controls.water_rate < current_rate)
            return Well::ProducerCMode::WRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::GRAT) && currentControl != Well::ProducerCMode::GRAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Vapour]];
        if (controls.gas_rate < current_rate)
            return Well::ProducerCMode::GRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::LRAT) && currentControl != Well::ProducerCMode::LRAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]];
        current_rate -= ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];

        bool skip = false;
        if (controls.liquid_rate == controls.oil_rate) {
            const double current_water_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];
            if (std::abs(current_water_rate) < 1e-12) {
                skip = true;
                deferred_logger.debug("LRAT_ORAT_WELL", "Well " + well_.name() + " The LRAT target is equal the ORAT target and the water rate is zero, skip checking LRAT");
            }
        }
        if (!skip && controls.liquid_rate < current_rate)
            return Well::ProducerCMode::LRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::RESV) && currentControl != Well::ProducerCMode::RESV) {
        double current_rate = 0.0;
        if (pu.phase_used[BlackoilPhases::Aqua])
            current_rate -= ws.reservoir_rates[pu.phase_pos[BlackoilPhases::Aqua]];

        if (pu.phase_used[BlackoilPhases::Liquid])
            current_rate -= ws.reservoir_rates[pu.phase_pos[BlackoilPhases::Liquid]];

        if (pu.phase_used[BlackoilPhases::Vapour])
            current_rate -= ws.reservoir_rates[pu.phase_pos[BlackoilPhases::Vapour]];

        if (controls.prediction_mode && controls.resv_rate < current_rate)
            return Well::ProducerCMode::RESV;

        if (!controls.prediction_mode) {
            const int fipreg = 0; // not considering the region for now
            const int np = well_.numPhases();

            std::vector<double> surface_rates(np, 0.0);
            if (pu.phase_used[BlackoilPhases::Aqua])
                surface_rates[pu.phase_pos[BlackoilPhases::Aqua]] = controls.water_rate;
            if (pu.phase_used[BlackoilPhases::Liquid])
                surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = controls.oil_rate;
            if (pu.phase_used[BlackoilPhases::Vapour])
                surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = controls.gas_rate;

            std::vector<double> voidage_rates(np, 0.0);
            calcReservoirVoidageRates(fipreg, well_.pvtRegionIdx(), surface_rates, voidage_rates);

            double resv_rate = 0.0;
            for (int p = 0; p < np; ++p)
                resv_rate += voidage_rates[p];

            if (resv_rate < current_rate)
                return Well::ProducerCMode::RESV;
        }
    }

    if (controls.hasControl(Well::ProducerCMode::THP) && currentControl != Well::ProducerCMode::THP) {
        const auto& thp = well_.getTHPConstraint(summaryState);
        double current_thp = ws.thp;
        if (thp > current_thp && !ws.trivial_target) {
            // If WVFPEXP item 4 is set to YES1 or YES2
            // switching to THP is prevented if the well will
            // produce at a higher rate with THP control
            const auto& wvfpexp = well_.wellEcl().getWVFPEXP();
            bool rate_less_than_potential = true;
            if (wvfpexp.prevent()) {
                for (int p = 0; p < well_.numPhases(); ++p) {
                    // Currently we use the well potentials here computed before the iterations.
                    // We may need to recompute the well potentials to get a more
                    // accurate check here.
                    rate_less_than_potential = rate_less_than_potential && (-ws.surface_rates[p]) <= ws.well_potentials[p];
                }
            }
            if (!wvfpexp.prevent() || !rate_less_than_potential) {
                thp_limit_violated_but_not_switched = false;
                return Well::ProducerCMode::THP;
            } else {
                thp_limit_violated_but_not_switched = true;
                deferred_logger.info("NOT_SWITCHING_TO_THP",
                "The THP limit is violated for producer " +
                well_.name() +
                ". But the rate will increase if switched to THP. " +
                "The well is therefore kept at " + Well::ProducerCMode2String(currentControl));
            }
        }
    }

    return currentControl;
}

Well::InjectorCMode WellConstraints::
activeInjectionConstraint(const SingleWellState& ws,
                          const SummaryState& summaryState,
                          bool& thp_limit_violated_but_not_switched,
                          DeferredLogger& deferred_logger) const
{
    const PhaseUsage& pu = well_.phaseUsage();

    const auto controls = well_.wellEcl().injectionControls(summaryState);
    const auto currentControl = ws.injection_cmode;

    if (controls.hasControl(Well::InjectorCMode::BHP) && currentControl != Well::InjectorCMode::BHP)
    {
        const auto& bhp = controls.bhp_limit;
        double current_bhp = ws.bhp;
        if (bhp < current_bhp)
            return Well::InjectorCMode::BHP;
    }

    if (controls.hasControl(Well::InjectorCMode::RATE) && currentControl != Well::InjectorCMode::RATE)
    {
        InjectorType injectorType = controls.injector_type;
        double current_rate = 0.0;

        switch (injectorType) {
        case InjectorType::WATER:
        {
            current_rate = ws.surface_rates[ pu.phase_pos[BlackoilPhases::Aqua] ];
            break;
        }
        case InjectorType::OIL:
        {
            current_rate = ws.surface_rates[ pu.phase_pos[BlackoilPhases::Liquid] ];
            break;
        }
        case InjectorType::GAS:
        {
            current_rate = ws.surface_rates[  pu.phase_pos[BlackoilPhases::Vapour] ];
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injectors " + well_.name());
        }

        if (controls.surface_rate < current_rate)
            return Well::InjectorCMode::RATE;
    }

    if (controls.hasControl(Well::InjectorCMode::RESV) && currentControl != Well::InjectorCMode::RESV)
    {
        double current_rate = 0.0;
        if( pu.phase_used[BlackoilPhases::Aqua] )
            current_rate += ws.reservoir_rates[ pu.phase_pos[BlackoilPhases::Aqua] ];

        if( pu.phase_used[BlackoilPhases::Liquid] )
            current_rate += ws.reservoir_rates[ pu.phase_pos[BlackoilPhases::Liquid] ];

        if( pu.phase_used[BlackoilPhases::Vapour] )
            current_rate += ws.reservoir_rates[ pu.phase_pos[BlackoilPhases::Vapour] ];

        if (controls.reservoir_rate < current_rate)
            return Well::InjectorCMode::RESV;
    }

    if (controls.hasControl(Well::InjectorCMode::THP) && currentControl != Well::InjectorCMode::THP)
    {
        const auto& thp = well_.getTHPConstraint(summaryState);
        double current_thp = ws.thp;
        if (thp < current_thp) {
            bool rate_less_than_potential = true;
            for (int p = 0; p < well_.numPhases(); ++p) {
                // Currently we use the well potentials here computed before the iterations.
                // We may need to recompute the well potentials to get a more
                // accurate check here.
                rate_less_than_potential = rate_less_than_potential && (ws.surface_rates[p]) <= ws.well_potentials[p];
            }
            if (!rate_less_than_potential) {
                thp_limit_violated_but_not_switched = false;
                return Well::InjectorCMode::THP;
            } else {
                thp_limit_violated_but_not_switched = true;
                deferred_logger.debug("NOT_SWITCHING_TO_THP",
                "The THP limit is violated for injector " +
                well_.name() +
                ". But the rate will increase if switched to THP. " +
                "The well is therefore kept at " + Well::InjectorCMode2String(currentControl));
            }
        }
    }

    return currentControl;
}

bool WellConstraints::
checkIndividualConstraints(SingleWellState& ws,
                           const SummaryState& summaryState,
                           const RateConvFunc& calcReservoirVoidageRates,
                           bool& thp_limit_violated_but_not_switched,
                           DeferredLogger& deferred_logger) const
{
    if (well_.isProducer()) {
        auto new_cmode = this->activeProductionConstraint(ws, summaryState,
                                                          calcReservoirVoidageRates,
                                                          thp_limit_violated_but_not_switched,
                                                          deferred_logger);
        if (new_cmode != ws.production_cmode) {
            ws.production_cmode = new_cmode;
            return true;
        }
    }

    if (well_.isInjector()) {
        auto new_cmode = this->activeInjectionConstraint(ws, summaryState,
                                                        thp_limit_violated_but_not_switched,
                                                        deferred_logger);
        if (new_cmode != ws.injection_cmode) {
            ws.injection_cmode = new_cmode;
            return true;
        }
    }

    return false;
}

bool WellConstraints::checkMaxRatioLimitWell(const SingleWellState& ws,
                                             const double max_ratio_limit,
                                             const RatioFunc& ratioFunc) const
{
    const int np = well_.numPhases();

    std::vector<double> well_rates(np, 0.0);
    for (int p = 0; p < np; ++p) {
        well_rates[p] = ws.surface_rates[p];
    }

    const double well_ratio = ratioFunc(well_rates, well_.phaseUsage());
    return (well_ratio > max_ratio_limit);
}

void WellConstraints::
checkMaxRatioLimitCompletions(const SingleWellState& ws,
                              const double max_ratio_limit,
                              const RatioFunc& ratioFunc,
                              const ParallelWellInfo& parallel_well_info,
                              RatioLimitCheckReport& report) const
{
    int worst_offending_completion = std::numeric_limits<int>::max();

    // the maximum water cut value of the completions
    // it is used to identify the most offending completion
    double max_ratio_completion = 0;
    const int np = well_.numPhases();

    const auto& perf_data = ws.perf_data;
    const auto& perf_phase_rates = perf_data.phase_rates;
    // look for the worst_offending_completion
    for (const auto& completion : well_.getCompletions()) {
        std::vector<double> completion_rates(np, 0.0);

        // looping through the connections associated with the completion
        const std::vector<int>& conns = completion.second;
        for (const int c : conns) {
            for (int p = 0; p < np; ++p) {
                const double connection_rate = perf_phase_rates[c * np + p];
                completion_rates[p] += connection_rate;
            }
        } // end of for (const int c : conns)

        parallel_well_info.communication().sum(completion_rates.data(), completion_rates.size());
        const double ratio_completion = ratioFunc(completion_rates, well_.phaseUsage());

        if (ratio_completion > max_ratio_completion) {
            worst_offending_completion = completion.first;
            max_ratio_completion = ratio_completion;
        }
    } // end of for (const auto& completion : completions_)

    const double violation_extent = max_ratio_completion / max_ratio_limit;

    if (violation_extent > report.violation_extent) {
        report.worst_offending_completion = worst_offending_completion;
        report.violation_extent = violation_extent;
    }
}

void WellConstraints::
checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                 const SingleWellState& ws,
                 const ParallelWellInfo& parallel_well_info,
                 RatioLimitCheckReport& report) const
{
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    // function to calculate gor based on rates
    auto gor = [](const std::vector<double>& rates,
                  const PhaseUsage& pu) {
        const double oil_rate = -rates[pu.phase_pos[Oil]];
        const double gas_rate = -rates[pu.phase_pos[Gas]];
        if (gas_rate <= 0.)
            return 0.;
        else if (oil_rate <= 0.)
            return 1.e100; // big value to mark it as violated
        else
            return (gas_rate / oil_rate);
    };

    const double max_gor_limit = econ_production_limits.maxGasOilRatio();
    assert(max_gor_limit > 0.);

    const bool gor_limit_violated = this->checkMaxRatioLimitWell(ws, max_gor_limit, gor);

    if (gor_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_gor_limit, gor, parallel_well_info, report);
    }
}

void WellConstraints::
checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                 const SingleWellState& ws,
                 const ParallelWellInfo& parallel_well_info,
                 RatioLimitCheckReport& report) const
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Gas = BlackoilPhases::Vapour;

    // function to calculate wgr based on rates
    auto wgr = [](const std::vector<double>& rates,
                  const PhaseUsage& pu) {

        const double water_rate = -rates[pu.phase_pos[Water]];
        const double gas_rate = -rates[pu.phase_pos[Gas]];
        if (water_rate <= 0.)
            return 0.;
        else if (gas_rate <= 0.)
            return 1.e100; // big value to mark it as violated
        else
            return (water_rate / gas_rate);
    };

    const double max_wgr_limit = econ_production_limits.maxWaterGasRatio();
    assert(max_wgr_limit > 0.);

    const bool wgr_limit_violated = this->checkMaxRatioLimitWell(ws, max_wgr_limit, wgr);

    if (wgr_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_wgr_limit, wgr, parallel_well_info, report);
    }
}

void WellConstraints::
checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                      const SingleWellState& ws,
                      const ParallelWellInfo& parallel_well_info,
                      RatioLimitCheckReport& report) const
{
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    // function to calculate water cut based on rates
    auto waterCut = [](const std::vector<double>& rates,
                       const PhaseUsage& pu) {
        const double oil_rate = -rates[pu.phase_pos[Oil]];
        const double water_rate = -rates[pu.phase_pos[Water]];
        const double liquid_rate = oil_rate + water_rate;
        if (liquid_rate <= 0.)
            return 0.;
        else if (water_rate < 0)
            return 0.;
        else if (oil_rate < 0)
            return 1.;
        else
            return (water_rate / liquid_rate);

    };

    const double max_water_cut_limit = econ_production_limits.maxWaterCut();
    assert(max_water_cut_limit > 0.);

    const bool watercut_limit_violated = this->checkMaxRatioLimitWell(ws, max_water_cut_limit, waterCut);

    if (watercut_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_water_cut_limit, waterCut,
                                            parallel_well_info, report);
    }
}

bool WellConstraints::
checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                    const std::vector<double>& rates_or_potentials,
                    DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    const PhaseUsage& pu = well_.phaseUsage();

    if (econ_production_limits.onMinOilRate()) {
        const double oil_rate = rates_or_potentials[pu.phase_pos[ Oil ] ];
        const double min_oil_rate = econ_production_limits.minOilRate();
        if (std::abs(oil_rate) < min_oil_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinGasRate() ) {
        const double gas_rate = rates_or_potentials[pu.phase_pos[ Gas ] ];
        const double min_gas_rate = econ_production_limits.minGasRate();
        if (std::abs(gas_rate) < min_gas_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinLiquidRate() ) {
        const double oil_rate = rates_or_potentials[pu.phase_pos[ Oil ] ];
        const double water_rate = rates_or_potentials[pu.phase_pos[ Water ] ];
        const double liquid_rate = oil_rate + water_rate;
        const double min_liquid_rate = econ_production_limits.minLiquidRate();
        if (std::abs(liquid_rate) < min_liquid_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinReservoirFluidRate()) {
        deferred_logger.warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
    }

    return false;
}

} // namespace Opm
