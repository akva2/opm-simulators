/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELL_CONSTRAINTS_HEADER_INCLUDED
#define OPM_WELL_CONSTRAINTS_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace Opm
{

class DeferredLogger;
class ParallelWellInfo;
class PhaseUsage;
using RegionId = int;
class Rates;
class SingleWellState;
class WellInterfaceGeneric;

struct RatioLimitCheckReport {
    bool ratio_limit_violated = false;
    int worst_offending_completion = std::numeric_limits<int>::max();
    double violation_extent = 0.0;
};

//! \brief Class for computing well group constraints.
class WellConstraints {
public:
    //! \brief Constructor sets reference to well.
    WellConstraints(const WellInterfaceGeneric& well) : well_(well) {}

    using RateConvFunc = std::function<void(const RegionId, const int,
                                            const std::vector<double>&,
                                            std::vector<double>&)>;

    using RatioFunc = std::function<double(const std::vector<double>& rates,
                                           const PhaseUsage& pu)>;

    bool
    checkIndividualConstraints(SingleWellState& ws,
                               const SummaryState& summaryState,
                               const RateConvFunc& calcReservoirVoidageRates,
                               bool& thp_limit_violated_but_not_switched,
                               DeferredLogger& deferred_logger) const;

    bool checkMaxRatioLimitWell(const SingleWellState& ws,
                                const double max_ratio_limit,
                                const RatioFunc& ratioFunc) const;

    void checkMaxRatioLimitCompletions(const SingleWellState& ws,
                                       const double max_ratio_limit,
                                       const RatioFunc& ratioFunc,
                                       const ParallelWellInfo& parallel_well_info,
                                       RatioLimitCheckReport& report) const;

    bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                            const std::vector<double>& rates_or_potentials,
                            DeferredLogger& deferred_logger) const;

    RatioLimitCheckReport checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                               const SingleWellState& ws,
                                               const ParallelWellInfo& parallel_well_info,
                                               DeferredLogger& deferred_logger) const;

private:
    Well::ProducerCMode
    activeProductionConstraint(const SingleWellState& ws,
                               const SummaryState& summaryState,
                               const RateConvFunc& calcReservoirVoidageRates,
                               bool& thp_limit_violated_but_not_switched,
                               DeferredLogger& deferred_logger) const;

    Well::InjectorCMode
    activeInjectionConstraint(const SingleWellState& ws,
                              const SummaryState& summaryState,
                              bool& thp_limit_violated_but_not_switched,
                              DeferredLogger& deferred_logger) const;

    void checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState& ws,
                          const ParallelWellInfo& parallel_well_info,
                          RatioLimitCheckReport& report) const;

    void checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState& ws,
                          const ParallelWellInfo& parallel_well_info,
                          RatioLimitCheckReport& report) const;

    void checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                               const SingleWellState& ws,
                               const ParallelWellInfo& parallel_well_info,
                               RatioLimitCheckReport& report) const;

    const WellInterfaceGeneric& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_CONSTRAINTS_HEADER_INCLUDED
