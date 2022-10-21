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
#include <utility>
#include <vector>

namespace Opm
{

class DeferredLogger;
class PhaseUsage;
using RegionId = int;
class SingleWellState;
class WellInterfaceGeneric;

//! \brief Class for well constraints.
class WellConstraints {
public:
    //! \brief Constructor sets reference to well.
    WellConstraints(const WellInterfaceGeneric& well) : well_(well) {}

    using RateConvFunc = std::function<void(const RegionId, const int,
                                            const std::vector<double>&,
                                            std::vector<double>&)>;

    bool
    checkIndividualConstraints(SingleWellState& ws,
                               const SummaryState& summaryState,
                               const RateConvFunc& calcReservoirVoidageRates,
                               bool& thp_limit_violated_but_not_switched,
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

    const WellInterfaceGeneric& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_CONSTRAINTS_HEADER_INCLUDED
