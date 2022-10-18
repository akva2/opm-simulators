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


#ifndef OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED
#define OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED

#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Opm
{

class DeferredLogger;
class SummaryState;
struct ThrowOnError;
class VFPProperties;
class Well;
class WellState;

//! \brief Class for computing BHP limits.
class WellBhpThpCalculator {
public:
    //! \brief Set to true to enable extra debug output.
    static constexpr bool extraBhpAtThpLimitOutput = false;

    //! \brief Checks if a well has THP constraints.
    static bool wellHasTHPConstraints(const SummaryState& summaryState,
                                      const Well& well_ecl);

    //! \brief Get THP constraint for a well.
    static double getTHPConstraint(const SummaryState& summaryState,
                                   const Well& well_ecl);

    //! \brief Obtain BHP limit for a well.
    static double mostStrictBhpFromBhpLimits(const SummaryState& summaryState,
                                             const Well& well_ecl);


    //! \brief Calculate THP from BHP.
    static double calculateThpFromBhp(const WellState &well_state,
                                      const std::vector<double>& rates,
                                      const double bhp,
                                      DeferredLogger& deferred_logger);

    //! \brief Compute BHP at THP limit for an injector.
    static std::optional<double>
    computeBhpAtThpLimitInj(const std::function<std::vector<double>(const double)>& frates,
                            const SummaryState& summary_state,
                            const VFPProperties& vfpProperties,
                            const Well& well_ecl,
                            const double rho,
                            const double refDepth,
                            const double gravity,
                            const bool throwOnError,
                            const int max_iteration,
                            const double flo_rel_tol,
                            DeferredLogger& deferred_logger);

    //! \brief Compute BHP at THP limit for a producer.
    static std::optional<double>
    computeBhpAtThpLimitProd(const std::function<std::vector<double>(const double)>& frates,
                             const SummaryState& summary_state,
                             const VFPProperties& vfpProperties,
                             const Well& well_ecl,
                             const double maxPerfPress,
                             const double rho,
                             const double alq_value,
                             const double refDepth,
                             const double gravity,
                             const int indexOfWell,
                             const bool useVfpExplicit,
                             DeferredLogger& deferred_logger);

protected:
    //! \brief Find the bhp-point where production becomes nonzero.
    static std::optional<double>
    bhpMax(const std::function<double(const double)>& fflo,
           const std::string& name,
           const double bhp_limit,
           const double maxPerfPress,
           const double vfp_flo_front,
           DeferredLogger& deferred_logger);

    //! \brief Compute BHP at THP limit.
    static std::optional<double>
    computeBhpAtThpLimit(const std::function<std::vector<double>(const double)>& frates,
                         const std::function<double(const std::vector<double>)>& fbhp,
                         const std::array<double, 2>& range,
                         const std::string& name,
                         DeferredLogger& deferred_logger);

    template<class ErrorPolicy>
    static std::optional<double>
    computeBhpAtThpLimitInjImpl(const std::function<std::vector<double>(const double)>& frates,
                                const SummaryState& summary_state,
                                const VFPProperties& vfpProperties,
                                const Well& well_ecl,
                                const double rho,
                                const double refDepth,
                                const double gravity,
                                const int max_iteration,
                                const double flo_rel_tol,
                                DeferredLogger& deferred_logger);

    //! \brief Brute-force solve for limits.
    static bool bruteForceBracket(const std::function<double(const double)>& eq,
                                  const std::array<double, 2>& range,
                                  double& low, double& high,
                                  DeferredLogger& deferred_logger);

    //! \brief Bisection solve for limits.
    static bool bisectBracket(const std::function<double(const double)>& eq,
                              const std::array<double, 2>& range,
                              const std::string& name,
                              double& low, double& high,
                              std::optional<double>& approximate_solution,
                              DeferredLogger& deferred_logger);
};

}

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
