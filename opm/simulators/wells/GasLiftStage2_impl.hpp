/*
  Copyright 2021 Equinor ASA.

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

#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>

#include <cmath>
#include <optional>
#include <string>

#include <fmt/format.h>

namespace Opm {

template<typename TypeTag>
GasLiftStage2<TypeTag>::
GasLiftStage2(const BlackoilWellModel<TypeTag>& well_model,
              const Simulator &ebos_simulator,
              DeferredLogger &deferred_logger,
              WellState &well_state,
              GLiftProdWells &prod_wells,
              GLiftOptWells &glift_wells,
              GLiftWellStateMap &state_map)
    : GasLiftStage2Generic(deferred_logger,
                           well_state,
                           state_map,
                           prod_wells,
                           glift_wells,
                           ebos_simulator.vanguard().schedule(),
                           ebos_simulator.vanguard().summaryState(),
                           ebos_simulator.vanguard().grid().comm(),
                           ebos_simulator.episodeIndex())
    , ebos_simulator_{ebos_simulator}
    , well_model_{well_model}
    , phase_usage_{well_model_.phaseUsage()}
{
    this->report_step_idx_ = ebos_simulator_.episodeIndex();
//    this->time_step_idx_
//        = this->ebos_simulator_.model().newtonMethod().currentTimeStep();
    this->nonlinear_iteration_idx_
        = this->ebos_simulator_.model().newtonMethod().numIterations();

}

/********************************************
 * Public methods in alphabetical order
 ********************************************/

// runOptimize():
//
// If a group has any production rate constraints, and/or a limit on
// its total rate of lift gas supply, allocates lift gas
// preferentially to the wells that gain the most benefit from
// it. Lift gas increments are allocated in turn to the well that
// currently has the largest weighted incremental gradient. The
// procedure takes account of any limits on the group production rate
// or lift gas supply applied to any level of group, including the FIELD level group.
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
runOptimize()
{
    const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);

    optimizeGroupsRecursive_(group);

}


/********************************************
 * Private methods in alphabetical order
 ********************************************/

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
recalculateGradientAndUpdateData_(
    GradPairItr &grad_itr, bool increase,

    //incremental and decremental gradients, if 'grads' are incremental, then
    // 'other_grads' are decremental, or conversely, if 'grads' are decremental, then
    // 'other_grads' are incremental
    std::vector<GradPair> &grads,  std::vector<GradPair> &other_grads)
{
    // NOTE: We make a copy of the name string instead of taking a reference
    //   since we may have to erase grad_itr (in the "else" condition below)
    const std::string name = grad_itr->first;
    std::optional<GradInfo> old_grad = std::nullopt;

    // only applies to wells in the well_state_map (i.e. wells on this rank)
    // the grads and other grads are synchronized later
    if(this->stage1_wells_.count(name) > 0) {
        GasLiftSingleWellGeneric &gs_well = *(this->stage1_wells_.at(name).get());
        auto grad = calcIncOrDecGrad_(name, gs_well, increase);
        if (grad) {
            grad_itr->second = grad->grad;
            old_grad = updateGrad_(name, *grad, increase);
        }
        else {
            grads.erase(grad_itr); // NOTE: this invalidates grad_itr
            old_grad = deleteGrad_(name, increase);
        }
    }

    if (old_grad) {
        // NOTE: Either creates a new item or reassigns
        // The old incremental gradient becomes the new decremental gradient
        //   or the old decremental gradient becomes the new incremental gradient
        updateGrad_(name, *old_grad, !increase);
        // NOTE: This may invalidate any iterator into 'other_grads' since
        //   updateGradVector_() will do a push_back() if 'name' is not found..
        updateGradVector_(name, other_grads, old_grad->grad);
    }
}

} // namespace Opm
