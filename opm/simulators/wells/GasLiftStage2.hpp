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

#ifndef OPM_GASLIFT_STAGE2_HEADER_INCLUDED
#define OPM_GASLIFT_STAGE2_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/simulators/wells/GasLiftStage2Generic.hpp>

#include <cassert>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#include <fmt/format.h>

namespace Opm
{

template<class TypeTag> class BlackoilWellModel;
template<class TypeTag> class WellInterface;

    template<class TypeTag>
    class GasLiftStage2 : public GasLiftStage2Generic
    {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        using GLiftOptWells = std::map<std::string,std::unique_ptr<GasLiftSingleWellGeneric>>;
        using GLiftProdWells = std::map<std::string,const WellInterfaceGeneric*>;
        using GLiftWellStateMap = std::map<std::string,std::unique_ptr<GasLiftWellState>>;
        using GradPair = std::pair<std::string, double>;
        using GradPairItr = std::vector<GradPair>::iterator;
        using GradInfo = GasLiftSingleWellGeneric::GradInfo;
        using GradMap = std::map<std::string, GradInfo>;
        using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
        using Communication = Dune::Communication<MPIComm>;
#else
        using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif
    public:
        GasLiftStage2(
            const BlackoilWellModel<TypeTag> &well_model,
            const Simulator &ebos_simulator,
            DeferredLogger &deferred_logger,
            WellState &well_state,
            GLiftProdWells &prod_wells,
            GLiftOptWells &glift_wells,
            GLiftWellStateMap &state_map
        );
        void runOptimize();
    private:
        void recalculateGradientAndUpdateData_(
            GradPairItr &grad_itr, bool increase,
            std::vector<GradPair> &grads, std::vector<GradPair> &other_grads);

        const Simulator &ebos_simulator_;
        const BlackoilWellModel<TypeTag>& well_model_;

        const PhaseUsage &phase_usage_;
        //int time_step_idx_;
        int nonlinear_iteration_idx_;
    };

} // namespace Opm

#include "GasLiftStage2_impl.hpp"

#endif // OPM_GASLIFT_STAGE2_HEADER_INCLUDED
