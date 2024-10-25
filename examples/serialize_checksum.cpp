/*
  Copyright 2020 Equinor.

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

#define LOG_CHECKSUM 1

#include <opm/common/utility/MemPacker.hpp>
#include <opm/common/utility/Serializer.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ASTNode.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQASTNode.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/WDFAC.hpp>
#include <opm/input/eclipse/Schedule/Well/WellBrineProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMICPProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>

#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#include <functional>
#include <iostream>
#include <string>
#include <type_traits>

#include <fmt/format.h>

template<class L, class R = L>
using has_equality = std::is_invocable<std::equal_to<>, bool, L, R>;
template<class L, class R = L>
constexpr auto has_equality_v = has_equality<L, R>::value;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Need one parameter, the .DATA file to run serialization checksums for\n";
        return 1;
    }

    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    Opm::Parallel::Communication comm{MPI_COMM_SELF};
#else
    Opm::Parallel::Communication comm{};
#endif

    Opm::setupLogging(comm, argv[1], "", "false", false, "", false);

    std::shared_ptr<Opm::EclipseState>  eclipseState;
    std::shared_ptr<Opm::Schedule>      schedule;
    std::unique_ptr<Opm::UDQState>      udqState;
    std::unique_ptr<Opm::Action::State> actionState;
    std::unique_ptr<Opm::WellTestState> wtestState;
    std::shared_ptr<Opm::SummaryConfig> summaryConfig;
    std::shared_ptr<Opm::Python>        python;

    Opm::readDeck(comm, argv[1], eclipseState, schedule, udqState,
                  actionState, wtestState, summaryConfig, python,
                  "low", "low", "100", false, true, false, {});

    Opm::Serialization::MemPacker packer;
    Opm::Serializer<Opm::Serialization::MemPacker> ser(packer, true);

    auto serializeAndDeserialize =
        [&ser](const auto& orig, const std::string& label)
        {
            using T = std::remove_const_t<std::remove_reference_t<decltype(orig)>>;
            T copy{};
            const auto csum1 = ser.checksum(orig);
            ser.pack(orig);
            ser.unpack(copy);
            const auto csum2 = ser.checksum(copy);
//            if constexpr (has_equality_v<T>) {
                fmt::print("{:>14}  {:^10}  {:^10}  {:^5} ({})\n",
                           label, csum1, csum2,
                           csum1 == csum2 ? "ok" : "fail", orig == copy);
            // } else {
            //     fmt::print("{:>14}  {:^10}  {:^10}  {:^5}\n", label, csum1, csum2,
            //                csum1 == csum2 ? "ok" : "fail");
            // }
        };

    fmt::print("\n\n    Object       Before       After     Status\n"
               "-----------------------------------------------\n");
    serializeAndDeserialize(*actionState, "Action::State");
    //serializeAndDeserialize(*eclipseState, "EclipseState");
    serializeAndDeserialize(*schedule, "Schedule");
    serializeAndDeserialize(*summaryConfig, "SummaryConfig");
    serializeAndDeserialize(*udqState, "UDQState");
    serializeAndDeserialize(*wtestState, "WellTestState");

    return 0;
}
