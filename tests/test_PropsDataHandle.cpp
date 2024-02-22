/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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
#include <dune/alugrid/common/fromtogridfactory.hh>

#define BOOST_TEST_MODULE TestPropsDataHandle
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/simulators/utils/PropsDataHandle.hpp>
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#if HAVE_DUNE_ALUGRID
#define DISABLE_ALUGRID_SFC_ORDERING
#include <dune/alugrid/3d/alugrid.hh>
#endif

#include <fmt/format.h>

namespace {

struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
    std::string errorstring;
    int errorcode;
};

void MPI_err_handler(MPI_Comm*, int* err_code, ...)
{
    std::vector<char> err_string(MPI_MAX_ERROR_STRING);
    int err_length;
    MPI_Error_string(*err_code, err_string.data(), &err_length);
    std::string s(err_string.data(), err_length);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *err_code);
}

bool
init_unit_test_func()
{
    return true;
}

}

BOOST_AUTO_TEST_CASE(CpGrid)
{
    Opm::Parallel::Communication        comm;
    std::shared_ptr<Opm::EclipseState>  eclState;
    std::shared_ptr<Opm::Schedule>      schedule;
    std::unique_ptr<Opm::UDQState>      udqState;
    std::unique_ptr<Opm::Action::State> actionState;
    std::unique_ptr<Opm::WellTestState> wtestState;
    std::shared_ptr<Opm::SummaryConfig> summaryConfig;

    Opm::readDeck(comm, "PropsTest.data", eclState, schedule,
                  udqState, actionState, wtestState,
                  summaryConfig,nullptr, "normal", false, false, {});

    Dune::CpGrid grid(comm);

    Opm::ParallelEclipseState& pEclState = static_cast<Opm::ParallelEclipseState&>(*eclState);

    grid.processEclipseFormat(comm.rank() == 0 ? &eclState->getInputGrid() : nullptr,
                              eclState.get(), false, false, false);

    {
        Opm::PropsDataHandle handle(grid, pEclState);
        grid.loadBalance(handle, Dune::uniformEdgeWgt, nullptr, true,
                         nullptr, false, true);
        grid.switchToDistributedView();
    }

    pEclState.switchToDistributedProps();

    const auto& prop = eclState->fieldProps().get_double("PORO");
    const auto& gv = grid.levelGridView(0);
    const auto& idSet = grid.globalIdSet();
    Dune::MultipleCodimMultipleGeomTypeMapper mapper(gv, Dune::mcmgElementLayout());
    for (const auto elem : Dune::elements(gv, Dune::Partitions::interiorBorderOverlap)) {
        BOOST_CHECK_EQUAL(prop[mapper.index(elem)], idSet.id(elem) + 1);
    }
}

#if HAVE_DUNE_ALUGRID
BOOST_AUTO_TEST_CASE(AluGrid)
{
    Opm::Parallel::Communication        comm;
    std::shared_ptr<Opm::EclipseState>  eclState;
    std::shared_ptr<Opm::Schedule>      schedule;
    std::unique_ptr<Opm::UDQState>      udqState;
    std::unique_ptr<Opm::Action::State> actionState;
    std::unique_ptr<Opm::WellTestState> wtestState;
    std::shared_ptr<Opm::SummaryConfig> summaryConfig;

    Opm::readDeck(comm, "PropsTest.data", eclState, schedule,
                  udqState, actionState, wtestState,
                  summaryConfig,nullptr, "normal", false, false, {});

    using AluGridType = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming,Dune::ALUGridMPIComm>;
    Dune::CpGrid grid(comm);

    Opm::ParallelEclipseState& pEclState = static_cast<Opm::ParallelEclipseState&>(*eclState);

    grid.processEclipseFormat(comm.rank() == 0 ? &eclState->getInputGrid() : nullptr,
                              eclState.get(), false, false, false);

    using Factory = Dune::FromToGridFactory<AluGridType>;
    std::vector<int> cartesianCellId;
    std::vector<unsigned int> ordering;
    auto agrid = Factory().convert(grid, cartesianCellId, ordering);
    {
        Opm::PropsDataHandle handle(*agrid, pEclState);
        agrid->loadBalance(handle);
    }

    pEclState.switchToDistributedProps();

    const auto& prop = eclState->fieldProps().get_double("PORO");
    const auto& gv = agrid->levelGridView(0);
    const auto& mgv = agrid->macroGridView();
    Dune::MultipleCodimMultipleGeomTypeMapper mapper(gv, Dune::mcmgElementLayout());
    if( comm.rank() == 2) {
      for (const auto& p : prop)
          std::cout << p << " ";
      std::cout << std::endl;
    }
    for (const auto& elem : Dune::elements(gv, Dune::Partitions::all)) {
        //std::cout << "hoi " << mapper.index(elem) << " " << prop.size() << std::endl;
        //BOOST_CHECK_EQUAL(prop[mapper.index(elem)], mgv.macroId(elem) + 1);
    }

/*    for (int i = 0; i < comm.size(); ++i) {
        if (i == comm.rank()) {
            std::cout << comm.rank() << " size " << grid.size(0) << "\n";
            const auto& prop = eclState->fieldProps().get_double("PORO");
            const auto& idSet = grid.localIdSet();
            const auto& gv = grid.levelGridView(0);
            auto mapper = Dune::MultipleCodimMultipleGeomTypeMapper(gv,
                                                                    Dune::mcmgElementLayout());
//            Opm::FvBaseElementContext<Opm::Properties::TTag::EbosTypeTag> ctx(gv,mapper);
            for (const auto elem : Dune::elements(gv, Dune::Partitions::interiorBorder)) {
 //               ctx.updateStencil(elem);
                const std::size_t id = idSet.id(elem);
                std::cout << fmt::format("Local id: {}, index: {}, props: {}, mapper: {}, size: {}\n",
                                         id, elem.index(), prop[elem.index()], mapper.index(elem), prop.size());
            }
            for (const double& d : prop) {
                 std::cout << d << " ";
            }
            std::cout << "\n======================" << std::endl;
        }
        comm.barrier();
    }*/
}
#endif

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
