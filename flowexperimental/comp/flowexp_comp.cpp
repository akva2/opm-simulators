/*
  Copyright 2024, SINTEF Digital

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

#include <dune/common/parallel/mpihelper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <opm/models/utils/start.hh>

#include "flowexp_comp.hpp"

int
main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowExpCompProblem<0>;
    Opm::registerEclTimeSteppingParameters<double>();

    // This is a bit cumbersome, but we need to read the input file
    // first to figure out the number of components (I think) in order
    // to select the correct type tag.
    //
    // TODO: Do a more dynamic dispatch approach similar to the normal
    //       flow application
    auto comm = Dune::MPIHelper::instance(argc, argv).getCommunication();
    auto commPtr = std::make_unique<decltype(comm)>(comm);
    Opm::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), true);

    auto inputFilename
        = Opm::FlowGenericVanguard::canonicalDeckPath(Opm::Parameters::Get<Opm::Parameters::EclDeckFileName>());
    Opm::FlowGenericVanguard::setCommunication(std::move(commPtr));
    Opm::FlowGenericVanguard::readDeck(inputFilename);
    Opm::FlowGenericVanguard vanguard;
    const auto numComps = vanguard.eclState().compositionalConfig().numComps();

    OPM_ERROR_IF(numComps < 2 || numComps < 7,
                 "Deck has {} components, not supported. We support a maximum of 7 components, "
                 "and a minimum of 2.");

    switch (numComps) {
    case 2: return Opm::dispatchFlowExpComp<2>(argc, argv);
    case 3: return Opm::dispatchFlowExpComp<3>(argc, argv);
    case 4: return Opm::dispatchFlowExpComp<4>(argc, argv);
    case 5: return Opm::dispatchFlowExpComp<5>(argc, argv);
    case 6: return Opm::dispatchFlowExpComp<6>(argc, argv);
    case 7: return Opm::dispatchFlowExpComp<7>(argc, argv);
    }
}

