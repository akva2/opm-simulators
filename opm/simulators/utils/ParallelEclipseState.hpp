/*
  Copyright 2019 Equinor AS.

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
#ifndef PARALLEL_ECLIPSE_STATE_HPP
#define PARALLEL_ECLIPSE_STATE_HPP

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <dune/common/parallel/mpihelper.hh>

namespace Opm {

class ParallelEclipseState : public EclipseState {
public:
    ParallelEclipseState() = default;
    ParallelEclipseState(const Deck& deck,
                         const ParseContext& parseContext,
                         ErrorGuard& errors);

    void setFP(const FieldPropsManager& fp) { field_props = fp; }

    bool broadcast(Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm);
};
    
} // end namespace Opm
#endif // PARALLEL_ECLIPSE_STATE_HPP
