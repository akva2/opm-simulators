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
#include <config.h>

#include "ParallelEclipseState.hpp"
#include "ParallelRestart.hpp"

namespace Opm {


ParallelEclipseState::ParallelEclipseState(const Deck& deck,
                                           const ParseContext& parseContext,
                                           ErrorGuard& errors)
    : EclipseState(deck, parseContext, errors)
{
}


bool ParallelEclipseState::broadcast(Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm) 
{
    if (comm.size() == 1)
        return true;
  
    if (comm.rank() == 0) {
        std::size_t size = Mpi::packSize(m_tables, comm) +
                           Mpi::packSize(m_runspec, comm) +
                           Mpi::packSize(m_eclipseConfig, comm) +
                           Mpi::packSize(m_deckUnitSystem, comm) +
                           Mpi::packSize(m_inputNnc, comm) +
                           Mpi::packSize(m_inputEditNnc, comm) +
                           Mpi::packSize(m_simulationConfig, comm) +
                           Mpi::packSize(m_transMult, comm) +
                           Mpi::packSize(m_faults, comm) +
                           Mpi::packSize(m_title, comm);
        std::vector<char> buffer(size); 
        int position = 0;
        Mpi::pack(m_tables, buffer, position, comm);
        Mpi::pack(m_runspec, buffer, position, comm);
        Mpi::pack(m_eclipseConfig, buffer, position, comm);
        Mpi::pack(m_deckUnitSystem, buffer, position, comm);
        Mpi::pack(m_inputNnc, buffer, position, comm);
        Mpi::pack(m_inputEditNnc, buffer, position, comm);
        Mpi::pack(m_simulationConfig, buffer, position, comm);
        Mpi::pack(m_transMult, buffer, position, comm);
        Mpi::pack(m_faults, buffer, position, comm);
        Mpi::pack(m_title, buffer, position, comm);
        comm.broadcast(&position, 1, 0);
        comm.broadcast(buffer.data(), position, 0);
    } else {
        int size;
        comm.broadcast(&size, 1, 0);
        std::vector<char> buffer(size);
        comm.broadcast(buffer.data(), size, 0);
        int position = 0;
        Mpi::unpack(m_tables, buffer, position, comm);
        Mpi::unpack(m_runspec, buffer, position, comm);
        Mpi::unpack(m_eclipseConfig, buffer, position, comm);
        Mpi::unpack(m_deckUnitSystem, buffer, position, comm);
        Mpi::unpack(m_inputNnc, buffer, position, comm);
        Mpi::unpack(m_inputEditNnc, buffer, position, comm);
        Mpi::unpack(m_simulationConfig, buffer, position, comm);
        Mpi::unpack(m_transMult, buffer, position, comm);
        Mpi::unpack(m_faults, buffer, position, comm);
        Mpi::unpack(m_title, buffer, position, comm);
    }

    return true;    
}


} // end namespace Opm
