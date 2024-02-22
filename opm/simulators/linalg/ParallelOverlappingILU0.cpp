/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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

#include <opm/simulators/linalg/ParallelOverlappingILU0_impl.hpp>

#include <dune/istl/owneroverlapcopy.hh>

namespace Opm
{

#define INSTANCE_PAR(T, Dim, ...) \
  template class ParallelOverlappingILU0<Dune::BCRSMatrix<MatrixBlock<T,Dim,Dim>>, \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>, \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>, \
                                         __VA_ARGS__>; \
  template class ParallelOverlappingILU0<Dune::BCRSMatrix<Dune::FieldMatrix<T,Dim,Dim>>, \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>, \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>, \
                                         __VA_ARGS__>;

#if HAVE_MPI
#define INSTANCE(T, Dim) \
    INSTANCE_PAR(T, Dim, Dune::Amg::SequentialInformation) \
    INSTANCE_PAR(T, Dim, Dune::OwnerOverlapCopyCommunication<int,int>)
#else
#define INSTANCE(T, Dim) \
    INSTANCE_PAR(T, Dim, Dune::Amg::SequentialInformation)
#endif

#define INSTANCE_TYPE(T) \
    INSTANCE(T,1)        \
    INSTANCE(T,2)        \
    INSTANCE(T,3)        \
    INSTANCE(T,4)        \
    INSTANCE(T,5)        \
    INSTANCE(T,6)

INSTANCE_TYPE(double)

#if FLOW_INSTANCE_FLOAT
INSTANCE_TYPE(float)
#endif

} // end namespace Opm
