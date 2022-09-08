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
#ifndef MPI_SERIALIZER_HPP
#define MPI_SERIALIZER_HPP

#include <opm/common/utility/TimeService.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <bitset>
#include <cstddef>
#include <string>

namespace Opm {
namespace Mpi {
namespace detail {

template <bool pod, class T>
struct Packing
{
    static std::size_t packSize(const T&, Parallel::MPIComm);
    static void pack(const T&, std::vector<char>&, int&, Parallel::MPIComm);
    static void unpack(T&, std::vector<char>&, int&, Parallel::MPIComm);
};

//! \brief Packaging for pod data.
template<class T>
struct Packing<true,T>
{
    static std::size_t packSize(const T&, Parallel::MPIComm comm)
    {
        int size = 0;
        MPI_Pack_size(1, Dune::MPITraits<T>::getType(), comm, &size);
        return size;
    }

    static void pack(const T& data,
                     std::vector<char>& buffer,
                     int& position,
                     Parallel::MPIComm comm)
    {
        MPI_Pack(&data, 1, Dune::MPITraits<T>::getType(), buffer.data(),
                 buffer.size(), &position, comm);
    }

    static void unpack(T& data,
                       std::vector<char>& buffer,
                       int& position,
                       Parallel::MPIComm comm)
    {
        MPI_Unpack(buffer.data(), buffer.size(), &position, &data, 1,
                   Dune::MPITraits<T>::getType(), comm);
    }
};

//! \brief "Packging" for unsupported types.
template<class T>
struct Packing<false,T>
{
    static std::size_t packSize(const T&, Parallel::MPIComm)
    {
        static_assert(!std::is_same_v<T,T>, "Packing not supported for type");
        return 0;
    }

    static void pack(const T&, std::vector<char>&, int&,
                     Parallel::MPIComm)
    {
      static_assert(!std::is_same_v<T,T>, "Packing not supported for type");
    }

    static void unpack(T&, std::vector<char>&, int&,
                       Parallel::MPIComm)
    {
        static_assert(!std::is_same_v<T,T>, "Packing not supported for type");
    }
};

//! \brief Specialization for std::bitset
template <std::size_t Size>
struct Packing<false,std::bitset<Size>>
{
    static std::size_t packSize(const std::bitset<Size>&, Opm::Parallel::MPIComm);
    static void pack(const std::bitset<Size>&, std::vector<char>&, int&, Opm::Parallel::MPIComm);
    static void unpack(std::bitset<Size>&, std::vector<char>&, int&, Opm::Parallel::MPIComm);
};

#define ADD_PACK_SPECIALIZATION(T) \
    template<> \
    struct Packing<false,T> \
    { \
        static std::size_t packSize(const T&, Parallel::MPIComm); \
        static void pack(const T&, std::vector<char>&, int&, Parallel::MPIComm); \
        static void unpack(T&, std::vector<char>&, int&, Parallel::MPIComm); \
    };

ADD_PACK_SPECIALIZATION(std::string)
ADD_PACK_SPECIALIZATION(time_point)

}

template<class T>
std::size_t packSize(const T& data, Parallel::MPIComm comm)
{
    return detail::Packing<std::is_pod_v<T>,T>::packSize(data,comm);
}

template<class T>
void pack(const T& data,
          std::vector<char>& buffer,
          int& position,
          Parallel::MPIComm comm)
{
    detail::Packing<std::is_pod_v<T>,T>::pack(data, buffer, position, comm);
}

template<class T>
void unpack(T& data,
            std::vector<char>& buffer,
            int& position,
            Parallel::MPIComm comm)
{
    detail::Packing<std::is_pod_v<T>,T>::unpack(data, buffer, position, comm);
}

} // end namespace Mpi
} // end namespace Opm

#endif // MPI_SERIALIZER_HPP
