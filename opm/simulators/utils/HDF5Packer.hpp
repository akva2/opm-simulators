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
#ifndef HDF5_PACKER_HPP
#define HDF5_PACKER_HPP

#include <opm/common/utility/SimplePacker.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

namespace boost { namespace gregorian { class date; } }

namespace Opm {
namespace Serialization {
namespace detail {

//! \brief Specialization for Dune::FieldMatrix
template <class Scalar, int R, int C>
struct Packing<false,Dune::FieldMatrix<Scalar,R,C>>
{
    static std::size_t packSize(const Dune::FieldMatrix<Scalar,R,C>& data);

    static void pack(const Dune::FieldMatrix<Scalar,R,C>& data,
                     std::vector<char>& buffer, int& position);

    static void unpack(Dune::FieldMatrix<Scalar,R,C>& data,
                       std::vector<char>& buffer, int& position);
};

//! \brief Specialization for Dune:FieldVector:
template <class Scalar, int Size>
struct Packing<false,Dune::FieldVector<Scalar,Size>>
{
    static std::size_t packSize(const Dune::FieldVector<Scalar,Size>& data);

    static void pack(const Dune::FieldVector<Scalar,Size>& data,
                     std::vector<char>& buffer, int& position);

    static void unpack(Dune::FieldVector<Scalar,Size>& data,
                       std::vector<char>& buffer, int& position);
};

template<>
struct Packing<false,boost::gregorian::date>
{
    static std::size_t packSize(const boost::gregorian::date& data);

    static void pack(const boost::gregorian::date& data,
                     std::vector<char>& buffer, int& position);

    static void unpack(boost::gregorian::date& data,
                       std::vector<char>& buffer, int& position);
};

}

} // end namespace Serialization
} // end namespace Opm

#endif // HDF5_PACKER_HPP
