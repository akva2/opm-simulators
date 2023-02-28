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
#include <opm/simulators/utils/SerializationPackers.hpp>

#include <boost/date_time/gregorian/gregorian.hpp>

namespace Opm {
namespace Serialization {
namespace detail {

std::size_t Packing<false,boost::gregorian::date>::
packSize(const boost::gregorian::date& data)
{
    return Packing<false,std::string>::packSize(boost::gregorian::to_simple_string(data));
}

void Packing<false,boost::gregorian::date>::
pack(const boost::gregorian::date& data,
     std::vector<char>& buffer, int& position)
{
    Packing<false,std::string>::pack(boost::gregorian::to_simple_string(data), buffer, position);
}

void Packing<false,boost::gregorian::date>::
unpack(boost::gregorian::date& data,
       std::vector<char>& buffer, int& position)
{
    std::string date;
    Packing<false,std::string>::unpack(date, buffer, position);
    data = boost::gregorian::from_simple_string(date);
}

#if DUNE_VERSION_LT(DUNE_COMMON, 2, 7)
template<class Scalar, int Size>
std::size_t Packing<false,Dune::FieldVector<Scalar,Size>>::
packSize(const Dune::FieldVector<Scalar,Size>& data)
{
    return Size*Packing<true,Scalar>::packSize(data[0]);
}

template<class Scalar, int Size>
void Packing<false,Dune::FieldVector<Scalar,Size>>::
pack(const Dune::FieldVector<Scalar,Size>& data,
     std::vector<char>& buffer, int& position)
{
    for (int i = 0; i < Size; ++i)
        Packing<true,Scalar>::pack(data[i], buffer, position);
}

template<class Scalar, int Size>
void Packing<false,Dune::FieldVector<Scalar,Size>>::
unpack(Dune::FieldVector<Scalar,Size>& data,
       std::vector<char>& buffer, int& position)
{
    for (int i = 0; i < Size; ++i) {
        Packing<true,Scalar>::unpack(data[i], buffer, position);
    }
}

template class Packing<false,Dune::FieldVector<double,1>>;
template class Packing<false,Dune::FieldVector<double,3>>;

#endif

} // end namespace detail
} // end namespace Serialization
} // end namespace Opm
