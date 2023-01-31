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

template<class Scalar, int R, int C>
std::size_t Packing<false,Dune::FieldMatrix<Scalar,R,C>>::
packSize(const Dune::FieldMatrix<Scalar,R,C>& data)
{
    return R*C*Packing<true,Scalar>::packSize(data[0][0]);
}

template<class Scalar, int R, int C>
void Packing<false,Dune::FieldMatrix<Scalar,R,C>>::
pack(const Dune::FieldMatrix<Scalar,R,C>& data,
     std::vector<char>& buffer, int& position)
{
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
            Packing<true,Scalar>::pack(data[r][c], buffer, position);
}

template<class Scalar, int R, int C>
void Packing<false,Dune::FieldMatrix<Scalar,R,C>>::
unpack(Dune::FieldMatrix<Scalar,R,C>& data,
       std::vector<char>& buffer, int& position)
{
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c) {
            Packing<true,Scalar>::unpack(data[r][c], buffer, position);
        }
}

template struct Packing<false,Dune::FieldMatrix<double,3,3>>;

} // end namespace detail
} // end namespace Serialization
} // end namespace Opm
