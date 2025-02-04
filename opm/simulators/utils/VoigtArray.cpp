/*
  Copyright 2023 Equinor ASA

  This file is part of the Open Porous Media Project (OPM).

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
#include <opm/simulators/utils/VoigtArray.hpp>

#include <algorithm>

namespace Opm {

template<class Scalar>
VoigtArray<Scalar>::
VoigtArray(const std::size_t size)
{
    this->resize(size);
}

template<class Scalar>
void VoigtArray<Scalar>::
resize(const std::size_t size)
{
    std::for_each(data_.begin(), data_.end(),
                  [size](auto& d) { d.resize(size); });
}

template<class Scalar>
const std::vector<Scalar>&
VoigtArray<Scalar>::
operator[](const VoigtIndex idx) const
{
    return data_[static_cast<std::underlying_type_t<VoigtIndex>>(idx)];
}

template<class Scalar>
std::vector<Scalar>&
VoigtArray<Scalar>::
operator[](const VoigtIndex idx)
{
    return data_[static_cast<std::underlying_type_t<VoigtIndex>>(idx)];
}

template<class Scalar>
Scalar
VoigtArray<Scalar>::
operator()(const VoigtIndex idx, const std::size_t i) const
{
    return data_[static_cast<std::underlying_type_t<VoigtIndex>>(idx)].at(i);
}

template<class Scalar>
Scalar&
VoigtArray<Scalar>::
operator()(const VoigtIndex idx, const std::size_t i)
{
    return data_[static_cast<std::underlying_type_t<VoigtIndex>>(idx)].at(i);
}

template class VoigtArray<double>;

#if FLOW_INSTANTIATE_FLOAT
template class VoigtArray<float>;
#endif

} // namespace Opm
