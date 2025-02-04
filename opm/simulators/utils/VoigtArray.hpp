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

#ifndef OPM_UTIL_VOIGT_ARRAY_HPP
#define OPM_UTIL_VOIGT_ARRAY_HPP

#include <array>
#include <vector>

namespace Opm {

enum class VoigtIndex {
    XX =  0, XY =  3, XZ = 4,
    YX = XY, YY =  1, YZ = 5,
    ZX = XZ, ZY = YZ, ZZ = 2,
};

template<class Scalar>
class VoigtArray
{
public:
    VoigtArray() = default;
    explicit VoigtArray(const std::size_t size);

    void resize(const std::size_t size);

    const std::vector<Scalar>& operator[](const VoigtIndex idx) const;
    std::vector<Scalar>& operator [](const VoigtIndex idx);

    constexpr std::size_t size() const { return 6; }

    auto begin()
    { return data_.begin(); }
    auto end()
    { return data_.end(); }
    auto begin() const
    { return data_.begin(); }
    auto end() const
    { return data_.end(); }

    Scalar operator()(const VoigtIndex idx, const std::size_t i) const;
    Scalar& operator()(const VoigtIndex idx, const std::size_t i);

private:
    std::array<std::vector<Scalar>, 6> data_{};
};

} // namespace Opm

#endif // OPM_UTIL_VOIGT_ARRAY_HPP
