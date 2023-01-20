/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS

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
#include <opm/simulators/aquifers/AquiferAnalyticalRestart.hpp>

#include <opm/output/data/Aquifer.hpp>

namespace Opm {

template<class Scalar>
std::optional<Scalar> AquiferAnalyticalRestart<Scalar>::
initRestart_(const std::map<int,data::AquiferData>& aquiferSoln,
             const int aquiferID)
{
    auto xaqPos = aquiferSoln.find(aquiferID);
    if (xaqPos == aquiferSoln.end())
        return std::nullopt;

    this->assignRestartData(xaqPos->second);

    this->pa0_ = xaqPos->second.initPressure;
    this->solution_set_from_restart_ = true;
    return xaqPos->second.volume;

}

template class AquiferAnalyticalRestart<double>;

} // namespace Opm
