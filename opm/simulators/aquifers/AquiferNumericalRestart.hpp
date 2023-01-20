/*
  Copyright (C) 2020 Equinor ASA
  Copyright (C) 2020 SINTEF Digital

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

#ifndef OPM_AQUIFERNUMERICAL_RESTART_HEADER_INCLUDED
#define OPM_AQUIFERNUMERICAL_RESTART_HEADER_INCLUDED

#include <cstddef>
#include <map>
#include <vector>

namespace Opm {

namespace data { class AquiferData; }

template<class Scalar>
class AquiferNumericalRestart {
protected:
    AquiferNumericalRestart(const std::size_t size);

    void initFromRestart_(const std::map<int,data::AquiferData>& aquiferSoln,
                          const int aquiferID);
    data::AquiferData aquiferData_(const int aquiferID) const;

    std::vector<Scalar> init_pressure_;
    bool solution_set_from_restart_ {false};
    bool connects_to_reservoir_ {false};
    Scalar cumulative_flux_{0.0}; // cumulative aquifer influx
    Scalar flux_rate_{0.0}; // aquifer influx rate
    Scalar pressure_{0.0}; // aquifer pressure
};

} // namespace Opm

#endif
