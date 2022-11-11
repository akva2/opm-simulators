/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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


#ifndef OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED
#define OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED

#include <vector>

namespace Opm
{

class WellInterfaceGeneric;

template<class Scalar>
class StandardWellConnections
{
public:
    StandardWellConnections(const WellInterfaceGeneric& well);

    void computeConnectionPressureDelta();

    Scalar getRho() const
    {
        return this->perf_densities_.empty() ? 0.0 : perf_densities_[0];
    }

    // densities of the fluid in each perforation
    std::vector<Scalar> perf_densities_;
    // pressure drop between different perforations
    std::vector<Scalar> perf_pressure_diffs_;

private:
    // Base interface reference
    const WellInterfaceGeneric& well_;
};

}

#endif // OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED
