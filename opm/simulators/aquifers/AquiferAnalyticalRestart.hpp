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

#ifndef OPM_AQUIFERANALYTICAL_RESTART_HEADER_INCLUDED
#define OPM_AQUIFERANALYTICAL_RESTART_HEADER_INCLUDED

#include <optional>

namespace Opm {

namespace data { struct AquiferData; }
class RestartValue;

template<class Scalar>
class AquiferAnalyticalRestart {
protected:
    virtual void assignRestartData(const data::AquiferData& xaq) = 0;

    std::optional<Scalar> initRestart_(const RestartValue& aquiferSoln,
                                       const int aquiferID);

    Scalar pa0_{}; // initial aquifer pressure
    bool solution_set_from_restart_ {false};
};

} // namespace Opm

#endif
