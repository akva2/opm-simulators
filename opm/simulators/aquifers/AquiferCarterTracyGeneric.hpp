/*
  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
  Copyright 2017 Statoil ASA.

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

#ifndef OPM_AQUIFERCT_GENERIC_HEADER_INCLUDED
#define OPM_AQUIFERCT_GENERIC_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferCT.hpp>

#include <utility>

namespace Opm {

namespace data { class AquiferData; }

template<class Scalar>
class AquiferCarterTracyGeneric {
protected:
    AquiferCarterTracyGeneric(const AquiferCT::AQUCT_data& aquct_data)
        : aquct_data_(aquct_data)
    {}

    Scalar assignRestartData_(const data::AquiferData& xaq);
    data::AquiferData aquiferData_(const int aquiferID,
                                   const Scalar pa0,
                                   const Scalar flux,
                                   const Scalar volume,
                                   const Scalar Tc,
                                   const Scalar rhow,
                                   const bool co2store) const;

    std::pair<Scalar, Scalar> getInfluenceTableValues(const Scalar td_plus_dt);
    std::pair<Scalar, Scalar> calculateEqnConstants_(const Scalar timeStepSize,
                                                     const Scalar time,
                                                     const Scalar Tc,
                                                     const Scalar dpai);

    template<class FluidSystem>
    Scalar calculateAquiferConstants_(const bool co2store);
    template<class FluidSystem>
    Scalar waterDensity_(const bool co2store);

    std::size_t pvtRegionIdx() const
    {
        return this->aquct_data_.pvttableID - 1;
    }

    // Variables constants
    AquiferCT::AQUCT_data aquct_data_;

    Scalar beta_{0}; // Influx constant
    // TODO: it is possible it should be a AD variable
    Scalar fluxValue_{0}; // value of flux

    Scalar dimensionless_time_{0};
    Scalar dimensionless_pressure_{0};
};

} // namespace Opm

#endif
