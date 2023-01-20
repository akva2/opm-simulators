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

#include <config.h>
#include <opm/simulators/aquifers/AquiferCarterTracyGeneric.hpp>

#include <opm/common/utility/numeric/linearInterpolation.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/output/data/Aquifer.hpp>

namespace Opm {

template<class Scalar>
Scalar AquiferCarterTracyGeneric<Scalar>::
assignRestartData_(const data::AquiferData& xaq)
{
    this->fluxValue_ = xaq.volume;
    return this->aquct_data_.waterDensity();
}

template<class Scalar>
std::pair<Scalar, Scalar> AquiferCarterTracyGeneric<Scalar>::
getInfluenceTableValues(const Scalar td_plus_dt)
{
    // We use the opm-common numeric linear interpolator
    this->dimensionless_pressure_ =
        linearInterpolation(this->aquct_data_.dimensionless_time,
                            this->aquct_data_.dimensionless_pressure,
                            this->dimensionless_time_);

    const auto PItd =
        linearInterpolation(this->aquct_data_.dimensionless_time,
                            this->aquct_data_.dimensionless_pressure,
                            td_plus_dt);

    const auto PItdprime =
        linearInterpolationDerivative(this->aquct_data_.dimensionless_time,
                                      this->aquct_data_.dimensionless_pressure,
                                      td_plus_dt);

    return std::make_pair(PItd, PItdprime);
}

template<class Scalar>
std::pair<Scalar, Scalar> AquiferCarterTracyGeneric<Scalar>::
calculateEqnConstants_(const Scalar timeStepSize,
                       const Scalar time,
                       const Scalar Tc,
                       const Scalar dpai)
{
    const Scalar td_plus_dt = (timeStepSize + time) / Tc;
    this->dimensionless_time_ = time / Tc;

    const auto [PItd, PItdprime] = this->getInfluenceTableValues(td_plus_dt);

    const auto denom = Tc * (PItd - this->dimensionless_time_ * PItdprime);
    const auto a = (this->beta_ * dpai - this->fluxValue_ * PItdprime) / denom;
    const auto b = this->beta_ / denom;

    return std::make_pair(a, b);
}

template<class Scalar>
data::AquiferData AquiferCarterTracyGeneric<Scalar>::
aquiferData_(const int aquiferID,
             const Scalar pa0,
             const Scalar flux,
             const Scalar volume,
             const Scalar Tc,
             const Scalar rhow,
             const bool co2store) const
{
    data::AquiferData data;
    data.aquiferID = aquiferID;
    // TODO: not sure how to get this pressure value yet
    data.pressure = pa0;
    data.fluxRate = flux;
    data.volume = volume;
    data.initPressure = pa0;

    auto* aquCT = data.typeData.template create<data::AquiferType::CarterTracy>();

    aquCT->dimensionless_time = this->dimensionless_time_;
    aquCT->dimensionless_pressure = this->dimensionless_pressure_;
    aquCT->influxConstant = this->aquct_data_.influxConstant();

    if (co2store) {
        aquCT->waterDensity = rhow;
        aquCT->timeConstant = Tc;
        const auto x = this->aquct_data_.porosity *
                       this->aquct_data_.total_compr *
                       this->aquct_data_.inner_radius *
                       this->aquct_data_.inner_radius;
        aquCT->waterViscosity = Tc *  this->aquct_data_.permeability / x;
    } else {
        aquCT->timeConstant = this->aquct_data_.timeConstant();
        aquCT->waterDensity = this->aquct_data_.waterDensity();
        aquCT->waterViscosity = this->aquct_data_.waterViscosity();
    }

    return data;
}

template<class Scalar>
template<class FluidSystem>
Scalar AquiferCarterTracyGeneric<Scalar>::
calculateAquiferConstants_(const bool co2store)
{
    this->beta_ = this->aquct_data_.influxConstant();
    if (co2store) {
         const auto press = this->aquct_data_.initial_pressure.value();
         Scalar temp = FluidSystem::reservoirTemperature();
         if (this->aquct_data_.initial_temperature.has_value())
             temp = this->aquct_data_.initial_temperature.value();

         Scalar rs = 0.0; // no dissolved CO2
         Scalar waterViscosity = FluidSystem::oilPvt().viscosity(this->pvtRegionIdx(),
                                                                 temp, press, rs);
         const auto x = this->aquct_data_.porosity *
                        this->aquct_data_.total_compr *
                        this->aquct_data_.inner_radius *
                        this->aquct_data_.inner_radius;
         return waterViscosity * x / this->aquct_data_.permeability;
    } else {
         return this->aquct_data_.timeConstant();
    }
}

template<class Scalar>
template<class FluidSystem>
Scalar AquiferCarterTracyGeneric<Scalar>::
waterDensity_(const bool co2store)
{
    if (co2store) {
         const auto press = this->aquct_data_.initial_pressure.value();

         Scalar temp = FluidSystem::reservoirTemperature();
         if (this->aquct_data_.initial_temperature.has_value())
             temp = this->aquct_data_.initial_temperature.value();

         constexpr Scalar rs = 0.0; // no dissolved CO2
         const Scalar waterDensity =
             FluidSystem::oilPvt().inverseFormationVolumeFactor(this->pvtRegionIdx(),
                                                                temp, press, rs) *
             FluidSystem::oilPvt().oilReferenceDensity(this->pvtRegionIdx());
         return waterDensity;
    } else {
         return this->aquct_data_.waterDensity();
    }
}

template class AquiferCarterTracyGeneric<double>;
using FS = BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>;
template double AquiferCarterTracyGeneric<double>::calculateAquiferConstants_<FS>(const bool);
template double AquiferCarterTracyGeneric<double>::waterDensity_<FS>(const bool);

} // namespace Opm
