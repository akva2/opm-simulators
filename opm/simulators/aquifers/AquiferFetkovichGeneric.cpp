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
#include <opm/simulators/aquifers/AquiferFetkovichGeneric.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <stdexcept>

namespace Opm {

template<class Scalar>
Scalar AquiferFetkovichGeneric<Scalar>::
assignRestartData_(const data::AquiferData& xaq)
{
    if (!xaq.typeData.is<data::AquiferType::Fetkovich>()) {
        throw std::invalid_argument {
            "Analytic aquifer data for unexpected aquifer "
            "type passed to Fetkovich aquifer"
        };
    }

    this->aquifer_pressure_ = xaq.pressure;
    return this->aqufetp_data_.waterDensity();
}

template<class Scalar>
data::AquiferData AquiferFetkovichGeneric<Scalar>::
aquiferData_(const int aquiferID,
             const Scalar flux,
             const Scalar volume,
             const Scalar pa0) const
{
    // TODO: how to unify the two functions?
    auto data = data::AquiferData{};

    data.aquiferID = aquiferID;
    data.pressure = this->aquifer_pressure_;
    data.fluxRate = flux;
    data.volume = volume;
    data.initPressure = pa0;

    auto* aquFet = data.typeData.template create<data::AquiferType::Fetkovich>();
    aquFet->initVolume = this->aqufetp_data_.initial_watvolume;
    aquFet->prodIndex = this->aqufetp_data_.prod_index;
    aquFet->timeConstant = this->aqufetp_data_.timeConstant();

    return data;
}

template class AquiferFetkovichGeneric<double>;

} // namespace Opm
