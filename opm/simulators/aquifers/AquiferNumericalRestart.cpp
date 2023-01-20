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

#include <config.h>
#include <opm/simulators/aquifers/AquiferNumericalRestart.hpp>

#include <opm/output/eclipse/RestartValue.hpp>

namespace Opm {

template<class Scalar>
AquiferNumericalRestart<Scalar>::
AquiferNumericalRestart(const std::size_t size, const int aquiferID)
    : AquiferInterfaceRestart(aquiferID)
    , init_pressure_(size, 0.0)
{}

template<class Scalar>
void AquiferNumericalRestart<Scalar>::
initFromRestart_(const RestartValue& aquiferSoln)
{
    auto xaqPos = aquiferSoln.aquifer.find(aquiferID_);
    if (xaqPos == aquiferSoln.aquifer.end())
        return;

    if (this->connects_to_reservoir_) {
        this->cumulative_flux_ = xaqPos->second.volume;
    }

    if (const auto* aqData = xaqPos->second.typeData.template get<data::AquiferType::Numerical>();
        aqData != nullptr)
    {
        this->init_pressure_ = aqData->initPressure;
    }

    this->solution_set_from_restart_ = true;
}

template<class Scalar>
data::AquiferData AquiferNumericalRestart<Scalar>::
aquiferData() const
{
    data::AquiferData data;
    data.aquiferID = aquiferID_;
    data.pressure = this->pressure_;
    data.fluxRate = this->flux_rate_;
    data.volume = this->cumulative_flux_;

    auto* aquNum = data.typeData.template create<data::AquiferType::Numerical>();
    aquNum->initPressure = this->init_pressure_;

    return data;
}

template class AquiferNumericalRestart<double>;

} // namespace Opm
