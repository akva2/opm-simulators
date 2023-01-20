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

#ifndef OPM_AQUIFETP_GENERIC_HEADER_INCLUDED
#define OPM_AQUIFETP_GENERIC_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Aquifer/Aquifetp.hpp>

#include <opm/simulators/aquifers/AquiferInterfaceRestart.hpp>

namespace Opm {

namespace data { struct AquiferData; }

template<class Scalar>
class AquiferFetkovichGeneric : public AquiferInterfaceRestart {
protected:
    AquiferFetkovichGeneric(const Aquifetp::AQUFETP_data& aqufetp_data,
                            const int aquiferID)
        : AquiferInterfaceRestart(aquiferID)
        , aqufetp_data_(aqufetp_data)
    {}

    data::AquiferData aquiferData() const override;

    Scalar assignRestartData_(const data::AquiferData& xaq);

    virtual Scalar getFlux() const = 0;
    virtual Scalar getVolumeFlux() const = 0;
    virtual Scalar getInitialPressure() const = 0;

    // Aquifer Fetkovich Specific Variables
    Aquifetp::AQUFETP_data aqufetp_data_;
    Scalar aquifer_pressure_; // aquifer
};

} // namespace Opm

#endif
