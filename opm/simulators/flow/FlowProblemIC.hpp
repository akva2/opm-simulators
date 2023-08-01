// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::FlowProblem
 */
#ifndef OPM_FLOW_PROBLEM_IC_HPP
#define OPM_FLOW_PROBLEM_IC_HPP

#include <opm/simulators/flow/EquilInitializer.hpp>

#include <vector>

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief Handling of initial conditions for FlowProblem.
 */
template <class TypeTag>
class FlowProblemIC
{
public:
    using EclMaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using InitialFluidState = typename EquilInitializer<TypeTag>::ScalarFluidState;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    static constexpr bool enableBrine = getPropValue<TypeTag, Properties::EnableBrine>();
    static constexpr bool enableDisgasInWater = getPropValue<TypeTag, Properties::EnableDisgasInWater>();
    static constexpr bool enableDissolvedGas = Indices::compositionSwitchIdx >= 0;
    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();
    static constexpr bool enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>();
    static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
    static constexpr Scalar smallSaturationTolerance_ = 1.e-6;

    //! \brief Returns a const reference to initial fluid state for an element.
    const InitialFluidState& initialFluidState(const unsigned idx) const
    { return initialFluidStates_[idx]; }

    //! \brief Sets up equilibrium initial conditions.
    void equilInitialCondition_(EclMaterialLawManager& materialLawManager,
                                const Simulator& simulator,
                                const std::size_t numElems)
    {
        // initial condition corresponds to hydrostatic conditions.
        using EquilInitializer = EquilInitializer<TypeTag>;
        EquilInitializer equilInitializer(simulator, materialLawManager);

        initialFluidStates_.resize(numElems);
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }
    }

    //! \brief Read explicitly specified initial conditions from eclipse state.
    void explicitInitialCondition_(const FieldPropsManager& fp,
                                   const EclMaterialLawManager& materialLawManager,
                                   const std::size_t numDof,
                                   std::function<int(int)> pvtRegionIndex)
    {
        bool has_swat     = fp.has_double("SWAT");
        bool has_sgas     = fp.has_double("SGAS");
        bool has_rs       = fp.has_double("RS");
        bool has_rsw      = fp.has_double("RSW");
        bool has_rv       = fp.has_double("RV");
        bool has_rvw      = fp.has_double("RVW");
        bool has_pressure = fp.has_double("PRESSURE");
        bool has_salt     = fp.has_double("SALT");
        bool has_saltp    = fp.has_double("SALTP");

        // make sure all required quantities are enables
        if (Indices::numPhases > 1) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && !has_swat) {
                throw std::runtime_error("The ECL input file requires the presence of the SWAT keyword if "
                                         "the water phase is active");
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                !has_sgas && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
            {
                throw std::runtime_error("The ECL input file requires the presence of the SGAS keyword if "
                                         "the gas phase is active");
            }
        }
        if (!has_pressure) {
            throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                     "keyword if the model is initialized explicitly");
        }
        if (FluidSystem::enableDissolvedGas() && !has_rs) {
            throw std::runtime_error("The ECL input file requires the RS keyword to be present if"
                                     " dissolved gas is enabled and the model is initialized explicitly");
        }
        if (FluidSystem::enableDissolvedGasInWater() && !has_rsw) {
            OpmLog::warning("The model is initialized explicitly and the RSW keyword is not present in the"
                            " ECL input file. The RSW values are set equal to 0");
        }
        if (FluidSystem::enableVaporizedOil() && !has_rv) {
            throw std::runtime_error("The ECL input file requires the RV keyword to be present if"
                                     " vaporized oil is enabled and the model is initialized explicitly");
        }
        if (FluidSystem::enableVaporizedWater() && !has_rvw) {
            throw std::runtime_error("The ECL input file requires the RVW keyword to be present if"
                                     " vaporized water is enabled and the model is initialized explicitly");
        }
        if (enableBrine && !has_salt) {
            throw std::runtime_error("The ECL input file requires the SALT keyword to be present if"
                                     " brine is enabled and the model is initialized explicitly");
        }
        if (enableSaltPrecipitation && !has_saltp)
            throw std::runtime_error("The ECL input file requires the SALTP keyword to be present if"
                                     " salt precipitation is enabled and the model is initialized explicitly");

        this->initialFluidStates_.resize(numDof);

        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> pressureData;
        std::vector<double> rsData;
        std::vector<double> rswData;
        std::vector<double> rvData;
        std::vector<double> rvwData;
        std::vector<double> tempiData;
        std::vector<double> saltData;
        std::vector<double> saltpData;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && Indices::numPhases > 1) {
            waterSaturationData = fp.get_double("SWAT");
        }
        else {
            waterSaturationData.resize(numDof);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
        {
            gasSaturationData = fp.get_double("SGAS");
        }
        else {
            gasSaturationData.resize(numDof);
        }

        pressureData = fp.get_double("PRESSURE");
        if (FluidSystem::enableDissolvedGas()) {
            rsData = fp.get_double("RS");
        }

        if (FluidSystem::enableDissolvedGasInWater() && has_rsw) {
            rswData = fp.get_double("RSW");
        }

        if (FluidSystem::enableVaporizedOil()) {
            rvData = fp.get_double("RV");
        }

        if (FluidSystem::enableVaporizedWater()) {
            rvwData = fp.get_double("RVW");
        }

        // initial reservoir temperature
        tempiData = fp.get_double("TEMPI");

        // initial salt concentration data
        if constexpr (enableBrine) {
            saltData = fp.get_double("SALT");
        }

        // initial precipitated salt saturation data
        if constexpr (enableSaltPrecipitation) {
            saltpData = fp.get_double("SALTP");
        }

        // calculate the initial fluid states
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = this->initialFluidStates_[dofIdx];

            dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            //////
            // set temperature
            //////
            if constexpr (energyModuleType != EnergyModules::NoTemperature) {
                Scalar temperatureLoc = tempiData[dofIdx];
                if (!std::isfinite(temperatureLoc) || temperatureLoc <= 0) {
                    temperatureLoc = FluidSystem::surfaceTemperature;
                }
                dofFluidState.setTemperature(temperatureLoc);
            }

            //////
            // set salt concentration
            //////
            if constexpr (enableBrine) {
                dofFluidState.setSaltConcentration(saltData[dofIdx]);
            }

            //////
            // set precipitated salt saturation
            //////
            if constexpr (enableSaltPrecipitation) {
                dofFluidState.setSaltSaturation(saltpData[dofIdx]);
            }

            //////
            // set saturations
            //////
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            waterSaturationData[dofIdx]);
            }

            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                1.0
                                                - waterSaturationData[dofIdx]);
                }
                else {
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                gasSaturationData[dofIdx]);
                }
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const Scalar soil = 1.0 - waterSaturationData[dofIdx] - gasSaturationData[dofIdx];
                if (soil < smallSaturationTolerance_) {
                    dofFluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
                }
                else {
                    dofFluidState.setSaturation(FluidSystem::oilPhaseIdx, soil);
                }
            }

            //////
            // set phase pressures
            //////
            Scalar pressure = pressureData[dofIdx]; // oil pressure (or gas pressure for water-gas system or water pressure for single phase)

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            std::array<Scalar, FluidSystem::numPhases> pc = {0};
            const auto& matParams = materialLawManager.materialLawParams(dofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, dofFluidState);
            Valgrind::CheckDefined(pressure);
            Valgrind::CheckDefined(pc);
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                if (Indices::oilEnabled) {
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::oilPhaseIdx]));
                }
                else if (Indices::gasEnabled) {
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::gasPhaseIdx]));
                }
                else if (Indices::waterEnabled) {
                    //single (water) phase
                    dofFluidState.setPressure(phaseIdx, pressure);
                }
            }

            if constexpr (enableDissolvedGas) {
                if (FluidSystem::enableDissolvedGas()) {
                    dofFluidState.setRs(rsData[dofIdx]);
                }
                else if (Indices::gasEnabled && Indices::oilEnabled) {
                    dofFluidState.setRs(0.0);
                }
                if (FluidSystem::enableVaporizedOil()) {
                    dofFluidState.setRv(rvData[dofIdx]);
                }
                else if (Indices::gasEnabled && Indices::oilEnabled) {
                    dofFluidState.setRv(0.0);
                }
            }

            if constexpr (enableDisgasInWater) {
                if (FluidSystem::enableDissolvedGasInWater() && has_rsw) {
                    dofFluidState.setRsw(rswData[dofIdx]);
                }
            }

            if constexpr (enableVapwat) {
                if (FluidSystem::enableVaporizedWater()) {
                    dofFluidState.setRvw(rvwData[dofIdx]);
                }
            }

            //////
            // set invB_
            //////
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const auto& b = FluidSystem::inverseFormationVolumeFactor(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setInvB(phaseIdx, b);

                const auto& rho = FluidSystem::density(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setDensity(phaseIdx, rho);

            }
        }
    }

    std::vector<InitialFluidState> initialFluidStates_; //!< Vector of initial fluid states

};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_HPP
