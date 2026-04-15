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
#ifndef OPM_FLOW_PROBLEM_BLACKOIL_IC_HPP
#define OPM_FLOW_PROBLEM_BLACKOIL_IC_HPP

#include <opm/simulators/flow/EquilInitializer.hpp>
#include <opm/simulators/flow/FlowProblemIC.hpp>

namespace Opm {

template<class TypeTag> class FlowProblemBlackoil;

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief Handling of initial conditions for FlowProblemBlackoil.
 */
template <class TypeTag>
class FlowProblemBlackoilIC : public FlowProblemIC<TypeTag>
{
public:
    FlowProblemBlackoilIC(FlowProblemBlackoil<TypeTag>& problem)
        : problem_(problem)
    {}

protected:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr bool enableBrine = getPropValue<TypeTag, Properties::EnableBrine>();
    static constexpr bool enableDisgasInWater = getPropValue<TypeTag, Properties::EnableDisgasInWater>();
    static constexpr bool enableDissolvedGas = Indices::compositionSwitchIdx >= 0;
    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();
    static constexpr bool enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>();
    static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
    static constexpr Scalar smallSaturationTolerance_ = 1.e-6;

    //! \brief Sets up equilibrium initial conditions.
    void equil_() override
    {
        // initial condition corresponds to hydrostatic conditions.
        EquilInitializer<TypeTag> equilInitializer(problem_.simulator(),
                                                   *problem_.materialLawManager());
        const std::size_t numElems = problem_.model().numGridDof();
        this->initialFluidStates_.resize(numElems);
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = this->initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }
    }

    void explicit_() override
    {
        const auto& fp = problem_.simulator().vanguard().eclState().fieldProps();
        const std::size_t numDof = problem_.model().numGridDof();

        struct FieldInfo
        {
            bool present;
            std::vector<double> data;
        };

        auto createField = [&fp, numDof](const std::string& name, bool cond, bool addzero = false)
        {
            const bool p = fp.has_double(name);
            return cond && p
                ? std::make_pair(name, FieldInfo{p, fp.get_double(name)})
                : std::make_pair(name, FieldInfo{p, (addzero ? std::vector<double>(numDof, 0.0)
                                                             : std::vector<double>{})});
        };

        const auto fields = std::map{
            createField("PRESSURE", true),
            createField("RS",   FluidSystem::enableDissolvedGas()),
            createField("RSW", FluidSystem::enableDissolvedGasInWater()),
            createField("RV", FluidSystem::enableVaporizedOil()),
            createField("RVW", FluidSystem::enableVaporizedWater()),
            createField("SALT", enableBrine),
            createField("SALTP", enableSaltPrecipitation),
            createField("SGAS", FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                                FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx), true),
            createField("SWAT", FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                                Indices::numPhases > 1, true),
            createField("TEMPI", true),
        };

        // make sure all required quantities are enables
        if (Indices::numPhases > 1) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && !fields.at("SWAT").present) {
                throw std::runtime_error("The ECL input file requires the presence of the SWAT keyword if "
                                         "the water phase is active");
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && !fields.at("SGAS").present &&
                FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
            {
                throw std::runtime_error("The ECL input file requires the presence of the SGAS keyword if "
                                         "the gas phase is active");
            }
        }
        if (!fields.at("PRESSURE").present) {
            throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                     "keyword if the model is initialized explicitly");
        }
        if (FluidSystem::enableDissolvedGas() && !fields.at("RS").present) {
            throw std::runtime_error("The ECL input file requires the RS keyword to be present if"
                                     " dissolved gas is enabled and the model is initialized explicitly");
        }
        if (FluidSystem::enableDissolvedGasInWater() && !fields.at("RSW").present) {
            OpmLog::warning("The model is initialized explicitly and the RSW keyword is not present in the"
                            " ECL input file. The RSW values are set equal to 0");
        }
        if (FluidSystem::enableVaporizedOil() && !fields.at("RV").present) {
            throw std::runtime_error("The ECL input file requires the RV keyword to be present if"
                                     " vaporized oil is enabled and the model is initialized explicitly");
        }
        if (FluidSystem::enableVaporizedWater() && !fields.at("RVW").present) {
            throw std::runtime_error("The ECL input file requires the RVW keyword to be present if"
                                     " vaporized water is enabled and the model is initialized explicitly");
        }
        if (enableBrine && !fields.at("SALT").present) {
            throw std::runtime_error("The ECL input file requires the SALT keyword to be present if"
                                     " brine is enabled and the model is initialized explicitly");
        }
        if (enableSaltPrecipitation && !fields.at("SALTP").present) {
            throw std::runtime_error("The ECL input file requires the SALTP keyword to be present if"
                                     " salt precipitation is enabled and the model is initialized explicitly");
        }

        this->initialFluidStates_.resize(numDof);

        // calculate the initial fluid states
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = this->initialFluidStates_[dofIdx];

            dofFluidState.setPvtRegionIndex(problem_.pvtRegionIndex(dofIdx));

            //////
            // set temperature
            //////
            if constexpr (energyModuleType != EnergyModules::NoTemperature) {
                Scalar temperatureLoc = fields.at("TEMPI").data[dofIdx];
                if (!std::isfinite(temperatureLoc) || temperatureLoc <= 0) {
                    temperatureLoc = FluidSystem::surfaceTemperature;
                }
                dofFluidState.setTemperature(temperatureLoc);
            }

            //////
            // set salt concentration
            //////
            if constexpr (enableBrine) {
                dofFluidState.setSaltConcentration(fields.at("SALT").data[dofIdx]);
            }

            //////
            // set precipitated salt saturation
            //////
            if constexpr (enableSaltPrecipitation) {
                dofFluidState.setSaltSaturation(fields.at("SALTP").data[dofIdx]);
            }

            //////
            // set saturations
            //////
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            fields.at("SWAT").data[dofIdx]);
            }

            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                1.0
                                                - fields.at("SWAT").data[dofIdx]);
                }
                else {
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                fields.at("SGAS").data[dofIdx]);
                }
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const Scalar soil = 1.0 - fields.at("SWAT").data[dofIdx] - fields.at("SGAS").data[dofIdx];
                if (soil < smallSaturationTolerance_) {
                    dofFluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
                }
                else {
                    dofFluidState.setSaturation(FluidSystem::oilPhaseIdx, soil);
                }
            }

            //////
            // set phase pressures
            // oil pressure (or gas pressure for water-gas system or water pressure for single phase)
            const Scalar pressure = fields.at("PRESSURE").data[dofIdx];

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            std::array<Scalar, FluidSystem::numPhases> pc = {0};
            const auto& matParams = problem_.materialLawParams(dofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, dofFluidState);
            Valgrind::CheckDefined(pressure);
            Valgrind::CheckDefined(pc);
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                if constexpr (Indices::oilEnabled) {
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::oilPhaseIdx]));
                }
                else if constexpr (Indices::gasEnabled) {
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[FluidSystem::gasPhaseIdx]));
                }
                else if constexpr (Indices::waterEnabled) {
                    //single (water) phase
                    dofFluidState.setPressure(phaseIdx, pressure);
                }
            }

            if constexpr (enableDissolvedGas) {
                if (FluidSystem::enableDissolvedGas()) {
                    dofFluidState.setRs(fields.at("RS").data[dofIdx]);
                }
                else if (Indices::gasEnabled && Indices::oilEnabled) {
                    dofFluidState.setRs(0.0);
                }
                if (FluidSystem::enableVaporizedOil()) {
                    dofFluidState.setRv(fields.at("RV").data[dofIdx]);
                }
                else if (Indices::gasEnabled && Indices::oilEnabled) {
                    dofFluidState.setRv(0.0);
                }
            }

            if constexpr (enableDisgasInWater) {
                if (FluidSystem::enableDissolvedGasInWater()) {
                    dofFluidState.setRsw(fields.at("RSW").data[dofIdx]);
                }
            }

            if constexpr (enableVapwat) {
                if (FluidSystem::enableVaporizedWater()) {
                    dofFluidState.setRvw(fields.at("RVW").data[dofIdx]);
                }
            }

            //////
            // set invB_
            //////
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const auto& b = FluidSystem::inverseFormationVolumeFactor(dofFluidState, phaseIdx,
                                                                          problem_.pvtRegionIndex(dofIdx));
                dofFluidState.setInvB(phaseIdx, b);

                const auto& rho = FluidSystem::density(dofFluidState, phaseIdx,
                                                       problem_.pvtRegionIndex(dofIdx));
                dofFluidState.setDensity(phaseIdx, rho);
            }
        }
    }

private:
    FlowProblemBlackoil<TypeTag>& problem_;
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_HPP
