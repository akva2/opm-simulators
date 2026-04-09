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
#ifndef OPM_FLOW_PROBLEM_COMP_IC_HPP
#define OPM_FLOW_PROBLEM_COMP_IC_HPP

#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/simulators/flow/FlowProblemIC.hpp>

#include <vector>

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief Handling of initial conditions for FlowProblemComp.
 */
template <class TypeTag>
class FlowProblemCompIC : public FlowProblemIC<TypeTag>
{
public:
    using EclMaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    //! \brief Sets up equilibrium initial conditions.
    void equilInitialCondition_(EclMaterialLawManager&,
                                const Simulator&,
                                const std::size_t) override
    {
        throw std::logic_error("Equilibration is not supported by compositional modeling yet");
    }

    //! \brief Read explicitly specified initial conditions from eclipse state.
    void explicitInitialCondition_(const FieldPropsManager& fp,
                                   const EclMaterialLawManager&,
                                   const std::size_t numDof,
                                   std::function<int(int)>) override
    {
        const bool has_pressure = fp.has_double("PRESSURE");
        if (!has_pressure) {
            throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                     "keyword if the model is initialized explicitly");
        }

        const bool has_xmf = fp.has_double("XMF");
        const bool has_ymf = fp.has_double("YMF");
        const bool has_zmf = fp.has_double("ZMF");
        if ( !has_zmf && !(has_xmf && has_ymf) ) {
            throw std::runtime_error("The ECL input file requires the presence of ZMF or XMF and YMF "
                                     "keyword if the model is initialized explicitly");
        }

        if (has_zmf && (has_xmf || has_ymf)) {
            throw std::runtime_error("The ECL input file can not handle explicit initialization "
                                     "with both ZMF and XMF or YMF");
        }

        if (has_xmf != has_ymf) {
            throw std::runtime_error("The ECL input file needs XMF and YMF combined to do the explicit "
                                     "initializtion when using XMF or YMF");
        }

        const bool has_temp = fp.has_double("TEMPI");

        // const bool has_gas = fp.has_double("SGAS");
        assert(fp.has_double("SGAS"));

        this->initialFluidStates_.resize(numDof);

        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> soilData;
        std::vector<double> pressureData;
        std::vector<double> tempiData;

        const bool water_active = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        const bool gas_active = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
        const bool oil_active = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);

        if (water_active && Indices::numPhases > 2) {
            waterSaturationData = fp.get_double("SWAT");
        }
        else {
            waterSaturationData.resize(numDof);
        }

        pressureData = fp.get_double("PRESSURE");

        if (has_temp) {
            tempiData = fp.get_double("TEMPI");
        }
        else {
            ; // TODO: throw?
        }

        if (gas_active) { // && FluidSystem::phaseIsActive(oilPhaseIdx))
            gasSaturationData = fp.get_double("SGAS");
        }
        else {
            gasSaturationData.resize(numDof);
        }

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = this->initialFluidStates_[dofIdx];
            // dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            Scalar temperatureLoc = tempiData[dofIdx];
            assert(std::isfinite(temperatureLoc) && temperatureLoc > 0);
            dofFluidState.setTemperature(temperatureLoc);

            if (gas_active) {
                dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                            gasSaturationData[dofIdx]);
            }
            if (oil_active) {
                dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                            1.0
                                            - waterSaturationData[dofIdx]
                                            - gasSaturationData[dofIdx]);
            }
            if (water_active) {
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            waterSaturationData[dofIdx]);
            }

            //////
            // set phase pressures
            //////
            const Scalar pressure = pressureData[dofIdx]; // oil pressure (or gas pressure for water-gas system or water pressure for single phase)

            // TODO: zero capillary pressure for now
            const std::array<Scalar, FluidSystem::numPhases> pc = {0};
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
                    // single (water) phase
                    dofFluidState.setPressure(phaseIdx, pressure);
                }
            }

            if (has_xmf && has_ymf) {
                const auto& xmfData = fp.get_double("XMF");
                const auto& ymfData = fp.get_double("YMF");
                for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx) {
                    const std::size_t data_idx = compIdx * numDof + dofIdx;
                    const Scalar xmf = xmfData[data_idx];
                    const Scalar ymf = ymfData[data_idx];

                    dofFluidState.setMoleFraction(FluidSystem::oilPhaseIdx, compIdx, xmf);
                    dofFluidState.setMoleFraction(FluidSystem::gasPhaseIdx, compIdx, ymf);
                }
            }

            if (has_zmf) {
                zmf_initialization_ = true;
                const auto& zmfData = fp.get_double("ZMF");
                for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx) {
                    const std::size_t data_idx = compIdx * numDof + dofIdx;
                    const Scalar zmf = zmfData[data_idx];
                    dofFluidState.setMoleFraction(compIdx, zmf);

                    if (gas_active) {
                        const auto ymf = (dofFluidState.saturation(FluidSystem::gasPhaseIdx) > 0.) ? zmf : Scalar{0};
                        dofFluidState.setMoleFraction(FluidSystem::gasPhaseIdx, compIdx, ymf);
                    }
                    if (oil_active) {
                        const auto xmf = (dofFluidState.saturation(FluidSystem::oilPhaseIdx) > 0.) ? zmf : Scalar{0};
                        dofFluidState.setMoleFraction(FluidSystem::oilPhaseIdx, compIdx, xmf);
                    }
                }
            }
        }
    }

    bool zmf_initialization_{false};
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_HPP
