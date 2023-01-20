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

#ifndef OPM_AQUIFERCT_HEADER_INCLUDED
#define OPM_AQUIFERCT_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferCT.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferAnalytical.hpp>
#include <opm/simulators/aquifers/AquiferCarterTracyGeneric.hpp>

#include <exception>
#include <memory>
#include <stdexcept>
#include <utility>

namespace Opm {

template <typename TypeTag>
class AquiferCarterTracy : public AquiferAnalytical<TypeTag>
                         , public AquiferCarterTracyGeneric<GetPropType<TypeTag,
                                                                        Properties::Scalar>>
{
public:
    using Base = AquiferAnalytical<TypeTag>;

    using typename Base::BlackoilIndices;
    using typename Base::ElementContext;
    using typename Base::Eval;
    using typename Base::FluidState;
    using typename Base::FluidSystem;
    using typename Base::IntensiveQuantities;
    using typename Base::RateVector;
    using typename Base::Scalar;
    using typename Base::Simulator;
    using typename Base::ElementMapper;

    AquiferCarterTracy(const std::vector<Aquancon::AquancCell>& connections,
                       const Simulator& ebosSimulator,
                       const AquiferCT::AQUCT_data& aquct_data)
        : Base(aquct_data.aquiferID, connections, ebosSimulator)
        , AquiferCarterTracyGeneric<Scalar>(aquct_data)
    {}

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->ebos_simulator_.timeStepSize();
        }
        this->fluxValue_ = this->W_flux_.value();
        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&this->fluxValue_, 1);
    }

    data::AquiferData aquiferData() const override
    {
        Scalar fluxRate = 0.;
        for (const auto& q : this->Qai_) {
            fluxRate += q.value();
        }

        return this->aquiferData_(this->aquiferID(),
                                  this->pa0_,
                                  fluxRate,
                                  this->W_flux_.value(),
                                  this->Tc_,
                                  this->rhow_,
                                  this->co2store_());
    }

protected:
    void assignRestartData(const data::AquiferData& xaq) override
    {
        this->rhow_ = this->assignRestartData_(xaq);
    }

    Scalar dpai(const int idx) const
    {
        const auto gdz =
            this->gravity_() * (this->cell_depth_.at(idx) - this->aquiferDepth());

        const auto dp = this->pa0_ + this->rhow_*gdz
            - this->pressure_previous_.at(idx);

        return dp;
    }

    // This function implements Eqs 5.8 and 5.9 of the EclipseTechnicalDescription
    std::pair<Scalar, Scalar>
    calculateEqnConstants(const int idx, const Simulator& simulator)
    {
        return this->calculateEqnConstants_(simulator.timeStepSize(),
                                            simulator.time(),
                                            this->Tc_,
                                            this->dpai(idx));
    }

    // This function implements Eq 5.7 of the EclipseTechnicalDescription
    void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        const auto [a, b] = this->calculateEqnConstants(idx, simulator);

        this->Qai_.at(idx) = this->alphai_.at(idx) *
            (a - b*(this->pressure_current_.at(idx) - this->pressure_previous_.at(idx)));
    }

    void calculateAquiferConstants() override
    {
        this->Tc_ = this->template calculateAquiferConstants_<FluidSystem>(this->co2store_());
    }

    inline void calculateAquiferCondition() override
    {
        if (this->solution_set_from_restart_) {
            return;
        }

        if (! this->aquct_data_.initial_pressure.has_value()) {
            this->aquct_data_.initial_pressure =
                this->calculateReservoirEquilibrium();

            const auto& tables = this->ebos_simulator_.vanguard()
                .eclState().getTableManager();

            this->aquct_data_.finishInitialisation(tables);
        }

        this->pa0_ = this->aquct_data_.initial_pressure.value();
        if (this->aquct_data_.initial_temperature.has_value())
            this->Ta0_ = this->aquct_data_.initial_temperature.value();

        this->rhow_ = this->template waterDensity_<FluidSystem>(this->co2store_());
    }

    virtual Scalar aquiferDepth() const override
    {
        return this->aquct_data_.datum_depth;
    }
}; // class AquiferCarterTracy

} // namespace Opm

#endif
