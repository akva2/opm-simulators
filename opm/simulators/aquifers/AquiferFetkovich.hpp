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

#ifndef OPM_AQUIFETP_HEADER_INCLUDED
#define OPM_AQUIFETP_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Aquifer/Aquifetp.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferAnalytical.hpp>
#include <opm/simulators/aquifers/AquiferFetkovichGeneric.hpp>

#include <numeric>
#include <vector>

namespace Opm {

template <typename TypeTag>
class AquiferFetkovich : public AquiferAnalytical<TypeTag>
                       , public AquiferFetkovichGeneric<GetPropType<TypeTag, Properties::Scalar>>
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

    AquiferFetkovich(const std::vector<Aquancon::AquancCell>& connections,
                     const Simulator& ebosSimulator,
                     const Aquifetp::AQUFETP_data& aqufetp_data)
        : Base(aqufetp_data.aquiferID, connections, ebosSimulator)
        , AquiferFetkovichGeneric<Scalar>(aqufetp_data)
    {
    }

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->ebos_simulator_.timeStepSize();
        }
        this->aquifer_pressure_ = aquiferPressure();
    }

    data::AquiferData aquiferData() const override
    {
        Scalar fluxRate = std::accumulate(this->Qai_.begin(), this->Qai_.end(), 0.0,
                                          [](const double flux, const auto& q) -> double
                                          {
                                              return flux + q.value();
                                          });
        return this->aquiferData_(this->aquiferID(),
                                  fluxRate,
                                  this->W_flux_.value(),
                                  this->pa0_);
    }

protected:
    void assignRestartData(const data::AquiferData& xaq) override
    {
        this->rhow_ = this->assignRestartData_(xaq);
    }

    Eval dpai(int idx)
    {
        const auto gdz =
            this->gravity_() * (this->cell_depth_[idx] - this->aquiferDepth());

        return this->aquifer_pressure_ + this->rhow_*gdz
            - this->pressure_current_.at(idx);
    }

    // This function implements Eq 5.12 of the EclipseTechnicalDescription
    Scalar aquiferPressure()
    {
        Scalar Flux = this->W_flux_.value();

        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&Flux, 1);

        const auto denom =
            this->aqufetp_data_.total_compr * this->aqufetp_data_.initial_watvolume;

        return this->pa0_ - (Flux / denom);
    }

    void calculateAquiferConstants() override
    {
        this->Tc_ = this->aqufetp_data_.timeConstant();
    }

    // This function implements Eq 5.14 of the EclipseTechnicalDescription
    void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        const Scalar td_Tc_ = simulator.timeStepSize() / this->Tc_;
        const Scalar coef = (1 - exp(-td_Tc_)) / td_Tc_;

        this->Qai_.at(idx) = coef * this->alphai_[idx] *
            this->aqufetp_data_.prod_index * dpai(idx);
    }

    void calculateAquiferCondition() override
    {
        if (this->solution_set_from_restart_) {
            return;
        }

        if (! this->aqufetp_data_.initial_pressure.has_value()) {
            this->aqufetp_data_.initial_pressure =
                this->calculateReservoirEquilibrium();

            const auto& tables = this->ebos_simulator_.vanguard()
                .eclState().getTableManager();

            this->aqufetp_data_.finishInitialisation(tables);
        }

        this->rhow_ = this->aqufetp_data_.waterDensity();
        this->pa0_ = this->aqufetp_data_.initial_pressure.value();
        if (this->aqufetp_data_.initial_temperature.has_value())
            this->Ta0_ = this->aqufetp_data_.initial_temperature.value();
        this->aquifer_pressure_ = this->pa0_;
    }

    Scalar aquiferDepth() const override
    {
        return this->aqufetp_data_.datum_depth;
    }
}; // Class AquiferFetkovich

} // namespace Opm

#endif
