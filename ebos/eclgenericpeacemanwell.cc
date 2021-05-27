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

#include <config.h>
#include <ebos/eclgenericpeacemanwell.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <iostream>

namespace Opm {

template<class FluidSystem, class Scalar, int numPhases>
EclGenericPeacemanWell<FluidSystem,Scalar,numPhases>::
EclGenericPeacemanWell(const Comm& comm)
    : comm_(comm)
{
    // set the initial status of the well
    wellType_ = Undefined;
    wellStatus_ = Shut;
    controlMode_ = BottomHolePressure;

    wellTotalVolume_ = 0.0;

    bhpLimit_ = 0.0;
    thpLimit_ = 0.0;

    targetBottomHolePressure_ = 0.0;
    actualBottomHolePressure_ = 0.0;
    maximumSurfaceRate_ = 0.0;
    maximumReservoirRate_ = 0.0;

    actualWeightedSurfaceRate_ = 0.0;
    actualWeightedResvRate_ = 0.0;
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        actualSurfaceRates_[phaseIdx] = 0.0;
        actualResvRates_[phaseIdx] = 0.0;

        volumetricWeight_[phaseIdx] = 0.0;
    }

    refDepth_ = 0.0;
    injectedPhaseIdx_ = oilPhaseIdx;
}

template<class FluidSystem, class Scalar, int numPhases>
void EclGenericPeacemanWell<FluidSystem,Scalar,numPhases>::
beginSpec()
{
    // this is going to be set to a real value by any realistic grid. Shall we bet?
    refDepth_ = 1e100;

    // By default, take the bottom hole pressure as a given
    controlMode_ = ControlMode::BottomHolePressure;

    // use one bar for the default bottom hole and tubing head
    // pressures. For the bottom hole pressure, this is probably
    // off by at least one magnitude...
    bhpLimit_ = 1e5;
    thpLimit_ = 1e5;

    // reset the actually observed bottom hole pressure
    actualBottomHolePressure_ = 0.0;

    // By default, all fluids exhibit the weight 1.0
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        volumetricWeight_[phaseIdx] = 1.0;

    wellType_ = Undefined;

    wellTotalVolume_ = 0.0;
}

template<class FluidSystem, class Scalar, int numPhases>
void EclGenericPeacemanWell<FluidSystem,Scalar,numPhases>::
endSpec(int nTotal)
{
    nTotal = comm_.sum(nTotal);
    if (nTotal == 0) {
        // well does not penetrate any active cell on any process. notify the
        // user about this.
        std::cout << "Well " << name() << " does not penetrate any active cell."
                  << " Assuming it to be shut!\n";
        setWellStatus(WellStatus::Shut);
    }

    // determine the maximum depth of the well over all processes
    refDepth_ = comm_.min(refDepth_);

    // the total volume of the well must also be summed over all processes
    wellTotalVolume_ = comm_.sum(wellTotalVolume_);
}

template<class FluidSystem, class Scalar, int numPhases>
void EclGenericPeacemanWell<FluidSystem,Scalar,numPhases>::
setVolumetricPhaseWeights(Scalar oilWeight, Scalar gasWeight, Scalar waterWeight)
{
    volumetricWeight_[oilPhaseIdx] = oilWeight;
    volumetricWeight_[gasPhaseIdx] = gasWeight;
    volumetricWeight_[waterPhaseIdx] = waterWeight;
}

template<class FluidSystem, class Scalar, int numPhases>
void EclGenericPeacemanWell<FluidSystem,Scalar,numPhases>::
beginTimeStep()
{
    if (wellStatus() == Shut)
        return;

    // calculate the bottom hole pressure to be actually used
    if (controlMode_ == ControlMode::TubingHeadPressure) {
        // assume a density of 650 kg/m^3 for the bottom hole pressure
        // calculation
        Scalar rho = 650.0;
        targetBottomHolePressure_ = thpLimit_ + rho*refDepth_;
    }
    else if (controlMode_ == ControlMode::BottomHolePressure)
        targetBottomHolePressure_ = bhpLimit_;
    else
        // TODO: also take the tubing head pressure limit into account...
        targetBottomHolePressure_ = bhpLimit_;

    // make it very likely that we screw up if we control for {surface,reservoir}
    // rate, but depend on the {reservoir,surface} rate somewhere...
    if (controlMode_ == ControlMode::VolumetricSurfaceRate)
        maximumReservoirRate_ = 1e100;
    else if (controlMode_ == ControlMode::VolumetricReservoirRate)
        maximumSurfaceRate_ = 1e100;

    // reset the iteration index
    iterationIdx_ = 0;
}


template<class FluidSystem, class Scalar, int numPhases>
void EclGenericPeacemanWell<FluidSystem,Scalar,numPhases>::
endTimeStep()
{
    if (wellStatus() == Shut)
        return;

    // we use a condition that is always false here to prevent the code below from
    // bitrotting. (i.e., at least it stays compileable)
    if (false && comm_.rank() == 0) {
        std::cout << "Well '" << name() << "':\n";
        std::cout << " Control mode: " << controlMode_ << "\n";
        std::cout << " BHP limit: " << bhpLimit_/1e5 << " bar\n";
        std::cout << " Observed BHP: " << actualBottomHolePressure_/1e5 << " bar\n";
        std::cout << " Weighted surface rate limit: " << maximumSurfaceRate_ << "\n";
        std::cout << " Weighted surface rate: " << std::abs(actualWeightedSurfaceRate_) << " (="
                  << 100*std::abs(actualWeightedSurfaceRate_)/maximumSurfaceRate_ << "%)\n";

        std::cout << " Surface rates:\n";
        std::cout << "  oil: "
                  << actualSurfaceRates_[oilPhaseIdx] << " m^3/s = "
                  << actualSurfaceRates_[oilPhaseIdx]*(24*60*60) << " m^3/day = "
                  << actualSurfaceRates_[oilPhaseIdx]*(24*60*60)/0.15898729 << " STB/day = "
                  << actualSurfaceRates_[oilPhaseIdx]*(24*60*60)
                     *FluidSystem::referenceDensity(oilPhaseIdx, /*pvtRegionIdx=*/0) << " kg/day"
                  << "\n";
        std::cout << "  gas: "
                  << actualSurfaceRates_[gasPhaseIdx] << " m^3/s = "
                  << actualSurfaceRates_[gasPhaseIdx]*(24*60*60) << " m^3/day = "
                  << actualSurfaceRates_[gasPhaseIdx]*(24*60*60)/28.316847 << " MCF/day = "
                  << actualSurfaceRates_[gasPhaseIdx]*(24*60*60)
                     *FluidSystem::referenceDensity(gasPhaseIdx, /*pvtRegionIdx=*/0) << " kg/day"
                  << "\n";
        std::cout << "  water: "
                  << actualSurfaceRates_[waterPhaseIdx] << " m^3/s = "
                  << actualSurfaceRates_[waterPhaseIdx]*(24*60*60) << " m^3/day = "
                  << actualSurfaceRates_[waterPhaseIdx]*(24*60*60)/0.15898729 << " STB/day = "
                  << actualSurfaceRates_[waterPhaseIdx]*(24*60*60)
                     *FluidSystem::referenceDensity(waterPhaseIdx, /*pvtRegionIdx=*/0) << " kg/day"
                  << "\n";
    }
}

template class EclGenericPeacemanWell<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                      double,
                                      3>;

} // namespace Opm
