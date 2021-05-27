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
/**
 * \file
 *
 * \copydoc Opm::EclPeacemanWell
 */
#ifndef EWOMS_ECL_GENERIC_PEACEMAN_WELL_HH
#define EWOMS_ECL_GENERIC_PEACEMAN_WELL_HH

#include <array>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/collectivecommunication.hh>

namespace Opm {

template<class FluidSystem, class Scalar, int numPhases>
class EclGenericPeacemanWell
{
public:
    // convenient access to the phase and component indices. If the compiler bails out
    // here, you're probably using an incompatible fluid system. This class has only been
    // tested with Opm::FluidSystems::BlackOil...
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned oilCompIdx = FluidSystem::oilCompIdx;
    static constexpr unsigned waterCompIdx = FluidSystem::waterCompIdx;
    static constexpr unsigned gasCompIdx = FluidSystem::gasCompIdx;

    using Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>;

    enum ControlMode {
        BottomHolePressure,
        TubingHeadPressure,
        VolumetricSurfaceRate,
        VolumetricReservoirRate
    };

    enum WellType {
        Undefined,
        Injector,
        Producer
    };

    enum WellStatus {
        // production/injection is ongoing
        Open,

        // no production/injection, but well is only closed above the reservoir, so cross
        // flow is possible
        Closed,

        // well is completely separated from the reservoir, e.g. by filling it with
        // concrete.
        Shut
    };

    EclGenericPeacemanWell(const Comm& comm);

    /*!
     * \brief Begin the specification of the well.
     *
     * The specification process is the following:
     *
     * beginSpec()
     * setName("FOO");
     * // add degrees of freedom to the well
     * for (dof in wellDofs)
     *    addDof(dof);
     * endSpec()
     *
     * // set the radius of the well at the dof [m].
     * // optional, if not specified, it is assumed to be 0.1524m
     * setRadius(dof, someRadius);
     *
     * // set the skin factor of the well.
     * // optional, if not specified, it is assumed to be 0
     * setSkinFactor(dof, someSkinFactor);
     *
     * // specify the phase which is supposed to be injected. (Optional,
     * // if unspecified, the well will throw an
     * // exception if it would inject something.)
     * setInjectedPhaseIndex(phaseIdx);
     *
     * // set maximum production rate at reservoir conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // shut for production)
     * setMaximumReservoirRate(someMassRate);
     *
     * // set maximum injection rate at reservoir conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // shut for injection)
     * setMinmumReservoirRate(someMassRate);
     *
     * // set the relative weight of the mass rate of a fluid phase.
     * // (Optional, if unspecified each phase exhibits a weight of 1)
     * setPhaseWeight(phaseIdx, someWeight);
     *
     * // set maximum production rate at surface conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // not limited by the surface rate)
     * setMaximumSurfaceRate(someMassRate);
     *
     * // set maximum production rate at surface conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // not limited by the surface rate)
     * setMinimumSurfaceRate(someMassRate);
     *
     * // set the minimum pressure at the bottom of the well (Pa,
     * // optional, if not specified, the well is assumes it estimates
     * // the bottom hole pressure based on the tubing head pressure
     * // assuming hydrostatic conditions.)
     * setMinimumBottomHolePressure(somePressure);
     *
     * // set the pressure at the top of the well (Pa,
     * // optional, if not specified, the tubing head pressure is
     * // assumed to be 1 bar)
     * setTubingHeadPressure(somePressure);
     *
     * // set the control mode of the well [m].
     * // optional, if not specified, it is assumed to be "BottomHolePressure"
     * setControlMode(Well::TubingHeadPressure);
     *
     * // set the tubing head pressure of the well [Pa]
     * // only require  if the control mode is "TubingHeadPressure"
     * setTubingHeadPressure(1e5);
     */
    void beginSpec();

    /*!
     * \brief Finalize the specification of the borehole.
     */
    void endSpec(int nTotal);

    /*!
     * \brief Set the relative weight of the volumetric phase rates.
     */
    void setVolumetricPhaseWeights(Scalar oilWeight, Scalar gasWeight, Scalar waterWeight);

    /*!
     * \brief Return the human-readable name of the well
     *
     * Well, let's say "readable by some humans".
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Set the human-readable name of the well
     */
    void setName(const std::string& newName)
    { name_ = newName; }

    /*!
     * \brief Set the control mode of the well.
     *
     * This specifies which quantities are assumed to be externally
     * given and which must be calculated based on those.
     */
    void setControlMode(ControlMode controlMode)
    { controlMode_ = controlMode; }

    /*!
     * \brief Set the temperature of the injected fluids [K]
     */
    void setTemperature(Scalar value)
    { wellTemperature_ = value; }

    /*!
     * \brief Set the type of the well (i.e., injector or producer).
     */
    void setWellType(WellType wellType)
    { wellType_ = wellType; }

    /*!
     * \brief Returns the type of the well (i.e., injector or producer).
     */
    WellType wellType() const
    { return wellType_; }

    /*!
     * \brief Set the index of fluid phase to be injected.
     *
     * This is only relevant if the well type is an injector.
     */
    void setInjectedPhaseIndex(unsigned injPhaseIdx)
    { injectedPhaseIdx_ = injPhaseIdx; }

    /*!
     * \brief Sets the reference depth for the bottom hole pressure [m]
     */
    void setReferenceDepth(Scalar value)
    { refDepth_ = value; }

    /*!
     * \brief The reference depth for the bottom hole pressure [m]
     */
    Scalar referenceDepth() const
    { return refDepth_; }

    /*!
     * \brief Set whether the well is open,closed or shut
     */
    void setWellStatus(WellStatus status)
    { wellStatus_ = status; }

    /*!
     * \brief Return whether the well is open,closed or shut
     */
    WellStatus wellStatus() const
    { return wellStatus_; }

    /*!
     * \brief Set the maximum/minimum bottom hole pressure [Pa] of the well.
     */
    void setTargetBottomHolePressure(Scalar val)
    { bhpLimit_ = val; }

    /*!
     * \brief Return the maximum/minimum bottom hole pressure [Pa] of the well.
     *
     * For injectors, this is the maximum, for producers it's the minimum.
     */
    Scalar targetBottomHolePressure() const
    { return bhpLimit_; }

    /*!
     * \brief Return the maximum/minimum bottom hole pressure [Pa] of the well.
     */
    Scalar bottomHolePressure() const
    { return actualBottomHolePressure_; }

    /*!
     * \brief Set the tubing head pressure [Pa] of the well.
     */
    void setTargetTubingHeadPressure(Scalar val)
    { thpLimit_ = val; }

    /*!
     * \brief Return the maximum/minimum tubing head pressure [Pa] of the well.
     *
     * For injectors, this is the maximum, for producers it's the minimum.
     */
    Scalar targetTubingHeadPressure() const
    { return thpLimit_; }

    /*!
     * \brief Return the maximum/minimum tubing head pressure [Pa] of the well.
     */
    Scalar tubingHeadPressure() const
    {
        // warning: this is a bit hacky...
        Scalar rho = 650; // kg/m^3
        Scalar g = 9.81; // m/s^2
        return actualBottomHolePressure_ + rho*refDepth_*g;
    }

    /*!
     * \brief Set the maximum combined rate of the fluids at the surface.
     */
    void setMaximumSurfaceRate(Scalar value)
    { maximumSurfaceRate_ = value; }

    /*!
     * \brief Return the weighted maximum surface rate [m^3/s] of the well.
     */
    Scalar maximumSurfaceRate() const
    { return maximumSurfaceRate_; }

    /*!
     * \brief Set the maximum combined rate of the fluids at the surface.
     */
    void setMaximumReservoirRate(Scalar value)
    { maximumReservoirRate_ = value; }

    /*!
     * \brief Return the weighted maximum reservoir rate [m^3/s] of the well.
     */
    Scalar maximumReservoirRate() const
    { return maximumReservoirRate_; }

    /*!
     * \brief Return the reservoir rate [m^3/s] actually seen by the well in the current time
     *        step.
     */
    Scalar reservoirRate() const
    { return actualWeightedResvRate_; }

    /*!
     * \brief Return the weighted surface rate [m^3/s] actually seen by the well in the current time
     *        step.
     */
    Scalar surfaceRate() const
    { return actualWeightedSurfaceRate_; }

    /*!
     * \brief Return the reservoir rate [m^3/s] of a given fluid which is actually seen
     *        by the well in the current time step.
     */
    Scalar reservoirRate(unsigned phaseIdx) const
    { return actualResvRates_[phaseIdx]; }

    /*!
     * \brief Return the weighted surface rate [m^3/s] of a given fluid which is actually
     *        seen by the well in the current time step.
     */
    Scalar surfaceRate(unsigned phaseIdx) const
    { return actualSurfaceRates_[phaseIdx]; }

    /*!
     * \brief Informs the well that a time step has just begun.
     */
    void beginTimeStep();

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    { ++iterationIdx_; }

    /*!
     * \brief Called by the simulator after each time step.
     */
    void endTimeStep();

protected:
    const Comm& comm_;

    std::string name_;

    // The assumed bottom hole and tubing head pressures as specified by the user
    Scalar bhpLimit_;
    Scalar thpLimit_;

    // The bottom hole pressure to be targeted by the well model. This may be computed
    // from the tubing head pressure (if the control mode is TubingHeadPressure), or it may be
    // just the user-specified bottom hole pressure if the control mode is
    // BottomHolePressure.
    Scalar targetBottomHolePressure_;

    // The bottom hole pressure which is actually observed in the well
    Scalar actualBottomHolePressure_;

    // the sum of the total volumes of all the degrees of freedoms that interact with the well
    Scalar wellTotalVolume_;

    // The volumetric surface rate which is actually observed in the well
    Scalar actualWeightedSurfaceRate_;
    std::array<Scalar, numPhases> actualSurfaceRates_;

    // The volumetric reservoir rate which is actually observed in the well
    Scalar actualWeightedResvRate_;
    std::array<Scalar, numPhases> actualResvRates_;

    // The relative weight of the volumetric rate of each fluid
    std::array<Scalar, numPhases> volumetricWeight_;

    // the type of the well (injector, producer or undefined)
    WellType wellType_;

    // Specifies whether the well is currently open, closed or shut. The difference
    // between "closed" and "shut" is that for the former, the well is assumed to be
    // closed above the reservoir so that cross-flow within the well is possible while
    // the well is completely separated from the reservoir if it is shut. (i.e., no
    // crossflow is possible in this case.)
    WellStatus wellStatus_;

    // specifies the quantities which are controlled for (i.e., which
    // should be assumed to be externally specified and which should
    // be computed based on those)
    ControlMode controlMode_;

    // the reference depth for the bottom hole pressure. if not specified otherwise, this
    // is the position of the _highest_ DOF in the well.
    Scalar refDepth_;

    // the temperature assumed for the fluid (in the case of an injector well)
    Scalar wellTemperature_;

    // The maximum weighted volumetric surface rates specified by the
    // user. This is used to apply rate limits and it is to be read as
    // the maximum absolute value of the rate, i.e., the well can
    // produce or inject the given amount.
    Scalar maximumSurfaceRate_;

    // The maximum weighted volumetric reservoir rates specified by
    // the user. This is used to apply rate limits and it is to be
    // read as the maximum absolute value of the rate, i.e., the well
    // can produce or inject the given amount.
    Scalar maximumReservoirRate_;

    // the number of times beginIteration*() was called for the current time step
    unsigned iterationIdx_;

    unsigned injectedPhaseIdx_;
};

} // namespace Opm

#endif
