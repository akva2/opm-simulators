/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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


#ifndef OPM_STANDARDWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
#define OPM_STANDARDWELL_PRIMARY_VARIABLES_HEADER_INCLUDED

#include <opm/material/densead/DynamicEvaluation.hpp>
#include <opm/simulators/wells/StandardWellEquations.hpp>

#include <vector>

namespace Opm
{

class DeferredLogger;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class StandardWellPrimaryVariables {
public:
    static constexpr bool has_wfrac_variable = Indices::waterEnabled && Indices::oilEnabled;
    static constexpr bool has_gfrac_variable = Indices::gasEnabled && Indices::numPhases > 1;

    static constexpr int WQTotal = 0;
    static constexpr int numStaticWellEq = StandardWellEquations<Indices,Scalar>::numStaticWellEq;
    static constexpr int Bhp = StandardWellEquations<Indices,Scalar>::Bhp;
    static constexpr int WFrac = has_wfrac_variable ? 1 : -1000;
    static constexpr int GFrac = has_gfrac_variable ? has_wfrac_variable + 1 : -1000;
    static constexpr int SFrac = Indices::enableSolvent ? 3 : -1000;

    using EvalWell = DenseAd::DynamicEvaluation<Scalar, numStaticWellEq + Indices::numEq + 1>;
    using BVectorWell = typename StandardWellEquations<Indices,Scalar>::BVectorWell;

    StandardWellPrimaryVariables(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well)
        : well_(well)
    {}

    // the values for the primary varibles
    // based on different solutioin strategies, the wells can have different primary variables
    std::vector<double> value_;

    // the Evaluation for the well primary variables, which contain derivatives and are used in AD calculation
    std::vector<EvalWell> evaluation_;

    //! \brief Initialize evaluations from values.
    void init(const int numWellEq);

    //! \brief Resize values and evaluations.
    void resize(const int numWellEq);

    //! \brief Copy values from well state.
    void update(const WellState& well_state, DeferredLogger& deferred_logger);

    //! \brief Update values from newton update vector.
    void updateNewton(const BVectorWell& dwells,
                      const double dFLimit,
                      const double dBHPLimit);

    //! \brief Update polymer molecular weight values from newton update vector.
    void updatePolyMW(const BVectorWell& dwells);

    //! \brief Copy values to well state.
    void copyToWellState(WellState& well_state, DeferredLogger& deferred_logger) const;

    //! \brief Copy polymer molecular weight values to well state.
    void copyToWellStatePolyMW(WellState& well_state) const;

    EvalWell wellVolumeFractionScaled(const int compIdx,
                                      const int numWellEq) const;
    EvalWell wellSurfaceVolumeFraction(const int compIdx,
                                       const int numWellEq) const;
    EvalWell getQs(const int compIdx,
                   const int numWellEq) const;

private:
    EvalWell wellVolumeFraction(const unsigned compIdx,
                                const int numWellEq) const;

    //! \brief Calculate a relaxation factor for producers.
    //! \details To avoid overshoot of the fractions which might result in negative rates.
    double relaxationFactorFractionsProducer(const std::vector<double>& primary_variables,
                                             const BVectorWell& dwells) const;

    //! \brief Handle non-reasonable fractions due to numerical overshoot.
    void processFractions();

    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_; //!< Reference to well interface
};

}

#endif // OPM_STANDARDWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
