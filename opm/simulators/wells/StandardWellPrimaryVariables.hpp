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

template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class StandardWellPrimaryVariables {
public:
    static constexpr int numStaticWellEq = StandardWellEquations<Indices,Scalar>::numStaticWellEq;
    static constexpr int Bhp = StandardWellEquations<Indices,Scalar>::Bhp;
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

    //! \brief Update polymer molecular weight values from solution vector.
    void updatePolyMW(const BVectorWell& dwells);

    //! \brief Copy values to well state.
    void copyToWellStatePolyMW(WellState& well_state) const;

private:
    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_; //!< Reference to well interface
};

}

#endif // OPM_STANDARDWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
