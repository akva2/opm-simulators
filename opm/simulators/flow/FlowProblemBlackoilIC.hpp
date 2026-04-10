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

private:
    FlowProblemBlackoil<TypeTag>& problem_;
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_HPP
