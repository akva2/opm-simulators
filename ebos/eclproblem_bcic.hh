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
 * \copydoc Opm::EclProblem
 */
#ifndef OPM_ECL_PROBLEM_BCIC_HH
#define OPM_ECL_PROBLEM_BCIC_HH

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <ebos/eclequilinitializer.hh>

#include <array>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Handling of boundary- and initial conditions for EclProblem.
 */
template <class TypeTag>
class EclProblemBCIC
{
public:
    using InitialFluidState = typename EclEquilInitializer<TypeTag>::ScalarFluidState;

    //! \brief Returns whether or not boundary conditions are trivial.
    bool nonTrivialBoundaryConditions() const
    {
        return nonTrivialBoundaryConditions_;
    }

    //! \brief Returns a const reference to initial fluid state for an element.
    const InitialFluidState& initialFluidState(const unsigned idx) const
    {
        return initialFluidStates_[idx];
    }

    //! \brief Container for boundary conditions.
    template<class T>
    struct BCData
    {
        std::array<std::vector<T>,6> data; //!< One vector per FaceDir::DirEnum entry

        //! \brief Resize vectors.
        void resize(std::size_t size, T defVal)
        {
            for (auto& d : data) {
                d.resize(size, defVal);
            }
        }

        //! \brief Returns a const reference to data for a direction.
        const std::vector<T>& operator()(FaceDir::DirEnum dir) const
        {
            if (dir == FaceDir::DirEnum::Unknown) {
                throw std::runtime_error("Tried to access BC data for the 'Unknown' direction");
            }
            int idx = 0;
            int div = static_cast<int>(dir);
            while ((div /= 2) >= 1) {
              ++idx;
            }
            assert(idx >= 0 && idx <= 5);
            return data[idx];
        }

        //! \brief Returns a reference to data for a direction.
        std::vector<T>& operator()(FaceDir::DirEnum dir)
        {
            return const_cast<std::vector<T>&>(std::as_const(*this)(dir));
        }
    };

    BCData<int> bcindex_; //!< Indices for boundary conditions
    bool nonTrivialBoundaryConditions_ = false; //!< Whether or not non-trivial boundary conditions are used
    std::vector<InitialFluidState> initialFluidStates_; //!< Vector of initial fluid states for elements
};

} // namespace Opm

#endif
