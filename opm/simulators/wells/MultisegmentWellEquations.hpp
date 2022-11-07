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


#ifndef OPM_MULTISEGMENTWELL_EQUATIONS_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_EQUATIONS_HEADER_INCLUDED

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <memory>

namespace Dune {
template<class M> class UMFPack;
}

namespace Opm
{

template<class T> class MultisegmentWellGeneric;
class WellContributions;

//! \brief Matrices and vectors for the well.
template<class Indices, class Scalar>
class MultisegmentWellEquations {
public:
    //  the number of well equations  TODO: it should have a more general strategy for it
    static constexpr int numWellEq = Indices::numPhases + 1;

    // sparsity pattern for the matrices
    // [A C^T    [x       =  [ res
    //  B  D ]   x_well]      res_well]

    // the vector type for the res_well and x_well
    using VectorBlockWellType = Dune::FieldVector<Scalar, numWellEq>;
    using BVectorWell = Dune::BlockVector<VectorBlockWellType>;

    using VectorBlockType = Dune::FieldVector<Scalar, Indices::numEq>;
    using BVector = Dune::BlockVector<VectorBlockType>;

    // the matrix type for the diagonal matrix D
    using DiagMatrixBlockWellType = Dune::FieldMatrix<Scalar, numWellEq, numWellEq>;
    using DiagMatWell = Dune::BCRSMatrix<DiagMatrixBlockWellType>;

    // the matrix type for the non-diagonal matrix B and C^T
    using OffDiagMatrixBlockWellType = Dune::FieldMatrix<Scalar, numWellEq, Indices::numEq>;
    using OffDiagMatWell = Dune::BCRSMatrix<OffDiagMatrixBlockWellType>;

    //! \brief Initialize matrices and vectors.
    void init(const int num_cells,
              const MultisegmentWellGeneric<Scalar>& well) const;

    //! \brief Clear equation system.
    void clear();

    //! \brief Ax = Ax - C D^-1 B x
    void apply(const BVector& x, BVector& Ax) const;

    //! \brief r = r - C D^-1 Rw
    void apply(BVector& r) const;

    //! \brief dx_well = D^-1*dx_well
    BVectorWell solve() const;

    // xw = inv(D)*(rw - C*x)
    void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

    //! \brief Add the contribution (C, D, B matrices) of this Well to the WellContributions object
    void addWellContribution(WellContributions& wellContribs) const;

    template<class SparseMatrixAdapter>
    void addWellContributions(SparseMatrixAdapter& jacobian) const;

    const BVectorWell& residual() const { return resWell_; }

    // TODO, the following should go to a class for computing purpose
    // two off-diagonal matrices
    mutable OffDiagMatWell duneB_;
    mutable OffDiagMatWell duneC_;
    // "diagonal" matrix for the well. It has offdiagonal entries for inlets and outlets.
    mutable DiagMatWell duneD_;

    // residuals of the well equations
    mutable BVectorWell resWell_;

protected:
    /// \brief solver for diagonal matrix
    ///
    /// This is a shared_ptr as MultisegmentWell is copied in computeWellPotentials...
    mutable std::shared_ptr<Dune::UMFPack<DiagMatWell> > duneDSolver_;
};

}

#endif // OPM_MULTISEGMENTWELL_EQUATIONS_HEADER_INCLUDED
