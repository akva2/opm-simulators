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

#include <config.h>
#include <opm/simulators/wells/MultisegmentWellEquations.hpp>

#include <dune/istl/umfpack.hh>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/linalg/bda/WellContributions.hpp>
#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/SmallDenseMatrixUtils.hpp>
#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <opm/simulators/wells/MultisegmentWellGeneric.hpp>

namespace Opm {

template<class Indices, class Scalar>
void MultisegmentWellEquations<Indices,Scalar>::
init(const int num_cells,
     const MultisegmentWellGeneric<Scalar>& well) const
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise);
    duneD_.setBuildMode(DiagMatWell::row_wise);

    // set the size and patterns for all the matrices and vectors
    // [A C^T    [x    =  [ res
    //  B D] x_well]      res_well]

    // calculatiing the NNZ for duneD_
    // NNZ = number_of_segments + 2 * (number_of_inlets / number_of_outlets)
    {
        int nnz_d = well.numberOfSegments();
        for (const std::vector<int>& inlets : well.inlets()) {
            nnz_d += 2 * inlets.size();
        }
        duneD_.setSize(well.numberOfSegments(), well.numberOfSegments(), nnz_d);
    }
    duneB_.setSize(well.numberOfSegments(), num_cells, well.numPerfs());
    duneC_.setSize(well.numberOfSegments(), num_cells, well.numPerfs());

    // we need to add the off diagonal ones
    for (auto row = duneD_.createbegin(),
              end = duneD_.createend(); row != end; ++row) {
        // the number of the row corrspnds to the segment now
        const int seg = row.index();
        // adding the item related to outlet relation
        const Segment& segment = well.segmentSet()[seg];
        const int outlet_segment_number = segment.outletSegment();
        if (outlet_segment_number > 0) { // if there is a outlet_segment
            const int outlet_segment_index = well.segmentNumberToIndex(outlet_segment_number);
            row.insert(outlet_segment_index);
        }

        // Add nonzeros for diagonal
        row.insert(seg);

        // insert the item related to its inlets
        for (const int& inlet : well.inlets()[seg]) {
            row.insert(inlet);
        }
    }

    // make the C matrix
    for (auto row = duneC_.createbegin(),
              end = duneC_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : well.segmentPerforations()[row.index()]) {
            const int cell_idx = well.cells()[perf];
            row.insert(cell_idx);
        }
    }

    // make the B^T matrix
    for (auto row = duneB_.createbegin(),
              end = duneB_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : well.segmentPerforations()[row.index()]) {
            const int cell_idx = well.cells()[perf];
            row.insert(cell_idx);
        }
    }

    resWell_.resize(well.numberOfSegments());
}

template<class Indices, class Scalar>
void MultisegmentWellEquations<Indices,Scalar>::clear()
{
    // clear all entries
    duneB_ = 0.0;
    duneC_ = 0.0;

    duneD_ = 0.0;
    resWell_ = 0.0;

    duneDSolver_.reset();
}

template<class Indices, class Scalar>
void MultisegmentWellEquations<Indices,Scalar>::
apply(const BVector& x, BVector& Ax) const
{
    BVectorWell Bx(duneB_.N());

    duneB_.mv(x, Bx);

    // invDBx = duneD^-1 * Bx_
    const BVectorWell invDBx = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, Bx);

    // Ax = Ax - duneC_^T * invDBx
    duneC_.mmtv(invDBx,Ax);
}

template<class Indices, class Scalar>
void MultisegmentWellEquations<Indices,Scalar>::
apply(BVector& r) const
{
    // invDrw_ = duneD^-1 * resWell_
    const BVectorWell invDrw = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell_);
    // r = r - duneC_^T * invDrw
    duneC_.mmtv(invDrw, r);
}

template<class Indices, class Scalar>
typename MultisegmentWellEquations<Indices,Scalar>::BVectorWell
MultisegmentWellEquations<Indices,Scalar>::solve() const
{
    return mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell_);
}

template<class Indices, class Scalar>
void MultisegmentWellEquations<Indices,Scalar>::
recoverSolutionWell(const BVector& x, BVectorWell& xw) const
{
    BVectorWell resWell = resWell_;
    // resWell = resWell - B * x
    duneB_.mmv(x, resWell);
    // xw = D^-1 * resWell
    xw = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell);
}

template<class Indices, class Scalar>
void MultisegmentWellEquations<Indices,Scalar>::
addWellContribution(WellContributions& wellContribs) const
{
    unsigned int Mb = duneB_.N();       // number of blockrows in duneB_, duneC_ and duneD_
    unsigned int BnumBlocks = duneB_.nonzeroes();
    unsigned int DnumBlocks = duneD_.nonzeroes();

    // duneC
    std::vector<unsigned int> Ccols;
    std::vector<double> Cvals;
    Ccols.reserve(BnumBlocks);
    Cvals.reserve(BnumBlocks * Indices::numEq * numWellEq);
    for (auto rowC = duneC_.begin(); rowC != duneC_.end(); ++rowC) {
        for (auto colC = rowC->begin(), endC = rowC->end(); colC != endC; ++colC) {
            Ccols.emplace_back(colC.index());
            for (int i = 0; i < numWellEq; ++i) {
                for (int j = 0; j < Indices::numEq; ++j) {
                    Cvals.emplace_back((*colC)[i][j]);
                }
            }
        }
    }

    // duneD
    Dune::UMFPack<DiagMatWell> umfpackMatrix(duneD_, 0);
    double *Dvals = umfpackMatrix.getInternalMatrix().getValues();
    auto *Dcols = umfpackMatrix.getInternalMatrix().getColStart();
    auto *Drows = umfpackMatrix.getInternalMatrix().getRowIndex();

    // duneB
    std::vector<unsigned int> Bcols;
    std::vector<unsigned int> Brows;
    std::vector<double> Bvals;
    Bcols.reserve(BnumBlocks);
    Brows.reserve(Mb+1);
    Bvals.reserve(BnumBlocks * Indices::numEq * numWellEq);
    Brows.emplace_back(0);
    unsigned int sumBlocks = 0;
    for (auto rowB = duneB_.begin(); rowB != duneB_.end(); ++rowB) {
        int sizeRow = 0;
        for (auto colB = rowB->begin(), endB = rowB->end(); colB != endB; ++colB) {
            Bcols.emplace_back(colB.index());
            for (int i = 0; i < numWellEq; ++i) {
                for (int j = 0; j < Indices::numEq; ++j) {
                    Bvals.emplace_back((*colB)[i][j]);
                }
            }
            sizeRow++;
        }
        sumBlocks += sizeRow;
        Brows.emplace_back(sumBlocks);
    }

    wellContribs.addMultisegmentWellContribution(Indices::numEq,
                                                 numWellEq,
                                                 Mb,
                                                 Bvals,
                                                 Bcols,
                                                 Brows,
                                                 DnumBlocks,
                                                 Dvals,
                                                 Dcols,
                                                 Drows,
                                                 Cvals);
}

template<class Indices, class Scalar>
template<class SparseMatrixAdapter>
void MultisegmentWellEquations<Indices,Scalar>::
addWellContributions(SparseMatrixAdapter& jacobian) const
{
    const auto invDuneD =
        mswellhelpers::invertWithUMFPack<DiagMatWell,BVectorWell>(duneD_, duneDSolver_);

    // We need to change matrix A as follows
    // A -= C^T D^-1 B
    // D is a (nseg x nseg) block matrix with (4 x 4) blocks.
    // B and C are (nseg x ncells) block matrices with (4 x 4 blocks).
    // They have nonzeros at (i, j) only if this well has a
    // perforation at cell j connected to segment i.  The code
    // assumes that no cell is connected to more than one segment,
    // i.e. the columns of B/C have no more than one nonzero.
    for (size_t rowC = 0; rowC < duneC_.N(); ++rowC) {
        for (auto colC = duneC_[rowC].begin(),
                  endC = duneC_[rowC].end(); colC != endC; ++colC) {
            const auto row_index = colC.index();
            for (size_t rowB = 0; rowB < duneB_.N(); ++rowB) {
                for (auto colB = duneB_[rowB].begin(),
                          endB = duneB_[rowB].end(); colB != endB; ++colB) {
                    const auto col_index = colB.index();
                    OffDiagMatrixBlockWellType tmp1;
                    detail::multMatrixImpl(invDuneD[rowC][rowB], (*colB), tmp1, std::true_type());
                    typename SparseMatrixAdapter::MatrixBlock tmp2;
                    detail::multMatrixTransposedImpl((*colC), tmp1, tmp2, std::false_type());
                    jacobian.addToBlock(row_index, col_index, tmp2);
                }
            }
        }
    }
}

#define INSTANCE(Block,...) \
template class MultisegmentWellEquations<__VA_ARGS__,double>; \
template void \
MultisegmentWellEquations<__VA_ARGS__,double>:: \
addWellContributions(Linear::IstlSparseMatrixAdapter<MatrixBlock<double,Block,Block>>&) const;

// One phase
INSTANCE(1, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(2, BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(6, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(1, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(3, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(3, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(4, BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(3, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(3, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)

// Blackoil
INSTANCE(3, BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(3, BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(4, BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(4, BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(3, BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)
INSTANCE(4, BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(4, BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(4, BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(4, BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(5, BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

}
