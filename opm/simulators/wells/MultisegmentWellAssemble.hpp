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


#ifndef OPM_MULTISEGMENTWELL_ASSEMBLE_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_ASSEMBLE_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <functional>
#include <vector>

namespace Opm
{

class DeferredLogger;
class GroupState;
template<class Indices, class Sclar> class MultisegmentWellEquations;
class Schedule;
class SummaryState;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class MultisegmentWellAssemble {
public:
    MultisegmentWellAssemble(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well)
        : well_(well)
    {}

    template<class EvalWell>
    void assembleControlEq(const WellState& well_state,
                           const GroupState& group_state,
                           const Schedule& schedule,
                           const SummaryState& summaryState,
                           const Well::InjectionControls& inj_controls,
                           const Well::ProductionControls& prod_controls,
                           const EvalWell& wqTotal,
                           const EvalWell& bhp,
                           const double rho,
                           const int SPres,
                           const std::function<EvalWell(int)>& getQs,
                           MultisegmentWellEquations<Indices,Scalar>& eqns,
                           DeferredLogger& deferred_logger) const;

    template<class EvalWell>
    void assemblePressureEq(const int seg,
                            const int seg_upwind,
                            const EvalWell& pressure_equation,
                            const EvalWell& outlet_pressure,
                            const int WFrac,
                            const int GFrac,
                            const int SPres,
                            const int WQTotal,
                            const int outlet_segment_index,
                            MultisegmentWellEquations<Indices,Scalar>& eqns) const;

    void assembleTrivialEq(const int seg,
                           Scalar value,
                           const int SPres,
                           const int WQTotal,
                           MultisegmentWellEquations<Indices,Scalar>& eqns) const;

    template<class EvalWell>
    void assembleAccelerationTerm(const int seg,
                                  const int comp_idx,
                                  const EvalWell& accumulation_term,
                                  MultisegmentWellEquations<Indices,Scalar>& eqns) const;

    template<class EvalWell>
    void assembleFlowTerm(const int seg,
                          const int seg_upwind,
                          const int comp_idx,
                          const int WFrac,
                          const int GFrac,
                          const int WQTotal,
                          const EvalWell& segment_rate,
                          MultisegmentWellEquations<Indices,Scalar>& eqns) const;

    template<class EvalWell>
    void assemblePerfTerm(const int seg,
                          const int cell_idx,
                          const int comp_idx,
                          const EvalWell& cq_s_effective,
                          MultisegmentWellEquations<Indices,Scalar>& eqns) const;

private:
    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_;
};

}

#endif // OPM_MULTISEGMENTWELL_ASSEMBLE_HEADER_INCLUDED
