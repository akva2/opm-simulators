/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELLINTERFACE_HEADER_INCLUDED
#define OPM_WELLINTERFACE_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
// NOTE: GasLiftSingleWell.hpp includes StandardWell.hpp which includes ourself
//   (WellInterface.hpp), so we need to forward declare GasLiftSingleWell
//   for it to be defined in this file. Similar for BlackoilWellModel
namespace Opm {
    template<typename TypeTag> class GasLiftSingleWell;
    template<typename TypeTag> class BlackoilWellModel;
}
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <string>
#include <memory>
#include <vector>
#include <cassert>

namespace Opm
{


    template<typename TypeTag>
    class WellInterface : public WellInterfaceGeneric
    {
    public:

        using WellState = WellStateFullyImplicitBlackoil;

        using ModelParameters = BlackoilModelParametersEbos<TypeTag>;

        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;

        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
        using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using RateVector = GetPropType<TypeTag, Properties::RateVector>;
        using GasLiftSingleWell = ::Opm::GasLiftSingleWell<TypeTag>;
        using GLiftOptWells = typename BlackoilWellModel<TypeTag>::GLiftOptWells;
        using GLiftProdWells = typename BlackoilWellModel<TypeTag>::GLiftProdWells;
        using GLiftWellStateMap =
            typename BlackoilWellModel<TypeTag>::GLiftWellStateMap;

        static const int numEq = Indices::numEq;
        static const int numPhases = Indices::numPhases;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;

        using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
        using MatrixBlockType = Dune::FieldMatrix<Scalar, numEq, numEq>;
        using BVector = Dune::BlockVector<VectorBlockType>;
        using Eval = DenseAd::Evaluation<Scalar, /*size=*/numEq>;

        static constexpr bool has_solvent = getPropValue<TypeTag, Properties::EnableSolvent>();
        static constexpr bool has_zFraction = getPropValue<TypeTag, Properties::EnableExtbo>();
        static constexpr bool has_polymer = getPropValue<TypeTag, Properties::EnablePolymer>();
        static constexpr bool has_energy = getPropValue<TypeTag, Properties::EnableEnergy>();
        static const bool has_temperature = getPropValue<TypeTag, Properties::EnableTemperature>();
        // flag for polymer molecular weight related
        static constexpr bool has_polymermw = getPropValue<TypeTag, Properties::EnablePolymerMW>();
        static constexpr bool has_foam = getPropValue<TypeTag, Properties::EnableFoam>();
        static constexpr bool has_brine = getPropValue<TypeTag, Properties::EnableBrine>();
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiZfracEqIdx = Indices::contiZfracEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        // index for the polymer molecular weight continuity equation
        static const int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
        static const int contiFoamEqIdx = Indices::contiFoamEqIdx;
        static const int contiBrineEqIdx = Indices::contiBrineEqIdx;

        // For the conversion between the surface volume rate and reservoir voidage rate
        using RateConverterType = RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
        static const bool compositionSwitchEnabled = Indices::gasEnabled;
        using FluidState = BlackOilFluidState<Eval,
                                              FluidSystem,
                                              has_temperature,
                                              has_energy,
                                              compositionSwitchEnabled,
                                              has_brine,
                                              Indices::numPhases >;
        /// Constructor
        WellInterface(const Well& well,
                      const ParallelWellInfo& pw_info,
                      const int time_step,
                      const ModelParameters& param,
                      const RateConverterType& rate_converter,
                      const int pvtRegionIdx,
                      const int num_components,
                      const int num_phases,
                      const int index_of_well,
                      const int first_perf_index,
                      const std::vector<PerforationData>& perf_data);

        /// Virtual destructor
        virtual ~WellInterface() = default;

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells,
                          const std::vector< Scalar >& B_avg);

        virtual void initPrimaryVariablesEvaluation() const = 0;

        virtual ConvergenceReport getWellConvergence(const WellState& well_state, const std::vector<double>& B_avg, DeferredLogger& deferred_logger, const bool relax_tolerance = false) const = 0;

        virtual void solveEqAndUpdateWellState(WellState& well_state, DeferredLogger& deferred_logger) = 0;

        void assembleWellEq(const Simulator& ebosSimulator,
                            const double dt,
                            WellState& well_state,
                            const GroupState& group_state,
                            DeferredLogger& deferred_logger);

        virtual void gasLiftOptimizationStage1 (
            WellState& well_state,
            const Simulator& ebosSimulator,
            DeferredLogger& deferred_logger,
            GLiftProdWells& prod_wells,
            GLiftOptWells& glift_wells,
            GLiftWellStateMap& state_map
        ) const = 0;

        void updateWellTestState(const WellState& well_state,
                                 const double& simulationTime,
                                 const bool& writeMessageToOPMLog,
                                 WellTestState& wellTestState,
                                 DeferredLogger& deferred_logger) const;

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state,
                                                           DeferredLogger& deferred_logger) const = 0;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const = 0;

        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const = 0;

        virtual std::optional<double> computeBhpAtThpLimitProdWithAlq(const Simulator&,
                                                                      const SummaryState&,
                                                                      DeferredLogger&,
                                                                      double) const
        {
            return {};
        }

        // TODO: before we decide to put more information under mutable, this function is not const
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials,
                                           DeferredLogger& deferred_logger) = 0;

        virtual void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                             const double& bhp,
                                             std::vector<double>& well_flux,
                                             Opm::DeferredLogger& deferred_logger) const = 0;

        void updateWellStateWithTarget(const Simulator& ebos_simulator,
                                       WellState& well_state,
                                       DeferredLogger& deferred_logger) const;

        enum class IndividualOrGroup { Individual, Group, Both };
        bool updateWellControl(const Simulator& ebos_simulator,
                               const IndividualOrGroup iog,
                               WellState& well_state,
                               const GroupState& group_state,
                               DeferredLogger& deferred_logger) /* const */;

        virtual void updatePrimaryVariables(const WellState& well_state, DeferredLogger& deferred_logger) const = 0;

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state,
                                                 DeferredLogger& deferred_logger) = 0; // should be const?

        virtual void updateProductivityIndex(const Simulator& ebosSimulator,
                                             const WellProdIndexCalculator& wellPICalc,
                                             WellState& well_state,
                                             DeferredLogger& deferred_logger) const = 0;

        /// \brief Wether the Jacobian will also have well contributions in it.
        virtual bool jacobianContainsWellContributions() const
        {
            return false;
        }

        // updating the voidage rates in well_state when requested
        void calculateReservoirRates(WellState& well_state) const;

        // Add well contributions to matrix
        virtual void addWellContributions(SparseMatrixAdapter&) const = 0;

        void addCellRates(RateVector& rates, int cellIdx) const;

        Scalar volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const;


        template <class EvalWell>
        Eval restrictEval(const EvalWell& in) const
        {
            Eval out = 0.0;
            out.setValue(in.value());
            for(int eqIdx = 0; eqIdx < numEq;++eqIdx) {
                out.setDerivative(eqIdx, in.derivative(eqIdx));
            }
            return out;
        }

        // TODO: theoretically, it should be a const function
        // Simulator is not const is because that assembleWellEq is non-const Simulator
        void wellTesting(const Simulator& simulator,
                         const double simulation_time, const int report_step,
                         const WellTestConfig::Reason testing_reason,
                         /* const */ WellState& well_state, const GroupState& group_state, WellTestState& welltest_state,
                         DeferredLogger& deferred_logger);

        void checkWellOperability(const Simulator& ebos_simulator, const WellState& well_state, DeferredLogger& deferred_logger);

        // check whether the well is operable under the current reservoir condition
        // mostly related to BHP limit and THP limit
        void updateWellOperability(const Simulator& ebos_simulator,
                                   const WellState& well_state,
                                   DeferredLogger& deferred_logger);


        // update perforation water throughput based on solved water rate
        virtual void updateWaterThroughput(const double dt, WellState& well_state) const = 0;

        /// Compute well rates based on current reservoir conditions and well variables.
        /// Used in updateWellStateRates().
        virtual std::vector<double> computeCurrentWellRates(const Simulator& ebosSimulator,
                                                            DeferredLogger& deferred_logger) const = 0;

        /// Modify the well_state's rates if there is only one nonzero rate.
        /// If so, that rate is kept as is, but the others are set proportionally
        /// to the rates returned by computeCurrentWellRates().
        void updateWellStateRates(const Simulator& ebosSimulator,
                                  WellState& well_state,
                                  DeferredLogger& deferred_logger) const;

        void solveWellEquation(const Simulator& ebosSimulator,
                               WellState& well_state,
                               const GroupState& group_state,
                               DeferredLogger& deferred_logger);

        virtual bool useInnerIterations() const = 0;

    protected:

        // to indicate a invalid completion
        static const int INVALIDCOMPLETION = INT_MAX;

        // simulation parameters
        const ModelParameters& param_;

        // For the conversion between the surface volume rate and resrevoir voidage rate
        const RateConverterType& rateConverter_;

        std::vector<RateVector> connectionRates_;

        std::vector< Scalar > B_avg_;

        bool changed_to_stopped_this_step_ = false;

        int flowPhaseToEbosCompIdx( const int phaseIdx ) const;

        int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const;

        int ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const;

        double wpolymer() const;

        double wfoam() const;

        double wsalt() const;

        bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                 const double * rates_or_potentials,
                                 DeferredLogger& deferred_logger) const;

        template <class ValueType>
        ValueType calculateBhpFromThp(const WellState& well_state, const std::vector<ValueType>& rates, const Well& well, const SummaryState& summaryState, DeferredLogger& deferred_logger) const;

        virtual double getRefDensity() const = 0;

        // Component fractions for each phase for the well
        const std::vector<double>& compFrac() const;

        struct RatioLimitCheckReport;

        void checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                  const WellState& well_state,
                                  RatioLimitCheckReport& report,
                                  DeferredLogger& deferred_logger) const;


        template <typename RatioFunc>
        bool checkMaxRatioLimitWell(const WellState& well_state,
                                    const double max_ratio_limit,
                                    const RatioFunc& ratioFunc) const;

        template <typename RatioFunc>
        void checkMaxRatioLimitCompletions(const WellState& well_state,
                                           const double max_ratio_limit,
                                           const RatioFunc& ratioFunc,
                                           RatioLimitCheckReport& report) const;

        double scalingFactor(const int comp_idx) const;

        std::vector<double> initialWellRateFractions(const Simulator& ebosSimulator, const std::vector<double>& potentials) const;

        // check whether the well is operable under BHP limit with current reservoir condition
        virtual void checkOperabilityUnderBHPLimitProducer(const WellState& well_state, const Simulator& ebos_simulator, DeferredLogger& deferred_logger) =0;

        // check whether the well is operable under THP limit with current reservoir condition
        virtual void checkOperabilityUnderTHPLimitProducer(const Simulator& ebos_simulator, const WellState& well_state, DeferredLogger& deferred_logger) =0;

        virtual void updateIPR(const Simulator& ebos_simulator, DeferredLogger& deferred_logger) const=0;


        void wellTestingEconomic(const Simulator& simulator,
                                 const double simulation_time, const WellState& well_state, const GroupState& group_state,
                                 WellTestState& welltest_state, DeferredLogger& deferred_logger);

        void wellTestingPhysical(const Simulator& simulator,
                                 const double simulation_time, const int report_step,
                                 WellState& well_state,
                                 const GroupState& group_state,
                                 WellTestState& welltest_state, DeferredLogger& deferred_logger);


        virtual void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                                    const double dt,
                                                    const Well::InjectionControls& inj_controls,
                                                    const Well::ProductionControls& prod_controls,
                                                    WellState& well_state,
                                                    const GroupState& group_state,
                                                    DeferredLogger& deferred_logger) = 0;

        // iterate well equations with the specified control until converged
        virtual bool iterateWellEqWithControl(const Simulator& ebosSimulator,
                                              const double dt,
                                              const Well::InjectionControls& inj_controls,
                                              const Well::ProductionControls& prod_controls,
                                              WellState& well_state,
                                              const GroupState& group_state,
                                              DeferredLogger& deferred_logger) = 0;

        bool iterateWellEquations(const Simulator& ebosSimulator,
                                  const double dt,
                                  WellState& well_state,
                                  const GroupState& group_state,
                                  DeferredLogger& deferred_logger);

        void updateWellTestStateEconomic(const WellState& well_state,
                                         const double simulation_time,
                                         const bool write_message_to_opmlog,
                                         WellTestState& well_test_state,
                                         DeferredLogger& deferred_logger) const;

        void solveWellForTesting(const Simulator& ebosSimulator, WellState& well_state, const GroupState& group_state,
                                 DeferredLogger& deferred_logger);

        bool checkConstraints(WellState& well_state,
                              const GroupState& group_state,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              DeferredLogger& deferred_logger) const;

        bool checkIndividualConstraints(WellState& well_state,
                                        const SummaryState& summaryState) const;

        bool checkGroupConstraints(WellState& well_state,
                                   const GroupState& group_state,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   DeferredLogger& deferred_logger) const;

        std::pair<bool, double> checkGroupConstraintsProd(const Group& group,
                                       const WellState& well_state,
                                                          const GroupState& group_state,
                                       const double efficiencyFactor,
                                       const Schedule& schedule,
                                       const SummaryState& summaryState,
                                       DeferredLogger& deferred_logger) const;

        std::pair<bool, double> checkGroupConstraintsInj(const Group& group,
                                      const WellState& well_state,
                                                         const GroupState& group_state,
                                      const double efficiencyFactor,
                                      const Schedule& schedule,
                                      const SummaryState& summaryState,
                                      DeferredLogger& deferred_logger) const;

        template <class EvalWell>
        void getGroupInjectionControl(const Group& group,
                                      const WellState& well_state,
                                      const GroupState& group_state,
                                      const Schedule& schedule,
                                      const SummaryState& summaryState,
                                      const InjectorType& injectorType,
                                      const EvalWell& bhp,
                                      const EvalWell& injection_rate,
                                      EvalWell& control_eq,
                                      double efficiencyFactor,
                                      DeferredLogger& deferred_logger);

        template <class EvalWell>
        void getGroupProductionControl(const Group& group,
                                       const WellState& well_state,
                                       const GroupState& group_state,
                                       const Schedule& schedule,
                                       const SummaryState& summaryState,
                                       const EvalWell& bhp,
                                       const std::vector<EvalWell>& rates,
                                       EvalWell& control_eq,
                                       double efficiencyFactor);

        template <class EvalWell, class BhpFromThpFunc>
        void assembleControlEqInj(const WellState& well_state,
                                  const GroupState& group_state,
                                  const Schedule& schedule,
                                  const SummaryState& summaryState,
                                  const Well::InjectionControls& controls,
                                  const EvalWell& bhp,
                                  const EvalWell& injection_rate,
                                  BhpFromThpFunc bhp_from_thp,
                                  EvalWell& control_eq,
                                  DeferredLogger& deferred_logger);

        template <class EvalWell, class BhpFromThpFunc>
        void assembleControlEqProd(const WellState& well_state,
                                   const GroupState& group_state,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   const Well::ProductionControls& controls,
                                   const EvalWell& bhp,
                                   const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                                   BhpFromThpFunc bhp_from_thp,
                                   EvalWell& control_eq,
                                   DeferredLogger& deferred_logger);
    };

    template<typename TypeTag>
    struct
    WellInterface<TypeTag>::
    RatioLimitCheckReport{
        bool ratio_limit_violated = false;
        int worst_offending_completion = INVALIDCOMPLETION;
        double violation_extent = 0.0;
    };

}

#include "WellInterface_impl.hpp"

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
