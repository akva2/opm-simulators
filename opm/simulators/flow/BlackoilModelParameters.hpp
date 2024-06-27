/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED
#define OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/SubDomain.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace Opm::Properties {

namespace TTag {
struct FlowModelParameters {};
}

struct EclDeckFileName {
    static constexpr auto value = "";
};
// template<class TypeTag, class MyTypeTag>
// struct DbhpMaxRel {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct DwellFractionMax {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxResidualAllowed {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct RelaxedMaxPvFraction {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ToleranceMb {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ToleranceMbRelaxed {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ToleranceCnv {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ToleranceCnvRelaxed {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ToleranceWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ToleranceWellControl {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxWelleqIter {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct UseMultisegmentWell {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxSinglePrecisionDays {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MinStrictCnvIter {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MinStrictMbIter {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct SolveWelleqInitially {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct UpdateEquationsScaling {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct UseUpdateStabilization {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MatrixAddWellContributions {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct EnableWellOperabilityCheck {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct EnableWellOperabilityCheckIter {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct DebugEmitCellPartition {
//     using type = UndefinedProperty;
// };
// parameters for multisegment wells
// template<class TypeTag, class MyTypeTag>
// struct TolerancePressureMsWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxPressureChangeMsWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxInnerIterMsWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct StrictInnerIterWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct RelaxedWellFlowTol {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct StrictOuterIterWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct RelaxedPressureTolMsw {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct RegularizationFactorWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxNewtonIterationsWithInnerWellIterations  {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ShutUnsolvableWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxInnerIterWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct AlternativeWellRateInit {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaximumNumberOfWellSwitches {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct UseAverageDensityMsWells {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LocalWellSolveControlSwitching {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct UseImplicitIpr {
//     using type = UndefinedProperty;
// };
// Network solver parameters
// template<class TypeTag, class MyTypeTag>
// struct NetworkMaxStrictIterations {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct NetworkMaxIterations {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct NonlinearSolver {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LocalSolveApproach {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MaxLocalSolveIterations {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LocalToleranceScalingMb {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LocalToleranceScalingCnv {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct NlddNumInitialNewtonIter {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct NumLocalDomains {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LocalDomainsPartitioningImbalance {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LocalDomainsPartitioningMethod {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LocalDomainsOrderingMeasure {
//     using type = UndefinedProperty;
// };
struct DbhpMaxRel {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
struct DwellFractionMax
{
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.2;
};
struct MaxResidualAllowed {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e7;
};
struct RelaxedMaxPvFraction {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.03;
};
struct ToleranceMb {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-6;
};
struct ToleranceMbRelaxed {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-6;
};
struct ToleranceCnv {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};
struct ToleranceCnvRelaxed {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1;
};
struct ToleranceWells {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-4;
};
struct ToleranceWellControl {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-7;
};
struct MaxWelleqIter {
    static constexpr int value = 30;
};
struct UseMultisegmentWell {
    static constexpr bool value = true;
};
struct MaxSinglePrecisionDays {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 20.0;
};
struct MinStrictCnvIter {
    static constexpr int value = 0;
};
struct MinStrictMbIter {
    static constexpr int value = -1;
};
struct SolveWelleqInitially {
    static constexpr bool value = true;
};
struct UpdateEquationsScaling { static constexpr bool value = false; };
struct UseUpdateStabilization {
    static constexpr bool value = true;
};
struct MatrixAddWellContributions {
    static constexpr bool value = false;
};
struct TolerancePressureMsWells {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.01*1e5;
};
struct MaxPressureChangeMsWells {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10*1e5;
};
struct MaxNewtonIterationsWithInnerWellIterations {
    static constexpr int value = 8;
};
struct MaxInnerIterMsWells
{
    static constexpr int value = 100;
};
struct MaxInnerIterWells
{
    static constexpr int value = 50;
};
struct ShutUnsolvableWells {
    static constexpr bool value = true;
};
struct AlternativeWellRateInit {
    static constexpr bool value = true;
};
struct StrictOuterIterWells {
    static constexpr int value = 6;
};
struct StrictInnerIterWells {
    static constexpr int value = 40;
};
struct RegularizationFactorWells {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100;
};
struct EnableWellOperabilityCheck {
    static constexpr bool value = true;
};
struct EnableWellOperabilityCheckIter {
    static constexpr bool value = false;
};
struct DebugEmitCellPartition {
    static constexpr bool value = false;
};
struct RelaxedWellFlowTol {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};
struct RelaxedPressureTolMsw {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0e4;
};
struct MaximumNumberOfWellSwitches {
    static constexpr int value = 3;
};
struct UseAverageDensityMsWells {
    static constexpr bool value = false;
};
struct LocalWellSolveControlSwitching {
    static constexpr bool value = false;
};
struct UseImplicitIpr {
    static constexpr bool value = false;
};

// Network solver parameters
struct NetworkMaxStrictIterations {
    static constexpr int value = 100;
};
struct NetworkMaxIterations {
    static constexpr int value = 200;
};
struct NonlinearSolver {
    static constexpr auto value = "newton";
};
struct LocalSolveApproach {
    static constexpr auto value = "gauss-seidel";
};
struct MaxLocalSolveIterations {
    static constexpr int value = 20;
};
struct LocalToleranceScalingMb {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
struct LocalToleranceScalingCnv {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.1;
};
struct NlddNumInitialNewtonIter {
    using type = int;
    static constexpr auto value = type{1};
};
struct NumLocalDomains {
    using type = int;
    static constexpr auto value = 0;
};
struct LocalDomainsPartitioningImbalance {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr auto value = type{1.03};
};
struct LocalDomainsPartitioningMethod {
    static constexpr auto value = "zoltan";
};
struct LocalDomainsOrderingMeasure {
    static constexpr auto value = "maxpressure";
};
// if openMP is available, determine the number threads per process automatically.
#if _OPENMP
// template<class TypeTag>
// struct ThreadsPerProcess<TypeTag, TTag::FlowModelParameters> {
//     static constexpr int value = -1;
// };
#endif

} // namespace Opm::Properties

namespace Opm
{
    /// Solver parameters for the BlackoilModel.
    template <class TypeTag>
    struct BlackoilModelParameters
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    public:
        /// Max relative change in bhp in single iteration.
        Scalar dbhp_max_rel_;
        /// Max absolute change in well volume fraction in single iteration.
        Scalar dwell_fraction_max_;
        /// Absolute max limit for residuals.
        Scalar max_residual_allowed_;
        //// Max allowed pore volume faction where CNV is violated. Below the
        //// relaxed tolerance tolerance_cnv_relaxed_ is used.
        Scalar relaxed_max_pv_fraction_;
        /// Relative mass balance tolerance (total mass balance error).
        Scalar tolerance_mb_;
        /// Relaxed mass balance tolerance (can be used when iter >= min_strict_mb_iter_).
        Scalar tolerance_mb_relaxed_;
        /// Local convergence tolerance (max of local saturation errors).
        Scalar tolerance_cnv_;
        /// Relaxed local convergence tolerance (can be used when iter >= min_strict_cnv_iter_ && cnvViolatedPV < relaxed_max_pv_fraction_).
        Scalar tolerance_cnv_relaxed_;
        /// Well convergence tolerance.
        Scalar tolerance_wells_;
        /// Tolerance for the well control equations
        //  TODO: it might need to distinguish between rate control and pressure control later
        Scalar tolerance_well_control_;
        /// Tolerance for the pressure equations for multisegment wells
        Scalar tolerance_pressure_ms_wells_;
        /// Relaxed tolerance for for the well flow residual
        Scalar relaxed_tolerance_flow_well_;

        /// Relaxed tolerance for the MSW pressure solution
        Scalar relaxed_tolerance_pressure_ms_well_;

        /// Maximum pressure change over an iteratio for ms wells
        Scalar max_pressure_change_ms_wells_;

        /// Maximum inner iteration number for ms wells
        int max_inner_iter_ms_wells_;

        /// Strict inner iteration number for wells
        int strict_inner_iter_wells_;

        /// Newton iteration where wells are stricly convergent
        int strict_outer_iter_wells_;

        /// Regularization factor for wells
        Scalar regularization_factor_wells_;

        /// Maximum newton iterations with inner well iterations
        int max_niter_inner_well_iter_;

        /// Whether to shut unsolvable well
        bool shut_unsolvable_wells_;

        /// Maximum inner iteration number for standard wells
        int max_inner_iter_wells_;

        /// Maximum iteration number of the well equation solution
        int max_welleq_iter_;

        /// Tolerance for time step in seconds where single precision can be used
        /// for solving for the Jacobian
        Scalar maxSinglePrecisionTimeStep_;

        /// Minimum number of Newton iterations before we can use relaxed CNV convergence criterion
        int min_strict_cnv_iter_;

        /// Minimum number of Newton iterations before we can use relaxed MB convergence criterion
        int min_strict_mb_iter_;

        /// Solve well equation initially
        bool solve_welleq_initially_;

        /// Update scaling factors for mass balance equations
        bool update_equations_scaling_;

        /// Try to detect oscillation or stagnation.
        bool use_update_stabilization_;

        /// Whether to use MultisegmentWell to handle multisegment wells
        /// it is something temporary before the multisegment well model is considered to be
        /// well developed and tested.
        /// if it is false, we will handle multisegment wells as standard wells, which will be
        /// the default behavoir for the moment. Later, we might set it to be true by default if necessary
        bool use_multisegment_well_;

        /// The file name of the deck
        std::string deck_file_name_;

        /// Whether to add influences of wells between cells to the matrix and preconditioner matrix
        bool matrix_add_well_contributions_;

        /// Whether to check well operability
        bool check_well_operability_;
        /// Whether to check well operability during iterations
        bool check_well_operability_iter_;

        /// Maximum number of times a well can switch to the same controt
        int max_number_of_well_switches_;

        /// Whether to approximate segment densities by averaging over segment and its outlet 
        bool use_average_density_ms_wells_;

        /// Whether to allow control switching during local well solutions 
        bool local_well_solver_control_switching_;

        /// Whether to use implicit IPR for thp stability checks and solution search
        bool use_implicit_ipr_;

        /// Maximum number of iterations in the network solver before relaxing tolerance
        int network_max_strict_iterations_;
        
        /// Maximum number of iterations in the network solver before giving up
        int network_max_iterations_;

        /// Nonlinear solver type: newton or nldd.
        std::string nonlinear_solver_;
        /// 'jacobi' and 'gauss-seidel' supported.
        DomainSolveApproach local_solve_approach_{DomainSolveApproach::Jacobi};

        int max_local_solve_iterations_;

        Scalar local_tolerance_scaling_mb_;
        Scalar local_tolerance_scaling_cnv_;

        int nldd_num_initial_newton_iter_{1};
        int num_local_domains_{0};
        Scalar local_domain_partition_imbalance_{1.03};
        std::string local_domain_partition_method_;
        DomainOrderingMeasure local_domain_ordering_{DomainOrderingMeasure::MaxPressure};

        bool write_partitions_{false};

        /// Construct from user parameters or defaults.
        BlackoilModelParameters()
        {
            dbhp_max_rel_=  Parameters::get<Properties::DbhpMaxRel>();
            dwell_fraction_max_ = Parameters::get<Properties::DwellFractionMax>();
            max_residual_allowed_ = Parameters::get<Properties::MaxResidualAllowed>();
            relaxed_max_pv_fraction_ = Parameters::get<Properties::RelaxedMaxPvFraction>();
            tolerance_mb_ = Parameters::get<Properties::ToleranceMb>();
            tolerance_mb_relaxed_ = std::max(tolerance_mb_, Parameters::get<Properties::ToleranceMbRelaxed>());
            tolerance_cnv_ = Parameters::get<Properties::ToleranceCnv>();
            tolerance_cnv_relaxed_ = std::max(tolerance_cnv_, Parameters::get<Properties::ToleranceCnvRelaxed>());
            tolerance_wells_ = Parameters::get<Properties::ToleranceWells>();
            tolerance_well_control_ = Parameters::get<Properties::ToleranceWellControl>();
            max_welleq_iter_ = Parameters::get<Properties::MaxWelleqIter>();
            use_multisegment_well_ = Parameters::get<Properties::UseMultisegmentWell>();
            tolerance_pressure_ms_wells_ = Parameters::get<Properties::TolerancePressureMsWells>();
            relaxed_tolerance_flow_well_ = Parameters::get<Properties::RelaxedWellFlowTol>();
            relaxed_tolerance_pressure_ms_well_ = Parameters::get<Properties::RelaxedPressureTolMsw>();
            max_pressure_change_ms_wells_ = Parameters::get<Properties::MaxPressureChangeMsWells>();
            max_inner_iter_ms_wells_ = Parameters::get<Properties::MaxInnerIterMsWells>();
            strict_inner_iter_wells_ = Parameters::get<Properties::StrictInnerIterWells>();
            strict_outer_iter_wells_ = Parameters::get<Properties::StrictOuterIterWells>();
            regularization_factor_wells_ = Parameters::get<Properties::RegularizationFactorWells>();
            max_niter_inner_well_iter_ = Parameters::get<Properties::MaxNewtonIterationsWithInnerWellIterations>();
            shut_unsolvable_wells_ = Parameters::get<Properties::ShutUnsolvableWells>();
            max_inner_iter_wells_ = Parameters::get<Properties::MaxInnerIterWells>();
            maxSinglePrecisionTimeStep_ = Parameters::get<Properties::MaxSinglePrecisionDays>() * 24 * 60 * 60;
            min_strict_cnv_iter_ = Parameters::get<Properties::MinStrictCnvIter>();
            min_strict_mb_iter_ = Parameters::get<Properties::MinStrictMbIter>();
            solve_welleq_initially_ = Parameters::get<Properties::SolveWelleqInitially>();
            update_equations_scaling_ = Parameters::get<Properties::UpdateEquationsScaling>();
            use_update_stabilization_ = Parameters::get<Properties::UseUpdateStabilization>();
            matrix_add_well_contributions_ = Parameters::get<Properties::MatrixAddWellContributions>();
            check_well_operability_ = Parameters::get<Properties::EnableWellOperabilityCheck>();
            check_well_operability_iter_ = Parameters::get<Properties::EnableWellOperabilityCheckIter>();
            max_number_of_well_switches_ = Parameters::get<Properties::MaximumNumberOfWellSwitches>();
            use_average_density_ms_wells_ = Parameters::get<Properties::UseAverageDensityMsWells>();
            local_well_solver_control_switching_ = Parameters::get<Properties::LocalWellSolveControlSwitching>();
            use_implicit_ipr_ = Parameters::get<Properties::UseImplicitIpr>();
            nonlinear_solver_ = Parameters::get<Properties::NonlinearSolver>();
            const auto approach = Parameters::get<Properties::LocalSolveApproach>();
            if (approach == "jacobi") {
                local_solve_approach_ = DomainSolveApproach::Jacobi;
            } else if (approach == "gauss-seidel") {
                local_solve_approach_ = DomainSolveApproach::GaussSeidel;
            } else {
                throw std::runtime_error("Invalid domain solver approach '" + approach + "' specified.");
            }

            max_local_solve_iterations_ = Parameters::get<Properties::MaxLocalSolveIterations>();
            local_tolerance_scaling_mb_ = Parameters::get<Properties::LocalToleranceScalingMb>();
            local_tolerance_scaling_cnv_ = Parameters::get<Properties::LocalToleranceScalingCnv>();
            nldd_num_initial_newton_iter_ = Parameters::get<Properties::NlddNumInitialNewtonIter>();
            num_local_domains_ = Parameters::get<Properties::NumLocalDomains>();
            local_domain_partition_imbalance_ = std::max(Scalar{1.0}, Parameters::get<Properties::LocalDomainsPartitioningImbalance>());
            local_domain_partition_method_ = Parameters::get<Properties::LocalDomainsPartitioningMethod>();
            deck_file_name_ = Parameters::get<Properties::EclDeckFileName>();
            network_max_strict_iterations_ = Parameters::get<Properties::NetworkMaxStrictIterations>();
            network_max_iterations_ = Parameters::get<Properties::NetworkMaxIterations>();
            local_domain_ordering_ = domainOrderingMeasureFromString(Parameters::get<Properties::LocalDomainsOrderingMeasure>());
            write_partitions_ = Parameters::get<Properties::DebugEmitCellPartition>();
        }

        static void registerParameters()
        {
            Parameters::registerParam<Properties::DbhpMaxRel>
                ("Maximum relative change of the bottom-hole pressure in a single iteration");
            Parameters::registerParam<Properties::DwellFractionMax>
                ("Maximum absolute change of a well's volume fraction in a single iteration");
            Parameters::registerParam<Properties::MaxResidualAllowed>
                ("Absolute maximum tolerated for residuals without cutting the time step size");
            Parameters::registerParam<Properties::RelaxedMaxPvFraction>
                ("The fraction of the pore volume of the reservoir "
                 "where the volumetric error (CNV) may be voilated "
                 "during strict Newton iterations.");
            Parameters::registerParam<Properties::ToleranceMb>
                ("Tolerated mass balance error relative to total mass present");
            Parameters::registerParam<Properties::ToleranceMbRelaxed>
                ("Relaxed tolerated mass balance error that applies for iterations "
                 "after the iterations with the strict tolerance");
            Parameters::registerParam<Properties::ToleranceCnv>
                ("Local convergence tolerance (Maximum of local saturation errors)");
            Parameters::registerParam<Properties::ToleranceCnvRelaxed>
                ("Relaxed local convergence tolerance that applies for iterations "
                 "after the iterations with the strict tolerance");
            Parameters::registerParam<Properties::ToleranceWells>
                ("Well convergence tolerance");
            Parameters::registerParam<Properties::ToleranceWellControl>
                ("Tolerance for the well control equations");
            Parameters::registerParam<Properties::MaxWelleqIter>
                ("Maximum number of iterations to determine solution the well equations");
            Parameters::registerParam<Properties::UseMultisegmentWell>
                ("Use the well model for multi-segment wells instead of the "
                 "one for single-segment wells");
            Parameters::registerParam<Properties::TolerancePressureMsWells>
                ("Tolerance for the pressure equations for multi-segment wells");
            Parameters::registerParam<Properties::RelaxedWellFlowTol>
                ("Relaxed tolerance for the well flow residual");
            Parameters::registerParam<Properties::RelaxedPressureTolMsw>
                ("Relaxed tolerance for the MSW pressure solution");
            Parameters::registerParam<Properties::MaxPressureChangeMsWells>
                ("Maximum relative pressure change for a single iteration "
                 "of the multi-segment well model");
            Parameters::registerParam<Properties::MaxInnerIterMsWells>
                ("Maximum number of inner iterations for multi-segment wells");
            Parameters::registerParam<Properties::StrictInnerIterWells>
                ("Number of inner well iterations with strict tolerance");
            Parameters::registerParam<Properties::StrictOuterIterWells>
                ("Number of newton iterations for which wells are checked with strict tolerance");
            Parameters::registerParam<Properties::MaxNewtonIterationsWithInnerWellIterations>
                ("Maximum newton iterations with inner well iterations");
            Parameters::registerParam<Properties::ShutUnsolvableWells>
                ("Shut unsolvable wells");
            Parameters::registerParam<Properties::MaxInnerIterWells>
                ("Maximum number of inner iterations for standard wells");
            Parameters::registerParam<Properties::AlternativeWellRateInit>
                ("Use alternative well rate initialization procedure");
            Parameters::registerParam<Properties::RegularizationFactorWells>
                ("Regularization factor for wells");
            Parameters::registerParam<Properties::MaxSinglePrecisionDays>
                ("Maximum time step size where single precision floating point "
                 "arithmetic can be used solving for the linear systems of equations");
            Parameters::registerParam<Properties::MinStrictCnvIter>
                ("Minimum number of Newton iterations before relaxed tolerances "
                 "can be used for the CNV convergence criterion");
            Parameters::registerParam<Properties::MinStrictMbIter>
                ("Minimum number of Newton iterations before relaxed tolerances "
                 "can be used for the MB convergence criterion. "
                 "Default -1 means that the relaxed tolerance is used when maximum "
                 "number of Newton iterations are reached.");
            Parameters::registerParam<Properties::SolveWelleqInitially>
                ("Fully solve the well equations before each iteration of the reservoir model");
            Parameters::registerParam<Properties::UpdateEquationsScaling>
                ("Update scaling factors for mass balance equations during the run");
            Parameters::registerParam<Properties::UseUpdateStabilization>
                ("Try to detect and correct oscillations or stagnation during the Newton method");
            Parameters::registerParam<Properties::MatrixAddWellContributions>
                ("Explicitly specify the influences of wells between cells in "
                 "the Jacobian and preconditioner matrices");
            Parameters::registerParam<Properties::EnableWellOperabilityCheck>
                ("Enable the well operability checking");
            Parameters::registerParam<Properties::EnableWellOperabilityCheckIter>
                ("Enable the well operability checking during iterations");
            Parameters::registerParam<Properties::MaximumNumberOfWellSwitches>
                ("Maximum number of times a well can switch to the same control");
            Parameters::registerParam<Properties::UseAverageDensityMsWells>
                ("Approximate segment densitities by averaging over segment and its outlet");
            Parameters::registerParam<Properties::LocalWellSolveControlSwitching>
                ("Allow control switching during local well solutions");
            Parameters::registerParam<Properties::UseImplicitIpr>
                ("Compute implict IPR for stability checks and stable solution search");
            Parameters::registerParam<Properties::NetworkMaxStrictIterations>
                ("Maximum iterations in network solver before relaxing tolerance");
            Parameters::registerParam<Properties::NetworkMaxIterations>
                ("Maximum number of iterations in the network solver before giving up");
            Parameters::registerParam<Properties::NonlinearSolver>
                ("Choose nonlinear solver. Valid choices are newton or nldd.");
            Parameters::registerParam<Properties::LocalSolveApproach>
                ("Choose local solve approach. Valid choices are jacobi and gauss-seidel");
            Parameters::registerParam<Properties::MaxLocalSolveIterations>
                ("Max iterations for local solves with NLDD nonlinear solver.");
            Parameters::registerParam<Properties::LocalToleranceScalingMb>
                ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
            Parameters::registerParam<Properties::LocalToleranceScalingCnv>
                ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
            Parameters::registerParam<Properties::NlddNumInitialNewtonIter>
                ("Number of initial global Newton iterations when running the NLDD nonlinear solver.");
            Parameters::registerParam<Properties::NumLocalDomains>
                ("Number of local domains for NLDD nonlinear solver.");
            Parameters::registerParam<Properties::LocalDomainsPartitioningImbalance>
                ("Subdomain partitioning imbalance tolerance. 1.03 is 3 percent imbalance.");
            Parameters::registerParam<Properties::LocalDomainsPartitioningMethod>
                ("Subdomain partitioning method. Allowed values are "
                 "'zoltan', "
                 "'simple', "
                 "and the name of a partition file ending with '.partition'.");
            Parameters::registerParam<Properties::LocalDomainsOrderingMeasure>
                ("Subdomain ordering measure. Allowed values are "
                 "'maxpressure', "
                 "'averagepressure' "
                 "and  'residual'.");
            Parameters::registerParam<Properties::DebugEmitCellPartition>
                ("Whether or not to emit cell partitions as a debugging aid.");

            Parameters::hideParam<Properties::DebugEmitCellPartition>();
        }
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELPARAMETERS_HEADER_INCLUDED
