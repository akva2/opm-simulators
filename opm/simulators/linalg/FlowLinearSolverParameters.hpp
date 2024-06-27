/*
  Copyright 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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

#ifndef OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED
#define OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED

#include <opm/simulators/linalg/MILU.hpp>

#include <opm/simulators/linalg/linalgproperties.hh>
#include <opm/models/utils/parametersystem.hh>

namespace Opm {
template <class TypeTag>
class ISTLSolverBda;
template <class TypeTag>
class ISTLSolver;
}



namespace Opm::Properties {

namespace TTag {
struct FlowIstlSolverParams {};
}

// template<class TypeTag, class MyTypeTag>
// struct LinearSolverReduction {
//     using type = UndefinedProperty;
// };

// template<class TypeTag, class MyTypeTag>
// struct RelaxedLinearSolverReduction {
//     using type = UndefinedProperty;
// };

// template<class TypeTag, class MyTypeTag>
// struct LinearSolverMaxIter {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LinearSolverRestart {
//     using type = UndefinedProperty;
// };
//
// LinearSolverVerbosity defined in opm-models
//
// template<class TypeTag, class MyTypeTag>
// struct IluRelaxation {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct IluFillinLevel {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct MiluVariant {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct IluRedblack {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct IluReorderSpheres {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct UseGmres {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LinearSolverIgnoreConvergenceFailure{
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct ScaleLinearSystem {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LinearSolver {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct LinearSolverPrintJsonDefinition {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct CprReuseSetup {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct CprReuseInterval {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct AcceleratorMode {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct BdaDeviceId {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct OpenclPlatformId {
//     using type = UndefinedProperty;
// };
// template<class TypeTag, class MyTypeTag>
// struct OpenclIluParallel {
//     using type = UndefinedProperty;
// };
struct LinearSolverReduction {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};
struct RelaxedLinearSolverReduction {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};
struct LinearSolverMaxIter {
    static constexpr int value = 200;
};
struct LinearSolverRestart {
    static constexpr int value = 40;
};
// template<class TypeTag>
// struct LinearSolverVerbosity<TypeTag, TTag::FlowIstlSolverParams> {
//     static constexpr int value = 0;
// };
struct IluRelaxation {
    using type = double;//GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.9;
};
struct IluFillinLevel {
    static constexpr int value = 0;
};
struct MiluVariant {
    static constexpr auto value = "ILU";
};
struct IluRedblack {
    static constexpr bool value = false;
};
struct IluReorderSpheres {
    static constexpr bool value = false;
};
struct UseGmres {
    static constexpr bool value = false;
};
struct LinearSolverIgnoreConvergenceFailure {
    static constexpr bool value = false;
};
struct ScaleLinearSystem {
    static constexpr bool value = false;
};
struct LinearSolver {
    static constexpr auto value = "ilu0";
};
struct LinearSolverPrintJsonDefinition {
    static constexpr auto value = true;
};
struct CprReuseSetup {
    static constexpr int value = 4;
};
struct CprReuseInterval {
    static constexpr int value = 30;
};
struct AcceleratorMode {
    static constexpr auto value = "none";
};
struct BdaDeviceId {
    static constexpr int value = 0;
};
struct OpenclPlatformId {
    static constexpr int value = 0;
};
struct OpenclIluParallel {
    static constexpr bool value = true; // note: false should only be used in debug
};

// Set the backend to be used.
template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::FlowIstlSolverParams> {
#if COMPILE_BDA_BRIDGE
    using type = ISTLSolverBda<TypeTag>;
#else
    using type = ISTLSolver<TypeTag>;
#endif
};
} // namespace Opm::Properties

namespace Opm
{


    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
    struct FlowLinearSolverParameters
    {
        double linear_solver_reduction_;
        double relaxed_linear_solver_reduction_;
        int    linear_solver_maxiter_;
        int    linear_solver_restart_;
        int    linear_solver_verbosity_;
        double ilu_relaxation_;
        int    ilu_fillin_level_;
        MILU_VARIANT   ilu_milu_;
        bool   ilu_redblack_;
        bool   ilu_reorder_sphere_;
        bool   newton_use_gmres_;
        bool   ignoreConvergenceFailure_;
        bool scale_linear_system_;
        std::string linsolver_;
        bool linear_solver_print_json_definition_;
        int cpr_reuse_setup_;
        int cpr_reuse_interval_;
        std::string accelerator_mode_;
        int bda_device_id_;
        int opencl_platform_id_;
        bool opencl_ilu_parallel_;

        template <class TypeTag>
        void init(bool cprRequestedInDataFile)
        {
            // TODO: these parameters have undocumented non-trivial dependencies
            linear_solver_reduction_ = Parameters::get<Properties::LinearSolverReduction>();
            relaxed_linear_solver_reduction_ = Parameters::get<Properties::RelaxedLinearSolverReduction>();
            linear_solver_maxiter_ = Parameters::get<Properties::LinearSolverMaxIter>();
            linear_solver_restart_ = Parameters::get<Properties::LinearSolverRestart>();
            linear_solver_verbosity_ = Parameters::get<Properties::LinearSolverVerbosity>();
            ilu_relaxation_ = Parameters::get<Properties::IluRelaxation>();
            ilu_fillin_level_ = Parameters::get<Properties::IluFillinLevel>();
            ilu_milu_ = convertString2Milu(Parameters::get<Properties::MiluVariant>());
            ilu_redblack_ = Parameters::get<Properties::IluRedblack>();
            ilu_reorder_sphere_ = Parameters::get<Properties::IluReorderSpheres>();
            newton_use_gmres_ = Parameters::get<Properties::UseGmres>();
            ignoreConvergenceFailure_ = Parameters::get<Properties::LinearSolverIgnoreConvergenceFailure>();
            scale_linear_system_ = Parameters::get<Properties::ScaleLinearSystem>();
            linsolver_ = Parameters::get<Properties::LinearSolver>();
            linear_solver_print_json_definition_ = Parameters::get<Properties::LinearSolverPrintJsonDefinition>();
            cpr_reuse_setup_  = Parameters::get<Properties::CprReuseSetup>();
            cpr_reuse_interval_  = Parameters::get<Properties::CprReuseInterval>();

            if (!Parameters::isSet<Properties::LinearSolver>() && cprRequestedInDataFile) {
                linsolver_ = "cpr";
            } else {
                linsolver_ = Parameters::get<Properties::LinearSolver>();
            }

            accelerator_mode_ = Parameters::get<Properties::AcceleratorMode>();
            bda_device_id_ = Parameters::get<Properties::BdaDeviceId>();
            opencl_platform_id_ = Parameters::get<Properties::OpenclPlatformId>();
            opencl_ilu_parallel_ = Parameters::get<Properties::OpenclIluParallel>();
        }

        template <class TypeTag>
        static void registerParameters()
        {
            Parameters::registerParam<TypeTag, Properties::LinearSolverReduction>
                ("The minimum reduction of the residual which the linear solver must achieve");
            Parameters::registerParam<TypeTag, Properties::RelaxedLinearSolverReduction>
                ("The minimum reduction of the residual which the linear solver need to "
                 "achieve for the solution to be accepted");
            Parameters::registerParam<TypeTag, Properties::LinearSolverMaxIter>
                ("The maximum number of iterations of the linear solver");
            Parameters::registerParam<TypeTag, Properties::LinearSolverRestart>
                ("The number of iterations after which GMRES is restarted");
            Parameters::registerParam<TypeTag, Properties::LinearSolverVerbosity>
                ("The verbosity level of the linear solver (0: off, 2: all)");
            Parameters::registerParam<TypeTag, Properties::IluRelaxation>
                ("The relaxation factor of the linear solver's ILU preconditioner");
            Parameters::registerParam<TypeTag, Properties::IluFillinLevel>
                ("The fill-in level of the linear solver's ILU preconditioner");
            Parameters::registerParam<TypeTag, Properties::MiluVariant>
                ("Specify which variant of the modified-ILU preconditioner ought to be used. "
                 "Possible variants are: ILU (default, plain ILU), "
                 "MILU_1 (lump diagonal with dropped row entries), "
                 "MILU_2 (lump diagonal with the sum of the absolute values of the dropped row entries), "
                 "MILU_3 (if diagonal is positive add sum of dropped row entries, otherwise subtract them), "
                 "MILU_4 (if diagonal is positive add sum of dropped row entries, otherwise do nothing");
            Parameters::registerParam<TypeTag, Properties::IluRedblack>
                ("Use red-black partitioning for the ILU preconditioner");
            Parameters::registerParam<TypeTag, Properties::IluReorderSpheres>
                ("Whether to reorder the entries of the matrix in the red-black "
                 "ILU preconditioner in spheres starting at an edge. "
                 "If false the original ordering is preserved in each color. "
                 "Otherwise why try to ensure D4 ordering (in a 2D structured grid, "
                 "the diagonal elements are consecutive).");
            Parameters::registerParam<TypeTag, Properties::UseGmres>
                ("Use GMRES as the linear solver");
            Parameters::registerParam<TypeTag, Properties::LinearSolverIgnoreConvergenceFailure>
                ("Continue with the simulation like nothing happened "
                 "after the linear solver did not converge");
            Parameters::registerParam<TypeTag, Properties::ScaleLinearSystem>
                ("Scale linear system according to equation scale and primary variable types");
            Parameters::registerParam<TypeTag, Properties::LinearSolver>
                ("Configuration of solver. Valid options are: ilu0 (default), "
                 "dilu, cprw, cpr (an alias for cprw), cpr_quasiimpes, "
                 "cpr_trueimpes, cpr_trueimpesanalytic, amg or hybrid (experimental). "
                 "Alternatively, you can request a configuration to be read from a "
                 "JSON file by giving the filename here, ending with '.json.'");
            Parameters::registerParam<TypeTag, Properties::LinearSolverPrintJsonDefinition>
                ("Write the JSON definition of the linear solver setup to the DBG file.");
            Parameters::registerParam<TypeTag, Properties::CprReuseSetup>
                ("Reuse preconditioner setup. Valid options are "
                 "0: recreate the preconditioner for every linear solve, "
                 "1: recreate once every timestep, "
                 "2: recreate if last linear solve took more than 10 iterations, "
                 "3: never recreate, "
                 "4: recreated every CprReuseInterval");
            Parameters::registerParam<TypeTag, Properties::CprReuseInterval>
                ("Reuse preconditioner interval. Used when CprReuseSetup is set to 4, "
                 "then the preconditioner will be fully recreated instead of reused "
                 "every N linear solve, where N is this parameter.");
            Parameters::registerParam<TypeTag, Properties::AcceleratorMode>
                ("Choose a linear solver, usage: "
                 "'--accelerator-mode=[none|cusparse|opencl|amgcl|rocalution|rocsparse]'");
            Parameters::registerParam<TypeTag, Properties::BdaDeviceId>
                ("Choose device ID for cusparseSolver or openclSolver, "
                 "use 'nvidia-smi' or 'clinfo' to determine valid IDs");
            Parameters::registerParam<TypeTag, Properties::OpenclPlatformId>
                ("Choose platform ID for openclSolver, use 'clinfo' "
                 "to determine valid platform IDs");
            Parameters::registerParam<TypeTag, Properties::OpenclIluParallel>
                ("Parallelize ILU decomposition and application on GPU");
        }

        FlowLinearSolverParameters() { reset(); }

        // set default values
        void reset()
        {
            relaxed_linear_solver_reduction_ = 1e-2;
            linear_solver_reduction_  = 1e-2;
            linear_solver_maxiter_    = 200;
            linear_solver_restart_    = 40;
            linear_solver_verbosity_  = 0;
            ilu_relaxation_           = 0.9;
            ilu_fillin_level_         = 0;
            ilu_milu_                 = MILU_VARIANT::ILU;
            ilu_redblack_             = false;
            ilu_reorder_sphere_       = false;
            newton_use_gmres_         = false;
            ignoreConvergenceFailure_ = false;
            scale_linear_system_      = false;
            linsolver_                = "ilu0";
            linear_solver_print_json_definition_ = true;
            cpr_reuse_setup_          = 4;
            cpr_reuse_interval_       = 30;
            accelerator_mode_         = "none";
            bda_device_id_            = 0;
            opencl_platform_id_       = 0;
            opencl_ilu_parallel_      = true;
        }
    };


} // namespace Opm




#endif // OPM_FLOWLINEARSOLVERPARAMETERS_HEADER_INCLUDED
