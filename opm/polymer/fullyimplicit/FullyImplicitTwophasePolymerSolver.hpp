/**/

#ifndef OPM_FULLYIMPLICITTWOPHASEPOLYMERSOLVER_HEADER_INCLUDED
#define OPM_FULLYIMPLICITTWOPHASEPOLYMERSOLVER_HEADER_INCLUDED

#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
#include <opm/polymer//fullyimplicit/AutoDiffHelpers.hpp>
#include <opm/polymer/fullyimplicit/IncompPropsAdInterface.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>


struct UnstructuredGrid;
namespace Opm {
    class LinearSolverInterface;
    class PolymerState;

    
    class FullyImplicitTwophasePolymerSolver
    {
    public:
        FullyImplicitTwophasePolymerSolver(const UnstructuredGrid&        grid,
                                           const IncompPropsAdInterface&  fluid,
                                           const PolymerPropsAd&          polymer_props_ad,
                                           const LinearSolverInterface&    linsolver);

        void step(const double   dt,
                  PolymerState& state,
                  const std::vector<double>& src,
                  const std::vector<double>& polymer_inflow
                  );
    private:
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;
        struct SolutionState {
            SolutionState(const int np);
            ADB             pressure;
            std::vector<ADB> saturation;
            ADB             concentration;
        };
        const UnstructuredGrid&         grid_;
        const IncompPropsAdInterface&   fluid_;
        const PolymerPropsAd&           polymer_props_ad_;
        const LinearSolverInterface&    linsolver_;
        const std::vector<int>          cells_;
        HelperOps                       ops_;
        std::vector<ADB>                residual_;
       

        SolutionState
        constantState(const PolymerState& x);
        SolutionState
        variableState(const PolymerState& x);
        void
        assemble(const V&               pvdt,
                 const SolutionState&   old_state,
                 const PolymerState&  x,
                 const std::vector<double>& src,
                 const std::vector<double>& polymer_inflow);
        V solveJacobianSystem() const;
        void updateState(const V&             dx,
                         PolymerState& x) const;
        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;
        V
        transmissibility() const;
        
        std::vector<ADB>
        computeMassFlux(const V&                trans,
                        const ADB&              mc,
                        const ADB&              kro,
                        const ADB&              krw_eff,
                        const SolutionState&    state ) const;
    
        std::vector<ADB>
        accumSource(const ADB&                 kro,
                    const ADB&                 krw_eff,
                    const ADB&                 c,
                    const std::vector<double>& src,
                    const std::vector<double>& polymer_inflow_c) const;

        
        std::vector<ADB>
        computeFracFlow(const ADB& kro,
                        const ADB& krw_eff,
                        const ADB& c) const;
        double
        residualNorm() const;
        ADB
        polymerSource(const std::vector<ADB>& kr,
                      const std::vector<double>& src,
                      const std::vector<double>& polymer_inflow_c,
                      const SolutionState& state) const;

        ADB
        computeCmax(const ADB& c) const;
        ADB 
        computeMc(const SolutionState& state) const;
        ADB
        rockPorosity(const ADB& p) const;
        ADB
        rockPermeability(const ADB& p) const;
        const double
        fluidDensity(const int phase) const;
        ADB
        transMult(const ADB& p) const;
    };


} // namespace Opm
#endif// OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
