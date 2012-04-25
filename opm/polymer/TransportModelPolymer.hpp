/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_TRANSPORTMODELPOLYMER_HEADER_INCLUDED
#define OPM_TRANSPORTMODELPOLYMER_HEADER_INCLUDED

#include <opm/polymer/PolymerProperties.hpp>
#include <opm/core/transport/reorder/TransportModelInterface.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <vector>

class UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;

    /// A transport model for two-phase flow with polymer in the
    /// water phase.
    /// \TODO Include permeability reduction effect.
    class TransportModelPolymer : public TransportModelInterface
    {
    public:

	enum SingleCellMethod { Bracketing, Newton };
        enum GradientMethod { Analytic, FinDif }; // Analytic is chosen (hard-coded)

	/// \TODO document me, especially method.
	TransportModelPolymer(const UnstructuredGrid& grid,
			      const double* porosity,
			      const double* porevolume,
			      const IncompPropertiesInterface& props,
			      const PolymerProperties& polyprops,
			      const SingleCellMethod method,
			      const double tol,
			      const int maxit);

	/// Solve transport eqn with implicit Euler scheme, reordered.
	/// \TODO Now saturation is expected to be one sw value per cell,
	/// change to [sw so] per cell.
	void solve(const double* darcyflux,
		   const double* source,
		   const double dt,
		   const double inflow_c,
		   double* saturation,
		   double* concentration,
		   double* cmax);

	virtual void solveSingleCell(const int cell);
	virtual void solveMultiCell(const int num_cells, const int* cells);
	void solveSingleCellBracketing(int cell);
	void solveSingleCellNewton(int cell);
	class ResidualEquation;


    private:
	const UnstructuredGrid& grid_;
	const double* porosity_;
	const double* porevolume_;  // one volume per cell
	const IncompPropertiesInterface& props_;
	const PolymerProperties& polyprops_;
	std::vector<double> smin_;
	std::vector<double> smax_;
	double tol_;
	double maxit_;

	const double* darcyflux_;   // one flux per grid face
	const double* source_;      // one source per cell
	double dt_;
	double inflow_c_;
	double* saturation_;        // one per cell
	double* concentration_;
	double* cmax_;
	std::vector<double> fractionalflow_;  // one per cell
	std::vector<double> mc_;  // one per cell
	const double* visc_;
	SingleCellMethod method_;

	struct ResidualC;
	struct ResidualS;


	void fracFlow(double s, double c, double cmax, int cell, double& ff) const;
	void fracFlowWithDer(double s, double c, double cmax, int cell, double& ff,
                               double* dff_dsdc) const;
	void fracFlowBoth(double s, double c, double cmax, int cell, double& ff,
                          double* dff_dsdc, bool if_with_der) const;
	void computeMc(double c, double& mc) const;
	void computeMcWithDer(double c, double& mc, double& dmc_dc) const;
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELPOLYMER_HEADER_INCLUDED
