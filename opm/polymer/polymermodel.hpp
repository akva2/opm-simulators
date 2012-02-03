/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */
/* Copyright 2012 (c) SINTEF */

#ifndef POLYMER_HPP_INCLUDED
#define POLYMER_HPP_INCLUDED

#include <opm/core/utility/linearInterpolation.hpp>

struct UnstructuredGrid;
namespace Opm
{
    class IncompPropertiesInterface;
}

struct PolymerData
{
    double c_max_limit;
    double omega;
    double viscMult(double c) const
    {
	return Opm::linearInterpolation(c_vals, visc_mult_vals, c);
    }

    std::vector<double> c_vals;
    std::vector<double> visc_mult_vals;
};


struct PolymerSolverData  {
    struct UnstructuredGrid *grid;
    const Opm::IncompPropertiesInterface* props;
    const PolymerData* polydata;
    const double            *darcyflux;   /* one flux per face  in cdata::grid*/
    const double            *porevolume;  /* one volume per cell */
    const double            *source;      /* one source per cell */
    double                   dt;
    double                  *saturation;      /* one per cell */
    double                  *concentration;   /* one per cell */
    double                  *fractionalflow;  /* one per cell */
};

struct NonlinearSolverCtrl;


void
polymer_solvecell(void *data, struct NonlinearSolverCtrl *ctrl, int cell);

void
destroy_solverdata(struct PolymerSolverData *d);

struct PolymerSolverData *
init_solverdata(struct UnstructuredGrid *grid,
		const Opm::IncompPropertiesInterface* props,
		const PolymerData* polydata,
		const double *darcyflux,
                const double *porevolume,
		const double *source,
                const double dt,
		double *saturation,
		double *concentration);


#endif /* POLYMER_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
