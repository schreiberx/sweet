/*
 * SWE_Sphere_TS_lg_irk.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */


#include "SWE_Sphere_TS_lg_irk.hpp"
#include <complex>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>


SWE_Sphere_TS_lg_irk::SWE_Sphere_TS_lg_irk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
	simVars(i_simVars),
	op(i_op),
	sphereDataConfig(op.sphereDataConfig)
{
}


/**
 * Setup the SWE REXI solver with SPH
 */
void SWE_Sphere_TS_lg_irk::setup(
		int i_timestep_order,
		double i_timestep_size
)
{
	if (i_timestep_order != 1)
		FatalError("Only 1st order IRK supported so far!");

	timestep_size = i_timestep_size;

	alpha = -1.0/timestep_size;
	beta = -1.0/timestep_size;

	r = simVars.sim.earth_radius;
	inv_r = 1.0/r;

	gh = simVars.sim.gravitation*simVars.sim.h0;
}



/**
 * Solve a REXI time step for the given initial conditions
 */

void SWE_Sphere_TS_lg_irk::run_timestep(
		SphereData &io_phi,		///< prognostic variables
		SphereData &io_vort,	///< prognostic variables
		SphereData &io_div,		///< prognostic variables

		double &o_dt,				///< time step restriction
		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
	{
		FatalError("TODO: Reduction of time in SPH REXI required, not yet implemented");
		i_fixed_dt = i_max_simulation_time-i_simulation_timestamp;
	}

	o_dt = i_fixed_dt;

	SphereData phi0 = io_phi;
	SphereData vort0 = io_vort;
	SphereData div0 = io_div;

	SphereData phi(sphereDataConfig);
	SphereData vort(sphereDataConfig);
	SphereData div(sphereDataConfig);

	{
		SphereData rhs = gh*div0 + alpha*phi0;
		phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

		vort = 1.0/alpha*(vort0);
		div = -1.0/gh*(phi0 - alpha*phi);
	}

	io_phi = phi * beta;
	io_vort = vort * beta;
	io_div = div * beta;
}


SWE_Sphere_TS_lg_irk::~SWE_Sphere_TS_lg_irk()
{
}
