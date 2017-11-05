/*
 * SWE_Sphere_TS_lg_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Sphere_TS_lg_erk.hpp"



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_lg_erk::euler_timestep_update(
		const SphereData &i_phi,	///< prognostic variables
		const SphereData &i_vort,	///< prognostic variables
		const SphereData &i_div,	///< prognostic variables

		SphereData &o_phi_t,	///< time updates
		SphereData &o_vort_t,	///< time updates
		SphereData &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	/*
	 * LINEAR
	 */
	double gh = simVars.sim.gravitation * simVars.sim.h0;

	o_phi_t = -gh*i_div;
	o_div_t = -op.laplace(i_phi);
	o_vort_t.spectral_set_zero();
}



void SWE_Sphere_TS_lg_erk::run_timestep(
		SphereData &io_phi,		///< prognostic variables
		SphereData &io_vort,	///< prognostic variables
		SphereData &io_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_lg_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_lg_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


SWE_Sphere_TS_lg_erk::SWE_Sphere_TS_lg_erk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_lg_erk::~SWE_Sphere_TS_lg_erk()
{
}

