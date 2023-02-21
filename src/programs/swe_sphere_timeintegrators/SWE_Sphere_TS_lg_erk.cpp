/*
 * SWE_Sphere_TS_lg_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "SWE_Sphere_TS_lg_erk.hpp"



void SWE_Sphere_TS_lg_erk::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_lg_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_lg_erk::euler_timestep_update(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

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
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(timestepping_order);
}



SWE_Sphere_TS_lg_erk::~SWE_Sphere_TS_lg_erk()
{
}

