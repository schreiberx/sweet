/*
 * SWE_Sphere_TS_l_erk_pvd.cpp
 *
 *  Created on: 31 Oct 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 * Implementation which solely uses phi, vort and div fields
 */

#include "SWE_Sphere_TS_l_erk_pvd.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_l_erk_pvd::euler_timestep_update(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (simVars.sim.sphere_use_fsphere)
		FatalError("f-sphere Not supported");


	double gh = simVars.sim.gravitation * simVars.sim.h0;

	o_phi_t = -gh*i_div;
	o_div_t = -op.laplace(i_phi);

	/*
	 * WARNING:
	 * This formulation has problems with the highest frequency which results
	 * likely in an aliased one
	 */
	SphereData_Spectral f(fg);
	o_vort_t = -f*i_div;
	o_div_t += f*i_vort;
}



void SWE_Sphere_TS_l_erk_pvd::run_timestep(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_l_erk_pvd::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_erk_pvd::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	if (simVars.sim.sphere_use_fsphere)
	{
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = simVars.sim.sphere_fsphere_f0;
			}
		);
	}
	else
	{
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu*2.0*simVars.sim.sphere_rotating_coriolis_omega;
			}
		);
	}

}


SWE_Sphere_TS_l_erk_pvd::SWE_Sphere_TS_l_erk_pvd(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_l_erk_pvd::~SWE_Sphere_TS_l_erk_pvd()
{
}

