/*
 * SWE_Sphere_TS_l_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_l_erk::euler_timestep_update(
		const SphereDataSpectral &i_phi,	///< prognostic variables
		const SphereDataSpectral &i_vort,	///< prognostic variables
		const SphereDataSpectral &i_div,	///< prognostic variables

		SphereDataSpectral &o_phi_t,	///< time updates
		SphereDataSpectral &o_vort_t,	///< time updates
		SphereDataSpectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (simVars.sim.sphere_use_fsphere)
	{
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		o_phi_t = -gh*i_div;
		o_div_t = -op.laplace(i_phi);

		o_vort_t = -simVars.sim.sphere_fsphere_f0*i_div;
		o_div_t += simVars.sim.sphere_fsphere_f0*i_vort;
	}
	else
	{
#if 0
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		o_phi_t = -gh*i_div;
		o_div_t = -op.laplace(i_phi);

		/*
		 * This doesn't converge to the reference solution
		 */
		SphereDataSpectral f(fg);
		o_vort_t = -f*i_div;
		o_div_t += f*i_vort;

#else

		double gh = simVars.sim.gravitation * simVars.sim.h0;

		/*
		 * Apply Coriolis Effect in physical VELOCITY space
		 */
		SphereDataPhysical ug(i_phi.sphereDataConfig);
		SphereDataPhysical vg(i_phi.sphereDataConfig);
		if (simVars.misc.sphere_use_robert_functions)
			op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
		else
			op.vortdiv_to_uv(i_vort, i_div, ug, vg);

		SphereDataPhysical tmpg1 = ug*fg;
		SphereDataPhysical tmpg2 = vg*fg;

		if (simVars.misc.sphere_use_robert_functions)
			op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);
		else
			op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

		o_vort_t *= -1.0;
		o_div_t += -op.laplace(i_phi);

		/*
		 * DIV on velocity field
		 */
		o_phi_t = (-gh)*i_div;
#endif
	}
}



void SWE_Sphere_TS_l_erk::run_timestep(
		SphereDataSpectral &io_phi,		///< prognostic variables
		SphereDataSpectral &io_vort,	///< prognostic variables
		SphereDataSpectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_l_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_erk::setup(
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


SWE_Sphere_TS_l_erk::SWE_Sphere_TS_l_erk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_l_erk::~SWE_Sphere_TS_l_erk()
{
}

