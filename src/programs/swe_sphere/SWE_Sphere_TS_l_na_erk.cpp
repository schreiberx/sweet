/*
 * SWE_Sphere_TS_l_na_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_na_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_l_na_erk::euler_timestep_update(
		const SphereData &i_phi,	///< prognostic variables
		const SphereData &i_vort,	///< prognostic variables
		const SphereData &i_div,	///< prognostic variables

		SphereData &o_phi_t,	///< time updates
		SphereData &o_vort_t,	///< time updates
		SphereData &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	SphereDataPhysical ug(i_phi.sphereDataConfig);
	SphereDataPhysical vg(i_phi.sphereDataConfig);

	SphereDataPhysical vrtg = i_vort.getSphereDataPhysical();
	SphereDataPhysical divg = i_div.getSphereDataPhysical();
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	else
		op.vortdiv_to_uv(i_vort, i_div, ug, vg);
	SphereDataPhysical phig = i_phi.getSphereDataPhysical();

	SphereDataPhysical tmpg1 = ug*(vrtg+fg);
	SphereDataPhysical tmpg2 = vg*(vrtg+fg);

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);
	else
		op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	double gh = simVars.sim.gravitation * simVars.sim.h0;
	tmpg1 = ug*gh;
	tmpg2 = vg*gh;

	SphereData tmpspec(i_phi.sphereDataConfig);
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);
	else
		op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	SphereDataPhysical tmpg = 0.5*(ug*ug+vg*vg);

	if (simVars.misc.sphere_use_robert_functions)
		tmpg = tmpg.robert_convertToNonRobertSquared();

	tmpspec = phig+tmpg;
	o_div_t += -op.laplace(tmpspec);
}



void SWE_Sphere_TS_l_na_erk::run_timestep(
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
			&SWE_Sphere_TS_l_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_na_erk::setup(
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


SWE_Sphere_TS_l_na_erk::SWE_Sphere_TS_l_na_erk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_l_na_erk::~SWE_Sphere_TS_l_na_erk()
{
}

