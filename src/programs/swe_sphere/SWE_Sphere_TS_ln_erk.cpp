/*
 * SWE_Sphere_TS_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_ln_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk::euler_timestep_update(
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
	 * NON-LINEAR SWE
	 *
	 * See
	 * 	Williamson, David L., Drake, John B., Hack, James J., Jakob, Rudiger, & Swarztrauber, Paul N. (1992).
	 * 	A standard test set for numerical approximations to the shallow water equations in spherical geometry.
	 * 	Journal of Computational Physics, 102(1), 211–224. https://doi.org/10.1016/S0021-9991(05)80016-6
	 *
	 * "2.3 Vorticity/Divergence Form"
	 */

	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 */
	SphereData_Physical phig = i_phi.getSphereDataPhysical();


	/*
	 * Step 1a
	 */
	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	else
		op.vortdiv_to_uv(i_vort, i_div, ug, vg);

	/*
	 * Step 1b
	 */
	SphereData_Physical vrtg = i_vort.getSphereDataPhysical();

	/*
	 * Step 1c
	 */
	// left part of eq. (19)
	SphereData_Physical tmp_u = ug*(vrtg+fg);

	// left part of eq. (20)
	SphereData_Physical tmp_v = vg*(vrtg+fg);

	/*
	 * Step 1d
	 */
	// Eq. (21) & left part of Eq. (22)
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmp_u, tmp_v, o_div_t, o_vort_t);
	else
		op.uv_to_vortdiv(tmp_u, tmp_v, o_div_t, o_vort_t);
	
	/*
	 * Step 1e
	 */
	o_vort_t *= -1.0;

	/*
	 * Step 1f
	 */
	// Right part of Eq. (22)
	SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	if (simVars.misc.sphere_use_robert_functions)
		tmpg = tmpg.robert_convertToNonRobertSquared();

	SphereData_Spectral e = phig+tmpg;

	/*
	 * Step 1g
	 */
	o_div_t -= op.laplace(e);

	/*
	 * Compute Phi geopotential tendencies
	 */

	/*
	 * Step 2a
	 */
	tmp_u = ug*phig;
	tmp_v = vg*phig;

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmp_u,tmp_v, e, o_phi_t);
	else
		op.uv_to_vortdiv(tmp_u,tmp_v, e, o_phi_t);
	o_phi_t *= -1.0;


}



void SWE_Sphere_TS_ln_erk::run_timestep(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_ln_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_ln_erk::setup(
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


SWE_Sphere_TS_ln_erk::SWE_Sphere_TS_ln_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_ln_erk::~SWE_Sphere_TS_ln_erk()
{
}

