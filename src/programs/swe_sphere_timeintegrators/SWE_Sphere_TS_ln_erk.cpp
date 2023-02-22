/*
 * SWE_Sphere_TS_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "SWE_Sphere_TS_ln_erk.hpp"



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk::euler_timestep_update_pert(
		const SphereData_Spectral &i_phi_pert,	///< prognostic variables
		const SphereData_Spectral &i_vrt,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_pert_t,	///< time updates
		SphereData_Spectral &o_vrt_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation * simVars.sim.h0;

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
	SphereData_Physical phi_pert_phys = i_phi_pert.toPhys();

	/*
	 * Step 1a
	 */
	SphereData_Physical ug, vg;
	op.vrtdiv_to_uv(i_vrt, i_div, ug, vg);

	/*
	 * Step 1b
	 */
	SphereData_Physical vrtg = i_vrt.toPhys();

	/*
	 * Step 1c
	 */
	// left part of eq. (19)
	SphereData_Physical u_nl = ug*(vrtg+op.fg);

	// left part of eq. (20)
	SphereData_Physical v_nl = vg*(vrtg+op.fg);

	/*
	 * Step 1d
	 */
	// Eq. (21) & left part of Eq. (22)
	op.uv_to_vrtdiv(u_nl, v_nl, o_div_t, o_vrt_t);


	/*
	 * Step 1e
	 */
	o_vrt_t *= -1.0;

	/*
	 * Step 1f
	 */
	// Right part of Eq. (22)
	SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	SphereData_Spectral e = phi_pert_phys+tmpg;

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
	u_nl = ug*(phi_pert_phys + gh0);
	v_nl = vg*(phi_pert_phys + gh0);

	op.uv_to_vrtdiv(u_nl,v_nl, e, o_phi_pert_t);

	o_phi_pert_t *= -1.0;


}



void SWE_Sphere_TS_ln_erk::run_timestep(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_ln_erk::euler_timestep_update_pert,	///< pointer to function to compute euler time step updates
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
}


SWE_Sphere_TS_ln_erk::SWE_Sphere_TS_ln_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(timestepping_order);
}



SWE_Sphere_TS_ln_erk::~SWE_Sphere_TS_ln_erk()
{
}

