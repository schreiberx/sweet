/*
 * Adv_Sphere_TS_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "Adv_Sphere_TS_na_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void Adv_Sphere_TS_na_erk::euler_timestep_update(
		const SphereData &i_phi,	///< prognostic variables
		const SphereData &i_vort,	///< prognostic variables
		const SphereData &i_div,	///< prognostic variables

		SphereData &o_phi_t,	///< time updates
		SphereData &o_vort_t,	///< time updates
		SphereData &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{


#if 0
	SphereDataPhysical phi_grad_lon(i_phi.sphereDataConfig);
	SphereDataPhysical phi_grad_lat(i_phi.sphereDataConfig);
	phi_grad_lon = op.robert_grad_lon(i_phi).getSphereDataPhysical()/simVars.sim.earth_radius;
	phi_grad_lat = op.robert_grad_lat(i_phi).getSphereDataPhysical()/simVars.sim.earth_radius;
//	op.robert_grad_to_vec(i_phi, phi_grad_lon, phi_grad_lat, simVars.sim.earth_radius);

	SphereDataPhysical ug(i_phi.sphereDataConfig);
	SphereDataPhysical vg(i_phi.sphereDataConfig);
	op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);

	o_phi_t = ug*phi_grad_lon + vg*phi_grad_lat;
	o_vort_t.spectral_set_zero();
	o_div_t.spectral_set_zero();


#else

	/**
	 * We simply compute
	 * 	-DIV(rho*U) = -rho DIV(U) - U.GRAD(rho) = - U.GRAD(rho)
	 * which is the Lagrangian contribution only.
	 *
	 * This is the case because the velocity field is divergence free!!!
	 */

	SphereDataPhysical ug(i_phi.sphereDataConfig);
	SphereDataPhysical vg(i_phi.sphereDataConfig);

	SphereDataPhysical vrtg = i_vort.getSphereDataPhysical();
	SphereDataPhysical divg = i_div.getSphereDataPhysical();
	op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	SphereDataPhysical phig = i_phi.getSphereDataPhysical();

	SphereDataPhysical tmpg1 = ug*phig;
	SphereDataPhysical tmpg2 = vg*phig;

	SphereData tmpspec(i_phi.sphereDataConfig);
	op.robert_uv_to_vortdiv(tmpg1, tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	o_vort_t.spectral_set_zero();
	o_div_t.spectral_set_zero();

#endif
}



void Adv_Sphere_TS_na_erk::run_timestep(
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
			&Adv_Sphere_TS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void Adv_Sphere_TS_na_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Adv_Sphere_TS_na_erk::Adv_Sphere_TS_na_erk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Adv_Sphere_TS_na_erk::~Adv_Sphere_TS_na_erk()
{
}

