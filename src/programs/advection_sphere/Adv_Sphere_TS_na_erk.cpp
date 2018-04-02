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

#if 0
	std::cout << phig.physical_reduce_max_abs() << std::endl;
	std::cout << ug.physical_reduce_max_abs() << std::endl;
	std::cout << o_phi_t.physical_reduce_max_abs() << std::endl;
	std::cout << std::endl;
#endif

	o_vort_t.spectral_set_zero();
	o_div_t.spectral_set_zero();
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

