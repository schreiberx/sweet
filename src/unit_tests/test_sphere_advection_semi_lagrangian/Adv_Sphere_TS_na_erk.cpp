/*
 * Adv_Sphere_TS_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "Adv_Sphere_TS_na_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void Adv_Sphere_TS_na_erk::euler_timestep_update(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

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

	if (simVars.benchmark.getExternalForcesCallback != nullptr)
	{
		SphereData_Spectral vort(i_phi.sphereDataConfig);
		SphereData_Spectral div(i_phi.sphereDataConfig);

		simVars.benchmark.getExternalForcesCallback(1, simVars.timecontrol.current_simulation_time, &vort, simVars.benchmark.getExternalForcesUserData);
		simVars.benchmark.getExternalForcesCallback(2, simVars.timecontrol.current_simulation_time, &div, simVars.benchmark.getExternalForcesUserData);

		SphereData_Physical ug(i_phi.sphereDataConfig);
		SphereData_Physical vg(i_phi.sphereDataConfig);

		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
		SphereData_Physical phig = i_phi.getSphereDataPhysical();

		SphereData_Physical tmpg1 = ug*phig;
		SphereData_Physical tmpg2 = vg*phig;

		SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, tmpspec, o_phi_t);

		o_phi_t *= -1.0;
	}
	else
	{
		SphereData_Physical ug(i_phi.sphereDataConfig);
		SphereData_Physical vg(i_phi.sphereDataConfig);

		SphereData_Physical vrtg = i_vort.getSphereDataPhysical();
		SphereData_Physical divg = i_div.getSphereDataPhysical();
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
		SphereData_Physical phig = i_phi.getSphereDataPhysical();

		SphereData_Physical tmpg1 = ug*phig;
		SphereData_Physical tmpg2 = vg*phig;

		SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, tmpspec, o_phi_t);

		o_phi_t *= -1.0;
	}

	o_vort_t.spectral_set_zero();
	o_div_t.spectral_set_zero();
}



void Adv_Sphere_TS_na_erk::run_timestep(
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
			&Adv_Sphere_TS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);

	if (simVars.benchmark.getExternalForcesCallback != nullptr)
	{
		// this is just called for cosmetic reasons to update the velocity field
		simVars.benchmark.getExternalForcesCallback(1, simVars.timecontrol.current_simulation_time+i_fixed_dt, &io_vort, simVars.benchmark.getExternalForcesUserData);
		simVars.benchmark.getExternalForcesCallback(2, simVars.timecontrol.current_simulation_time+i_fixed_dt, &io_div, simVars.benchmark.getExternalForcesUserData);
	}
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
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Adv_Sphere_TS_na_erk::~Adv_Sphere_TS_na_erk()
{
}

