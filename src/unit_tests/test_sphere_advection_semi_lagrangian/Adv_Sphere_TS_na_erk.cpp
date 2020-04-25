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
		const SphereData_Spectral &i_vrt,	///< prognostic variables
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

	SphereData_Spectral phi = i_phi;
	SphereData_Spectral vrt = i_vrt;
	SphereData_Spectral div = i_div;

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (sphereBenchmarks)
		sphereBenchmarks->update_time_varying_fields_pert(phi, vrt, div, i_simulation_timestamp);

	SphereData_Physical ug(phi.sphereDataConfig);
	SphereData_Physical vg(phi.sphereDataConfig);
	op.robert_vortdiv_to_uv(vrt, div, ug, vg);

	SphereData_Physical vrtg = vrt.getSphereDataPhysical();
	SphereData_Physical divg = div.getSphereDataPhysical();

	SphereData_Physical phig = phi.getSphereDataPhysical();

	SphereData_Physical tmpg1 = ug*phig;
	SphereData_Physical tmpg2 = vg*phig;

	SphereData_Spectral tmpspec(phi.sphereDataConfig);
	op.uv_to_vortdiv(tmpg1, tmpg2, tmpspec, o_phi_t, simVars.misc.sphere_use_robert_functions);

	o_phi_t *= -1.0;

	o_vort_t.spectral_set_zero();
	o_div_t.spectral_set_zero();
}



void Adv_Sphere_TS_na_erk::run_timestep(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const SWESphereBenchmarksCombined *i_sphereBenchmarks,
		SphereData_Physical &io_U_phi_phys
)
{
	sphereBenchmarks = i_sphereBenchmarks;

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

	io_U_phi_phys = io_phi.getSphereDataPhysical();
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

