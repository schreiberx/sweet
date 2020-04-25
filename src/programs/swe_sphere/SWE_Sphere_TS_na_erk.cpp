/*
 * SWE_Sphere_TS_na_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_na_erk.hpp"




/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_na_erk::euler_timestep_update(
		const SphereData_Spectral &i_phi_pert,	///< prognostic variables
		const SphereData_Spectral &i_vrt,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_pert_t,	///< time updates
		SphereData_Spectral &o_vrt_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	const SphereData_Config* sphereDataConfig = i_phi_pert.sphereDataConfig;

	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 */

	/*
	 * Step 1a
	 */
	SphereData_Physical ug(sphereDataConfig);
	SphereData_Physical vg(sphereDataConfig);

	op.vortdiv_to_uv(i_vrt, i_div, ug, vg, simVars.misc.sphere_use_robert_functions);

	/*
	 * Step 1b
	 */
	SphereData_Physical vrtg = i_vrt.getSphereDataPhysical();

	/*
	 * Step 1c
	 */
	// left part of eq. (19)

	// !!! na_erk change !!!
	//SphereData_Physical u_nl = ug*(vrtg+op.fg);
	SphereData_Physical u_nl = ug*vrtg;

	// left part of eq. (20)

	// !!! na_erk change !!!
	//SphereData_Physical v_nl = vg*(vrtg+op.fg);
	SphereData_Physical v_nl = vg*vrtg;

	/*
	 * Step 1d
	 */
	// Eq. (21) & left part of Eq. (22)
	op.uv_to_vortdiv(u_nl, v_nl, o_div_t, o_vrt_t, simVars.misc.sphere_use_robert_functions);


	/*
	 * Step 1e
	 */
	o_vrt_t *= -1.0;


	/*
	 * Step 1f
	 */
	// Right part of Eq. (22)
	SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	if (simVars.misc.sphere_use_robert_functions)
		tmpg = tmpg.robert_convertToNonRobertSquared();

	/*
	 * Step 1h
	 */
	o_div_t -= op.laplace(tmpg);


	/*
	 * Compute Phi geopotential tendencies
	 */

	/*
	 * Step 2a
	 */
	SphereData_Physical phig_pert = i_phi_pert.getSphereDataPhysical();

	/*
	 * Step 2b
	 */
	u_nl = ug*phig_pert;
	v_nl = vg*phig_pert;

	/*
	 * Step 2c
	 */
	SphereData_Spectral e(sphereDataConfig);
	op.uv_to_vortdiv(u_nl, v_nl, e, o_phi_pert_t, simVars.misc.sphere_use_robert_functions);

	/*
	 * Step 2d
	 */

	o_phi_pert_t *= -1.0;

	/*
	 * Step 2e
	 */
	SphereData_Physical divg = i_div.getSphereDataPhysical();

	/*
	 * Step 2f
	 */
	e = op.scalar_physical_to_spectral(divg*phig_pert);

	/*
	 * Step 2g
	 */
	o_phi_pert_t += e;
}



void SWE_Sphere_TS_na_erk::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,		///< prognostic variables
		SphereData_Spectral &io_vrt,		///< prognostic variables
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
			&SWE_Sphere_TS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_na_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


SWE_Sphere_TS_na_erk::SWE_Sphere_TS_na_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_na_erk::~SWE_Sphere_TS_na_erk()
{
}

