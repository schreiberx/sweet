/*
 * SWE_Sphere_TS_l_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_erk.hpp"



void SWE_Sphere_TS_l_erk::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation*simVars.sim.h0;
	io_phi_pert += gh0;
	run_timestep_nonpert(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
	io_phi_pert -= gh0;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_l_erk::euler_timestep_update(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (!simVars.sim.sphere_use_fsphere)
	{

		double gh0 = simVars.sim.gravitation * simVars.sim.h0;

		/*
		 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
		 */
		/*
		 * Step 1a
		 */
		SphereData_Physical ug(i_phi.sphereDataConfig);
		SphereData_Physical vg(i_phi.sphereDataConfig);
		op.vortdiv_to_uv(i_vort, i_div, ug, vg, simVars.misc.sphere_use_robert_functions);

		/*
		 * Step 1b
		 */
		SphereData_Physical tmpg1 = ug*op.fg;
		SphereData_Physical tmpg2 = vg*op.fg;

		/*
		 * Step 1c
		 */
		op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t, simVars.misc.sphere_use_robert_functions);

		/*
		 * Step 1d
		 */
		o_vort_t *= -1.0;

		/*
		 * Step 1e
		 */
		o_div_t += -op.laplace(i_phi);

		/*
		 * DIV on velocity field
		 */
		o_phi_t = (-gh0)*i_div;
	}
	else
	{
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		o_div_t = -op.laplace(i_phi);

		o_vort_t = -simVars.sim.sphere_fsphere_f0*i_div;
		o_div_t += simVars.sim.sphere_fsphere_f0*i_vort;

		o_phi_t = -gh*i_div;
	}
}





void SWE_Sphere_TS_l_erk::run_timestep_nonpert(
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
}


SWE_Sphere_TS_l_erk::SWE_Sphere_TS_l_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_l_erk::~SWE_Sphere_TS_l_erk()
{
}

