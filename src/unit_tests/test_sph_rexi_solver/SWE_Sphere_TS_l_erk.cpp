#error "This file might be obsolete!"
#error "This file might be obsolete!"
#error "This file might be obsolete!"


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
		const SphereData &i_phi,	///< prognostic variables
		const SphereData &i_vort,	///< prognostic variables
		const SphereData &i_div,	///< prognostic variables

		SphereData &o_phi_t,	///< time updates
		SphereData &o_vort_t,	///< time updates
		SphereData &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (simVars.sim.sphere_use_fsphere)
	{
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		o_phi_t = -gh*i_div;
		o_div_t = -op.laplace(i_phi);

		o_vort_t = -simVars.sim.plane_rotating_f0*i_div;
		o_div_t += simVars.sim.plane_rotating_f0*i_vort;
	}
	else
	{
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		SphereData_Physical ug(i_phi.sphereDataConfig);
		SphereData_Physical vg(i_phi.sphereDataConfig);
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
		SphereData_Physical phig = i_phi.getSphereDataPhysical();

		SphereData_Physical tmpg1 = ug*op.fg;
		SphereData_Physical tmpg2 = vg*op.fg;

		op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

		o_vort_t *= -1.0;

#if 0
		// Non-linear divergence
		tmpg1 = ug*phig;
		tmpg2 = vg*phig;
#else
		// Linearized divergence around avg. geopotential
		tmpg1 = ug*gh;
		tmpg2 = vg*gh;
#endif

		SphereData tmpspec(i_phi.sphereDataConfig);
		op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

		o_phi_t *= -1.0;

		tmpspec = phig;
		tmpspec.request_data_spectral();
		o_div_t += -op.laplace(tmpspec);
	}
}



void SWE_Sphere_TS_l_erk::run_timestep(
		SphereData &io_phi,		///< prognostic variables
		SphereData &io_vort,	///< prognostic variables
		SphereData &io_div,		///< prognostic variables

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

