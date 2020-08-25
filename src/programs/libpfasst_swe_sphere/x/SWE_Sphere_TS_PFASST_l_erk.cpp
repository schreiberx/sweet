/*
 * SWE_Sphere_TS_PFASST_l_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_PFASST_l_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_PFASST_l_erk::euler_timestep_update(
		const SphereData_Spectral &i_phi_pert,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_pert_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (simVars.sim.sphere_use_fsphere)
	{
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		o_phi_pert_t = -gh*i_div;
		o_div_t = -op.laplace(i_phi_pert);

		o_vort_t = -simVars.sim.plane_rotating_f0*i_div;
		o_div_t += simVars.sim.plane_rotating_f0*i_vort;
	}
	else
	{
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		SphereData_Physical ug(i_phi_pert.sphereDataConfig);
		SphereData_Physical vg(i_phi_pert.sphereDataConfig);
		op.robert_vrtdiv_to_uv(i_vort, i_div, ug, vg);

		SphereData_Physical tmpg1 = ug*fg;
		SphereData_Physical tmpg2 = vg*fg;

		op.robert_uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

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

		SphereData_Spectral tmpspec(i_phi_pert.sphereDataConfig);
		op.robert_uv_to_vrtdiv(tmpg1,tmpg2, tmpspec, o_phi_pert_t);

		o_phi_pert_t *= -1.0;

		SphereData_Physical phig = i_phi_pert.toPhys();
		tmpspec.loadSphereDataPhysical(phig);
		o_div_t += -op.laplace(tmpspec);
	}
}



void SWE_Sphere_TS_PFASST_l_erk::run_timestep_pert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_PFASST_l_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_PFASST_l_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	if (simVars.sim.sphere_use_fsphere)
	{
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = simVars.sim.plane_rotating_f0;
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


SWE_Sphere_TS_PFASST_l_erk::SWE_Sphere_TS_PFASST_l_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_PFASST_l_erk::~SWE_Sphere_TS_PFASST_l_erk()
{
}

