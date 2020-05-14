/*
 * SWE_Sphere_TS_PFASST_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_PFASST_ln_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_PFASST_ln_erk::euler_timestep_update_nonpert(
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
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	SphereData_Physical vrtg = i_vort.getSphereDataPhysical();
	SphereData_Physical divg = i_div.getSphereDataPhysical();
	op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	SphereData_Physical phig = i_phi.getSphereDataPhysical();

	SphereData_Physical tmpg1 = ug*(vrtg+fg);
	SphereData_Physical tmpg2 = vg*(vrtg+fg);

	op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	SphereData_Physical tmpg = o_div_t.getSphereDataPhysical();

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	tmpspec = (phig+0.5*(ug*ug+vg*vg));
//	tmpspec.request_data_spectral();
	o_div_t -= op.laplace(tmpspec);
}



void SWE_Sphere_TS_PFASST_ln_erk::run_timestep_nonpert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_PFASST_ln_erk::euler_timestep_update_nonpert,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_PFASST_ln_erk::setup(
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


SWE_Sphere_TS_PFASST_ln_erk::SWE_Sphere_TS_PFASST_ln_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_PFASST_ln_erk::~SWE_Sphere_TS_PFASST_ln_erk()
{
}

