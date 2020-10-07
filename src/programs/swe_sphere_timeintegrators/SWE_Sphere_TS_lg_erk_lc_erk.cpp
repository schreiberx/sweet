/*
 * SWE_Sphere_TS_lg_erk_lc_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_erk_lc_erk.hpp"



void SWE_Sphere_TS_lg_erk_lc_erk::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else
	{
		SWEETError("This time stepping order is not yet supported!");
	}
}





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (simVars.sim.sphere_use_fsphere)
	{
		double gh0 = simVars.sim.gravitation * simVars.sim.h0;

		o_phi_t = -gh0*i_div;
		o_div_t = -op.laplace(i_phi);

		o_vort_t = -simVars.sim.sphere_fsphere_f0*i_div;
		o_div_t += simVars.sim.sphere_fsphere_f0*i_vort;
	}
	else
	{
#if 0
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		o_phi_t = -gh*i_div;
		o_div_t = -op.laplace(i_phi);

		/*
		 * This doesn't converge to the reference solution
		 */
		SphereData_Spectral f(fg);
		o_vort_t = -f*i_div;
		o_div_t += f*i_vort;

#else

		double gh = simVars.sim.gravitation * simVars.sim.h0;

		/*
		 * Apply Coriolis Effect in physical VELOCITY space
		 */
		SphereData_Physical ug(i_phi.sphereDataConfig);
		SphereData_Physical vg(i_phi.sphereDataConfig);
		op.vrtdiv_to_uv(i_vort, i_div, ug, vg);

		SphereData_Physical tmpg1 = ug*op.fg;
		SphereData_Physical tmpg2 = vg*op.fg;

		op.uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

		o_vort_t *= -1.0;
		o_div_t += -op.laplace(i_phi);

		/*
		 * DIV on velocity field
		 */
		o_phi_t = (-gh)*i_div;
#endif
	}
}




void SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lg(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation*simVars.sim.h0;


	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 * Section "lg_erk"
	 */


	/*
	 * step 1a
	 */
	o_vort_t.spectral_set_zero();

	/*
	 * step 1b
	 */
	o_div_t = -op.laplace(i_phi);

	/*
	 * step 2a
	 */
	o_phi_t = -gh0*i_div;
}



void SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc(
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
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 * Section "lc_erk"
	 */

	/*
	 * step 1a
	 */
	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);
	op.vrtdiv_to_uv(i_vort, i_div, ug, vg);

	/*
	 * step 1b
	 */
	SphereData_Physical tmp_u = ug*op.fg;
	SphereData_Physical tmp_v = vg*op.fg;

	/*
	 * step 1c
	 */
	op.uv_to_vrtdiv(tmp_u, tmp_v, o_div_t, o_vort_t);

	/*
	 * step 1d
	 */
	o_vort_t *= -1.0;


	/*
	 * step 1e
	 * Nothing to do
	 */

	/*
	 * step 2a
	 * Zero tendencies
	 */
	o_phi_t.spectral_set_zero();

}



/**
 * This routine is used by other time step implementations
 */
void SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	SphereData_Spectral tmp_phi(io_phi.sphereDataConfig);
	SphereData_Spectral tmp_vort(io_vort.sphereDataConfig);
	SphereData_Spectral tmp_div(io_div.sphereDataConfig);

	euler_timestep_update_lc(
			io_phi,
			io_vort,
			io_div,

			tmp_phi,
			tmp_vort,
			tmp_div,

			i_simulation_timestamp
		);

	io_phi += i_dt*tmp_phi;
	io_vort += i_dt*tmp_vort;
	io_div += i_dt*tmp_div;
}



/*
 * Setup
 */
void SWE_Sphere_TS_lg_erk_lc_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


void SWE_Sphere_TS_lg_erk_lc_erk::setup_auto()
{
	setup(simVars.disc.timestepping_order);
}


SWE_Sphere_TS_lg_erk_lc_erk::SWE_Sphere_TS_lg_erk_lc_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
}



SWE_Sphere_TS_lg_erk_lc_erk::~SWE_Sphere_TS_lg_erk_lc_erk()
{
}

