/*
 * SWE_Sphere_TS_lg_erk_lf_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_erk_lc_erk.hpp"




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
		double gh = simVars.sim.gravitation * simVars.sim.h0;

		o_phi_t = -gh*i_div;
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
		if (simVars.misc.sphere_use_robert_functions)
			op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
		else
			op.vortdiv_to_uv(i_vort, i_div, ug, vg);

		SphereData_Physical tmpg1 = ug*fg;
		SphereData_Physical tmpg2 = vg*fg;

		if (simVars.misc.sphere_use_robert_functions)
			op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);
		else
			op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

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
	double gh = simVars.sim.gravitation*simVars.sim.h0;

#if 1

	o_phi_t = -gh*i_div;
	o_div_t = -op.laplace(i_phi);
	o_vort_t.spectral_set_zero();

#elif 0

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	else
		op.vortdiv_to_uv(i_vort, i_div, ug, vg);

	// Nonlinearity
	SphereData_Physical tmpg1 = ug*i_vort.getSphereDataPhysical();
	SphereData_Physical tmpg2 = vg*i_div.getSphereDataPhysical();

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, tmpspec, o_phi_t);
	else
		op.uv_to_vortdiv(tmpg1, tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

//	o_phi_t = -gh*i_div;
	o_div_t = -op.laplace(i_phi);
	o_vort_t.spectral_set_zero();


#else

	/*
	 * Apply Coriolis Effect in physical VELOCITY space
	 */
	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	else
		op.vortdiv_to_uv(i_vort, i_div, ug, vg);

	SphereData_Physical tmpg1 = ug*fg;
	SphereData_Physical tmpg2 = vg*fg;

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);
	else
		op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	o_vort_t.spectral_set_zero();
	o_div_t.spectral_set_zero();

	o_div_t += -op.laplace(i_phi);

	/*
	 * DIV on velocity field
	 */
	o_phi_t = (-gh)*i_div;
#endif
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
#if 1
	/*
	 * Apply Coriolis Effect in physical VELOCITY space
	 */
	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	else
		op.vortdiv_to_uv(i_vort, i_div, ug, vg);

	SphereData_Physical tmpg1 = ug*fg;
	SphereData_Physical tmpg2 = vg*fg;

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);
	else
		op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	o_phi_t.spectral_set_zero();


#else

	/*
	 * This doesn't converge to the reference implementation!
	 */
	o_div_t = f*i_vort;
	o_vort_t = -f*i_div;
	o_phi_t.spectral_set_zero();

#endif
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




void SWE_Sphere_TS_lg_erk_lc_erk::run_timestep(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{

#if 0
		/*
		 * This routine can be used to validate the correct splitting with the geostrophic balanced test case
		 */
		SphereData_Spectral t1_phi(io_phi.sphereDataConfig);
		SphereData_Spectral t1_vort(io_phi.sphereDataConfig);
		SphereData_Spectral t1_div(io_phi.sphereDataConfig);

		SphereData_Spectral t2_phi(io_phi.sphereDataConfig);
		SphereData_Spectral t2_vort(io_phi.sphereDataConfig);
		SphereData_Spectral t2_div(io_phi.sphereDataConfig);

#if 1
		euler_timestep_update_lg(
				io_phi,
				io_vort,
				io_div,

				t1_phi,
				t1_vort,
				t1_div,

				i_simulation_timestamp
		);
		euler_timestep_update_lc(
				io_phi,
				io_vort,
				io_div,

				t2_phi,
				t2_vort,
				t2_div,

				i_simulation_timestamp
		);

#else
		euler_timestep_update(
				io_phi,
				io_vort,
				io_div,

				t1_phi,
				t1_vort,
				t1_div,

				i_simulation_timestamp
		);

		t2_phi.spectral_set_zero();
		t2_vort.spectral_set_zero();
		t2_div.spectral_set_zero();
#endif


		double diff_phi = (t1_phi + t2_phi).getSphereDataPhysical().physical_reduce_max_abs();
		double diff_vort = (t1_vort + t2_vort).getSphereDataPhysical().physical_reduce_max_abs();
		double diff_div = (t1_div + t2_div).getSphereDataPhysical().physical_reduce_max_abs();

		std::cout << "DIFF PHI: " << diff_phi << std::endl;
		std::cout << "DIFF VORT: " << diff_vort << std::endl;
		std::cout << "DIFF DIV: " << diff_div << std::endl;

		if (std::abs(diff_phi) > 1e-10)
			FatalError("PHI error too high");

		if (std::abs(diff_vort) > 1e-10)
			FatalError("VORT error too high");

		if (std::abs(diff_div) > 1e-10)
			FatalError("DIV error too high");
#endif


		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
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
				io_phi, io_vort, io_div,
				i_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else
	{
		FatalError("Not yet supported!");
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_lg_erk_lc_erk::setup(
		int i_order	///< order of RK time stepping method
//		int i_order2
)
{
	timestepping_order = i_order;
//	timestepping_order2 = i_order2;

	if (simVars.sim.sphere_use_fsphere)
	{
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = simVars.sim.sphere_fsphere_f0;
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

	f = fg;
}


SWE_Sphere_TS_lg_erk_lc_erk::SWE_Sphere_TS_lg_erk_lc_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_lg_erk_lc_erk::~SWE_Sphere_TS_lg_erk_lc_erk()
{
}

