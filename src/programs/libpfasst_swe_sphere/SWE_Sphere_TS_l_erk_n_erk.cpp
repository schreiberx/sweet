/*
 * SWE_Sphere_TS_l_erk_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_erk_n_erk.hpp"



void SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear(
		const SphereData &i_phi,	///< prognostic variables
		const SphereData &i_vort,	///< prognostic variables
		const SphereData &i_div,	///< prognostic variables

		SphereData &o_phi_t,	///< time updates
		SphereData &o_vort_t,	///< time updates
		SphereData &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	double avgphi = simVars.sim.gravitation*simVars.sim.h0;

	SphereDataPhysical ug(i_phi.sphereDataConfig);
	SphereDataPhysical vg(i_phi.sphereDataConfig);

	SphereDataPhysical vrtg = i_vort.getSphereDataPhysical();
	SphereDataPhysical divg = i_div.getSphereDataPhysical();
	op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	SphereDataPhysical phig = i_phi.getSphereDataPhysical();

	SphereDataPhysical tmpg1 = ug*(/*vrtg+*/fg);
	SphereDataPhysical tmpg2 = vg*(/*vrtg+*/fg);

	op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

	/*
	tmpg1 = ug*phig;
	tmpg2 = vg*phig;
	*/
	tmpg1 = ug*avgphi;
	tmpg2 = vg*avgphi;

	SphereData tmpspec(i_phi.sphereDataConfig);
	op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	tmpspec = (phig/*+0.5*(ug*ug+vg*vg)*/);
	tmpspec.request_data_spectral();
	o_div_t += -op.laplace(tmpspec);
}


void SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_n(
		const SphereData &i_phi,	///< prognostic variables
		const SphereData &i_vort,	///< prognostic variables
		const SphereData &i_div,	///< prognostic variables

		SphereData &o_phi_dt,	///< time updates
		SphereData &o_vort_dt,	///< time updates
		SphereData &o_div_dt,	///< time updates

		double i_simulation_timestamp
)
{

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	double avgphi = simVars.sim.gravitation*simVars.sim.h0;

	SphereDataPhysical ug(i_phi.sphereDataConfig);
	SphereDataPhysical vg(i_phi.sphereDataConfig);

	SphereDataPhysical vrtg = i_vort.getSphereDataPhysical();
	SphereDataPhysical divg = i_div.getSphereDataPhysical();
	op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	SphereDataPhysical phig = i_phi.getSphereDataPhysical();

	SphereDataPhysical tmpg1 = ug*(vrtg/*+fg*/);
	SphereDataPhysical tmpg2 = vg*(vrtg/*+fg*/);

	op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_dt, o_vort_dt);

	o_vort_dt *= -1.0;

	SphereDataPhysical tmpg = o_div_dt.getSphereDataPhysical();

	tmpg1 = ug*(phig-avgphi);
	tmpg2 = vg*(phig-avgphi);

	SphereData tmpspec(i_phi.sphereDataConfig);
	op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_dt);

	o_phi_dt *= -1.0;

	tmpspec = (/*phig+*/0.5*(ug*ug+vg*vg));
	tmpspec.request_data_spectral();
	o_div_dt += -op.laplace(tmpspec);
}



/**
 * This routine is used by other time step implementations
 */
void SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_n(
		SphereData &io_phi,		///< prognostic variables
		SphereData &io_vort,	///< prognostic variables
		SphereData &io_div,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	SphereData tmp_phi(io_phi.sphereDataConfig);
	SphereData tmp_vort(io_vort.sphereDataConfig);
	SphereData tmp_div(io_div.sphereDataConfig);

	euler_timestep_update_n(
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




//#if 0
/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update(
		const SphereData &i_phi,	///< prognostic variables
		const SphereData &i_vort,	///< prognostic variables
		const SphereData &i_div,	///< prognostic variables

		SphereData &o_phi_t,	///< time updates
		SphereData &o_vort_t,	///< time updates
		SphereData &o_div_t,	///< time updates

		//double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	/*
	 * TIME STEP SIZE
	 */
        //if (i_fixed_dt <= 0)
	//	FatalError("Only fixed time step size allowed");

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	SphereDataPhysical ug(i_phi.sphereDataConfig);
	SphereDataPhysical vg(i_phi.sphereDataConfig);

	SphereDataPhysical vrtg = i_vort.getSphereDataPhysical();
	SphereDataPhysical divg = i_div.getSphereDataPhysical();
	op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	SphereDataPhysical phig = i_phi.getSphereDataPhysical();

	SphereDataPhysical tmpg1 = ug*(vrtg+fg);
	SphereDataPhysical tmpg2 = vg*(vrtg+fg);

	op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	SphereData tmpspec(i_phi.sphereDataConfig);
	op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	tmpspec = (phig+0.5*(ug*ug+vg*vg));
	tmpspec.request_data_spectral();
	o_div_t += -op.laplace(tmpspec);
}
//#endif


void SWE_Sphere_TS_l_erk_n_erk::run_timestep(
		SphereData &io_phi,		///< prognostic variables
		SphereData &io_vort,	///< prognostic variables
		SphereData &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_n,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order2,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_n,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
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
void SWE_Sphere_TS_l_erk_n_erk::setup(
		int i_order,	///< order of RK time stepping method for non-linear parts
		int i_order2	///< order of RK time stepping method for non-linear parts
)
{
	timestepping_order = i_order;
	timestepping_order2 = i_order2;

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


SWE_Sphere_TS_l_erk_n_erk::SWE_Sphere_TS_l_erk_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		fg(i_op.sphereDataConfig)
{
	setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2);
}



SWE_Sphere_TS_l_erk_n_erk::~SWE_Sphere_TS_l_erk_n_erk()
{
}

