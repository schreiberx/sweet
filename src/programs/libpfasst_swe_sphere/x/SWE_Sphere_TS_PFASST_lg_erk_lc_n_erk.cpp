/*
 * SWE_Sphere_TS_PFASST_lg_erk_lf_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk.hpp"



void SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_linear(
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

	double avgphi = simVars.sim.gravitation*simVars.sim.h0;

	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	SphereData_Physical vrtg = i_vort.toPhys();
	SphereData_Physical divg = i_div.toPhys();
	op.robert_vrtdiv_to_uv(i_vort, i_div, ug, vg);
	SphereData_Physical phig = i_phi.toPhys();

	// No Coriolis here
#if 0
	SphereData_Physical tmpg1 = ug*(/*vrtg+*/fg);
	SphereData_Physical tmpg2 = vg*(/*vrtg+*/fg);

	op.robert_uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	SphereData_Physical tmpg = o_div_t.toPhys();
#else
	o_vort_t.spectral_set_zero();
#endif

	/*
	tmpg1 = ug*phig;
	tmpg2 = vg*phig;
	*/
	SphereData_Physical tmpg1 = ug*avgphi;
	SphereData_Physical tmpg2 = vg*avgphi;

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	op.robert_uv_to_vrtdiv(tmpg1, tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	tmpspec.loadSphereDataPhysical(phig/*+0.5*(ug*ug+vg*vg)*/);

//	o_div_t += -op.laplace(tmpspec);
	o_div_t = -op.laplace(tmpspec);
}


void SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
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

	double avgphi = simVars.sim.gravitation*simVars.sim.h0;

	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	SphereData_Physical vrtg = i_vort.toPhys();
	SphereData_Physical divg = i_div.toPhys();
	op.robert_vrtdiv_to_uv(i_vort, i_div, ug, vg);
	SphereData_Physical phig = i_phi.toPhys();

	SphereData_Physical tmpg1 = ug*(vrtg+fg);
	SphereData_Physical tmpg2 = vg*(vrtg+fg);

	op.robert_uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	SphereData_Physical tmpg = o_div_t.toPhys();

	tmpg1 = ug*(phig-avgphi);
	tmpg2 = vg*(phig-avgphi);

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	op.robert_uv_to_vrtdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	tmpspec.loadSphereDataPhysical(/*phig+*/0.5*(ug*ug+vg*vg));

	o_div_t += -op.laplace(tmpspec);
}



/**
 * This routine is used by other time step implementations
 */
void SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
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

	euler_timestep_update_lc_n(
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




void SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::run_timestep_nonpert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{

		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);


			// FULL time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETError("Invalid verison id");
		}
	}
	else
	{
		SWEETError("Not yet supported!");
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::setup(
		int i_order,	///< order of RK time stepping method
		int i_version_id
)
{
	timestepping_order = i_order;

	version_id = i_version_id;

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



SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_order(-1),
		fg(i_op.sphereDataConfig)
{
}



SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk::~SWE_Sphere_TS_PFASST_lg_erk_lc_n_erk()
{
}

