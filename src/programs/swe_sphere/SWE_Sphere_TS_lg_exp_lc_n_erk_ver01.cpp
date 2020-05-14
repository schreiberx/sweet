/*
 * SWE_Sphere_TS_lg_irk_lf_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_exp_lc_n_erk_ver01.hpp"



void SWE_Sphere_TS_lg_exp_lc_n_erk::run_timestep_pert(
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


void SWE_Sphere_TS_lg_exp_lc_n_erk::run_timestep_nonpert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		SphereData_Spectral tmp_phi = io_phi;
		SphereData_Spectral tmp_vort = io_vort;
		SphereData_Spectral tmp_div = io_div;

		// first order IRK for linear
		timestepping_lg_rexi.run_timestep_nonpert(
				io_phi, io_vort, io_div,
				i_dt,
				i_simulation_timestamp
			);

		SphereData_Spectral phi_dt(io_phi.sphereDataConfig);
		SphereData_Spectral vort_dt(io_vort.sphereDataConfig);
		SphereData_Spectral div_dt(io_div.sphereDataConfig);

		// first order explicit for non-linear
		timestepping_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				tmp_phi, tmp_vort, tmp_div,
				phi_dt, vort_dt, div_dt,
				i_simulation_timestamp
			);


		io_phi += i_dt*phi_dt;
		io_vort += i_dt*vort_dt;
		io_div += i_dt*div_dt;
	}
	else if (timestepping_order == 2 || timestepping_order == 4)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_lg_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_n_erk,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute Euler time step updates
					io_phi, io_vort, io_div,
					i_dt,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_n_erk,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_n_erk,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
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
void SWE_Sphere_TS_lg_exp_lc_n_erk::setup(
		REXI_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size,
		int i_version_id
)
{
	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order2 == 1)
	{
		timestepping_lg_rexi.setup(
				i_rexiSimVars,
				"phi0",
				i_timestep_size,
				false,
				true
			);
	}
	else if (timestepping_order2 == 2 || timestepping_order2 == 4)
	{
		if (version_id == 0)
		{
			timestepping_lg_rexi.setup(
					i_rexiSimVars,
					"phi0",
					i_timestep_size*0.5,
					false,
					true
			);
		}
		else if (version_id == 1)
		{
			timestepping_lg_rexi.setup(
					i_rexiSimVars,
					"phi0",
					i_timestep_size,
					false,
					true
			);
		}
		else
		{
			SWEETError("Invalid version id");
		}
	}
	else
	{
		SWEETError("Invalid timestepping order");
	}


	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_lg_erk_lc_n_erk.setup(1, -1);
}



SWE_Sphere_TS_lg_exp_lc_n_erk::SWE_Sphere_TS_lg_exp_lc_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_order(-1),
		timestepping_lg_rexi(simVars, op),
		timestepping_lg_erk_lc_n_erk(simVars, op)
{
}



SWE_Sphere_TS_lg_exp_lc_n_erk::~SWE_Sphere_TS_lg_exp_lc_n_erk()
{
}

