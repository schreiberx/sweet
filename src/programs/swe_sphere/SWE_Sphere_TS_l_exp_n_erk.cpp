/*
 * SWE_Sphere_TS_l_rexi_n_erk.cpp
 *
 *  Created on: 1 Oct 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_exp_n_erk_ver01.hpp"



void SWE_Sphere_TS_l_exp_n_erk::run_timestep_pert(
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




void SWE_Sphere_TS_l_exp_n_erk::run_timestep_nonpert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order2 == 1)
	{
		if (version_id == 0)
		{
			// first order REXI for linear part
			timestepping_l_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_n(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
		else if (version_id == 1)
		{
			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_n(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order REXI for linear part
			timestepping_l_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETError("Invalid version id");
		}
	}
	else if (timestepping_order2 == 2 || timestepping_order2 == 4)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_l_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);

		}
		else if (version_id == 1)
		{

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_rexi.run_timestep_nonpert(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETError("Invalid version id");
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
void SWE_Sphere_TS_l_exp_n_erk::setup(
		REXI_SimulationVariables &i_rexiSimVars,
		int i_order,	///< order of RK time stepping method
		int i_order2,	///< order of RK time stepping method of non-linear parts
		double i_timestep_size,
		bool i_use_f_sphere,
		int i_version_id
)
{
	version_id = i_version_id;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	if (timestepping_order != timestepping_order2)
		SWEETError("Mismatch of orders, should be equal");

	timestep_size = i_timestep_size;

	if (timestepping_order2 == 1)
	{
		timestepping_l_rexi.setup(
				i_rexiSimVars,
				"phi0",
				i_timestep_size,
				i_use_f_sphere,
				false
		);
	}
	else if (timestepping_order2 == 2 || timestepping_order2 == 4)
	{
		if (version_id == 0)
		{
			timestepping_l_rexi.setup(
					i_rexiSimVars,
					"phi0",
					i_timestep_size*0.5,
					i_use_f_sphere,
					false
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_rexi.setup(
					i_rexiSimVars,
					"phi0",
					i_timestep_size,
					i_use_f_sphere,
					false
			);
		}
		else
		{
			SWEETError("Invalid version");
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
	timestepping_l_erk_n_erk.setup(1, 1);
}



void SWE_Sphere_TS_l_exp_n_erk::setup_auto()
{
	int version_id = 0;

	if (simVars.disc.timestepping_method == "l_exp_n_erk_ver1")
		version_id = 1;

	setup(
			simVars.rexi,
			simVars.disc.timestepping_order,
			simVars.disc.timestepping_order2,
			simVars.timecontrol.current_timestep_size,
			simVars.sim.sphere_use_fsphere,
			version_id
		);
}



SWE_Sphere_TS_l_exp_n_erk::SWE_Sphere_TS_l_exp_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_l_rexi(simVars, op),
		timestepping_l_erk_n_erk(simVars, op),
		version_id(0),
		timestepping_order(-1)
{
}



SWE_Sphere_TS_l_exp_n_erk::~SWE_Sphere_TS_l_exp_n_erk()
{
}

