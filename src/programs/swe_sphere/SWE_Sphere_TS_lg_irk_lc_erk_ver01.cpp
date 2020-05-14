/*
 * SWE_Sphere_TS_lg_irk_lf_n_erk.cpp
 *
 *  Created on: 11 Nov 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_irk_lc_erk_ver01.hpp"



void SWE_Sphere_TS_lg_irk_lc_erk::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		if (version_id == 0)
		{
			// first order IRK for linear
			timestepping_lg_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);


			// first order explicit for non-linear
			timestepping_lg_erk_lc_erk.euler_timestep_update_lc(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else if (version_id == 0)
		{
			// first order explicit for non-linear
			timestepping_lg_erk_lc_erk.euler_timestep_update_lc(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_lg_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_lg_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_erk,
					&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_erk,
					&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_erk,
					&SWE_Sphere_TS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
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
void SWE_Sphere_TS_lg_irk_lc_erk::setup(
		int i_timestepping_order,	///< order of RK time stepping method
		int i_version_id
)
{
	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_lg_irk.setup(
				1,
				timestep_size*0.5,
				0
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_lg_irk.setup(
					2,
					timestep_size*0.5,
					0,
					simVars.disc.timestepping_crank_nicolson_filter
			);
		}
		else if (version_id == 1)
		{
			timestepping_lg_irk.setup(
					2,
					timestep_size,
					0,
					simVars.disc.timestepping_crank_nicolson_filter
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
	timestepping_lg_erk_lc_erk.setup(1);
}


/*
 * Setup
 */
void SWE_Sphere_TS_lg_irk_lc_erk::setup_auto()
{
	if (simVars.disc.timestepping_method == "lg_irk_lc_erk_ver1")
		setup(simVars.disc.timestepping_order, 1);
	else
		setup(simVars.disc.timestepping_order, 0);
}



SWE_Sphere_TS_lg_irk_lc_erk::SWE_Sphere_TS_lg_irk_lc_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_order(-1),
		timestepping_lg_irk(simVars, op),
		timestepping_lg_erk_lc_erk(simVars, op)
{
}



SWE_Sphere_TS_lg_irk_lc_erk::~SWE_Sphere_TS_lg_irk_lc_erk()
{
}

