/*
 * SWE_Sphere_TS_l_irk_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_irk_n_erk_ver01.hpp"




void SWE_Sphere_TS_l_irk_n_erk::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{

	if (timestepping_order == 1 && timestepping_order2 == 1)
	{
		if (version_id == 0)
		{
			// first order IRK for linear
			timestepping_l_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_l_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
	}
	else if (timestepping_order == 2 && timestepping_order2 >= 2)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_l_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! However, we only have autonomous simulations so far */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_irk.run_timestep_pert(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
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
void SWE_Sphere_TS_l_irk_n_erk::setup(
		int i_order,	///< order of RK time stepping method
		int i_order2,	///< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	version_id = i_version_id;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;
	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_l_irk.setup(
				1,
				timestep_size,	// Half time step size for linear implicit part (applied 2x)
				simVars.rexi.use_sphere_extended_modes
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_l_irk.setup(
					2,
					timestep_size*0.5,	// Half time step size for linear implicit part (applied 2x at start/end of TS)
					simVars.rexi.use_sphere_extended_modes,
					simVars.disc.timestepping_crank_nicolson_filter
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_irk.setup(
					2,
					timestep_size,	// Full time step size for linear implicit part (applied once at center)
					simVars.rexi.use_sphere_extended_modes,
					simVars.disc.timestepping_crank_nicolson_filter
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



void SWE_Sphere_TS_l_irk_n_erk::setup_auto()
{
	if (
		simVars.disc.timestepping_method == "l_irk_n_erk" ||
		simVars.disc.timestepping_method == "l_irk_n_erk_ver0" ||
		simVars.disc.timestepping_method == "l_cn_n_erk" ||
		simVars.disc.timestepping_method == "l_cn_n_erk_ver0"
	)
		setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2, 0);
	else
		setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2, 1);
}


SWE_Sphere_TS_l_irk_n_erk::SWE_Sphere_TS_l_irk_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_l_irk(simVars, op),
		timestepping_l_erk_n_erk(simVars, op),
		version_id(0),
		timestepping_order(-1)
{

}



SWE_Sphere_TS_l_irk_n_erk::~SWE_Sphere_TS_l_irk_n_erk()
{
}

