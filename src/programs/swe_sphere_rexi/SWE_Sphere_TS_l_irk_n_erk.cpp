/*
 * SWE_Sphere_TS_l_irk_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Sphere_TS_l_irk_n_erk.hpp"



void SWE_Sphere_TS_l_irk_n_erk::run_timestep(
		SphereData &io_phi,		///< prognostic variables
		SphereData &io_vort,	///< prognostic variables
		SphereData &io_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		// first order IRK for linear
		timestepping_l_irk.run_timestep(
				io_phi, io_vort, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		// first order explicit for non-linear
		timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
				io_phi, io_vort, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
#if 1
		// HALF time step for linear part
		timestepping_l_cn.run_timestep(
				io_phi, io_vort, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				&timestepping_l_erk_n_erk,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_fixed_dt,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_l_cn.run_timestep(
				io_phi, io_vort, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
			);
#else
		// HALF time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				&timestepping_l_erk_n_erk,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for linear part
		timestepping_l_cn.run_timestep(
				io_phi, io_vort, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		// HALF time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				&timestepping_l_erk_n_erk,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

#endif
	}
	else
	{
		FatalError("Not yet supported!");
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_irk_n_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_l_irk.setup(
				1,
				timestep_size*0.5,	// Half time step size for linear implicit part (applied 2x)
				simVars.rexi.use_sphere_extended_modes
		);
	}
	else if (timestepping_order == 2)
	{
		timestepping_l_cn.setup(
				simVars.disc.crank_nicolson_filter,
				timestep_size*0.5,	// Half time step size for linear implicit part (applied 2x)
				simVars.rexi.use_sphere_extended_modes
		);
	}
	else
	{
		FatalError("Invalid timestepping method");
	}


	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_l_erk_n_erk.setup(1);
}



SWE_Sphere_TS_l_irk_n_erk::SWE_Sphere_TS_l_irk_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_l_irk(simVars, op),
		timestepping_l_cn(simVars, op),
		timestepping_l_erk_n_erk(simVars, op)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Sphere_TS_l_irk_n_erk::~SWE_Sphere_TS_l_irk_n_erk()
{
}

