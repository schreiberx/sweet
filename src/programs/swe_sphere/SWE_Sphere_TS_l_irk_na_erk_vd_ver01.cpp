/*
 * SWE_Sphere_TS_l_irk_na_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_irk_na_erk_vd_ver01.hpp"




void SWE_Sphere_TS_l_irk_na_erk_vd::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation*simVars.sim.h0;
	if (timestepping_order == 1 && timestepping_order2 == 1)
	{
		if (version_id == 0)
		{
			// first order IRK for linear
			timestepping_l_irk.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_ln_erk_split_vd.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_ln_erk_split_vd.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_l_irk.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
	}
	else if (timestepping_order == 2 && timestepping_order2 >= 2)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_l_cn.run_timestep_nonpert(
					io_phi_pert, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_ln_erk_split_vd.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_cn.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_ln_erk_split_vd.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_cn.run_timestep_nonpert(
					io_phi_pert, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_ln_erk_split_vd.run_timestep_pert(
					io_phi_pert, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);
		}
		else
		{
			FatalError("Invalid version id");
		}
	}
	else
	{
		FatalError("Not yet supported!");
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_irk_na_erk_vd::setup(
		int i_order,	///< order of RK time stepping method
		int i_order2,	///< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	i_order = i_order2;
	//if (i_order != i_order2)
	//	FatalError("Orders of 1st and 2nd one must match");

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
			timestepping_l_cn.setup(
					simVars.disc.timestepping_crank_nicolson_filter,
					timestep_size*0.5,	// Half time step size for linear implicit part (applied 2x at start/end of TS)
					simVars.rexi.use_sphere_extended_modes
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_cn.setup(
					simVars.disc.timestepping_crank_nicolson_filter,
					timestep_size,	// Full time step size for linear implicit part (applied once at center)
					simVars.rexi.use_sphere_extended_modes
			);
		}
		else
		{
			FatalError("Invalid version");
		}
	}
	else
	{
		FatalError("Invalid timestepping order");
	}


	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_ln_erk_split_vd.setup(i_order2, false, false, true, false);
}



void SWE_Sphere_TS_l_irk_na_erk_vd::setup_auto()
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


SWE_Sphere_TS_l_irk_na_erk_vd::SWE_Sphere_TS_l_irk_na_erk_vd(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_l_irk(simVars, op),
		timestepping_l_cn(simVars, op),
		timestepping_ln_erk_split_vd(simVars, op),
		version_id(0),
		timestepping_order(-1)
{

}



SWE_Sphere_TS_l_irk_na_erk_vd::~SWE_Sphere_TS_l_irk_na_erk_vd()
{
}

