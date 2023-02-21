/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_irk_lc_n_erk_ver01.hpp"



void SWE_Sphere_TS_lg_irk_lc_n_erk::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		if (version_id == 0)
		{
			// first order IRK for linear
			timestepping_lg_irk.run_timestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			SphereData_Spectral phi_fixed_dt(io_phi_pert.sphereDataConfig), vrt_dt(io_vrt.sphereDataConfig), div_dt(io_div.sphereDataConfig);

			// first order explicit for non-linear
			timestepping_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
					io_phi_pert, io_vrt, io_div,
					phi_fixed_dt, vrt_dt, div_dt,
					i_simulation_timestamp
				);

			io_phi_pert += i_fixed_dt*phi_fixed_dt;
			io_vrt += i_fixed_dt*vrt_dt;
			io_div += i_fixed_dt*div_dt;
		}
		else
		{
			SphereData_Spectral phi_fixed_dt(io_phi_pert.sphereDataConfig);
			SphereData_Spectral vrt_dt(io_vrt.sphereDataConfig);
			SphereData_Spectral div_dt(io_div.sphereDataConfig);

			timestepping_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
					io_phi_pert, io_vrt, io_div,
					phi_fixed_dt, vrt_dt, div_dt,
					i_simulation_timestamp
				);

			io_phi_pert += i_fixed_dt*phi_fixed_dt;
			io_vrt += i_fixed_dt*vrt_dt;
			io_div += i_fixed_dt*div_dt;

			// first order IRK for linear
			timestepping_lg_irk.run_timestep(
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
			timestepping_lg_irk.run_timestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_n_erk,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_irk.run_timestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_n_erk,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_irk.run_timestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_lg_erk_lc_n_erk,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
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
void SWE_Sphere_TS_lg_irk_lc_n_erk::setup(
		int i_timestepping_order,	///< order of RK time stepping method
		int i_timestepping_order2,
		int i_version_id
)
{
	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;
	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_lg_irk.setup(
				1,
				timestep_size*0.5
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_lg_irk.setup(
					2,
					timestep_size*0.5
			);
		}
		else if (version_id == 1)
		{
			timestepping_lg_irk.setup(
					2,
					timestep_size
			);
		}
		else
		{
			SWEETError("Invalid version id");
		}
	}
	else
	{
		SWEETError("Invalid time stepping order");
	}


	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_lg_erk_lc_n_erk.setup(1, -1);
}



void SWE_Sphere_TS_lg_irk_lc_n_erk::setup_auto()
{
	int version = 0;
	if (timestepping_method == "lg_irk_lc_n_erk_ver1")
		version = 1;

	setup(
			timestepping_order,
			timestepping_order2,
			version
		);
}



SWE_Sphere_TS_lg_irk_lc_n_erk::SWE_Sphere_TS_lg_irk_lc_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_order(-1),
		timestepping_lg_irk(simVars, op),
		timestepping_lg_erk_lc_n_erk(simVars, op)
{
}



SWE_Sphere_TS_lg_irk_lc_n_erk::~SWE_Sphere_TS_lg_irk_lc_n_erk()
{
}

