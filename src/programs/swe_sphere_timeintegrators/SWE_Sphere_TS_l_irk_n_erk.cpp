/*
 * SWE_Sphere_TS_l_irk_n_erk.cpp
 *
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_l_irk_n_erk.hpp"


bool SWE_Sphere_TS_l_irk_n_erk::implements_timestepping_method(const std::string &i_timestepping_method)
{
	if (
		i_timestepping_method == "l_irk_n_erk"	|| i_timestepping_method == "l_irk_n_erk_ver0"	||
		i_timestepping_method == "l_cn_n_erk"	|| i_timestepping_method == "l_cn_n_erk_ver0"		||
		i_timestepping_method == "l_irk_n_erk_ver1"
	)
		return true;

	return false;
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



std::string SWE_Sphere_TS_l_irk_n_erk::string_id()
{
	std::string s = "l_irk_n_erk_ver";

	if (version_id == 0)
		s += "0";
	else if (version_id == 1)
		s += "1";
	else
		SWEETError("Version ID");

	return s;
}


void SWE_Sphere_TS_l_irk_n_erk::run_timestep(
		SphereData_Spectral &io_phi,
		SphereData_Spectral &io_vrt,
		SphereData_Spectral &io_div,

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{

	if (timestepping_order == 1)
	{
		if (version_id == 0)
		{
			// first order IRK for linear
			timestepping_l_irk.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_l_irk.run_timestep(
					io_phi, io_vrt, io_div,
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
			timestepping_l_irk.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_irk.run_timestep(
					io_phi, io_vrt, io_div,
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
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_irk.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					&timestepping_l_erk_n_erk,
					&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi, io_vrt, io_div,
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
	if (i_order2 < 0)
		i_order2 = i_order;

	if (i_order != i_order2)
		SWEETError("Orders of 1st and 2nd one must match");

	version_id = i_version_id;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;
	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_l_irk.setup(
			1,
			timestep_size
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_l_irk.setup(
					2,
					timestep_size*0.5,
					simVars.disc.timestepping_crank_nicolson_filter,
					false
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_irk.setup(
					2,
					timestep_size,
					simVars.disc.timestepping_crank_nicolson_filter,
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

