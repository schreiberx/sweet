/*
 * SWE_Sphere_TS_l_irk_na_erk_vd.cpp
 *
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_irk_na_erk_vd.hpp"


bool SWE_Sphere_TS_l_irk_na_erk_vd::implements_timestepping_method(const std::string &i_timestepping_method
//////#if SWEET_PARAREAL
//////									,
//////									int &i_timestepping_order,
//////									int &i_timestepping_order2
//////#endif
									)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;
////////#if SWEET_PARAREAL
////////	timestepping_order = i_timestepping_order;
////////	timestepping_order2 = i_timestepping_order2;
////////#endif
	if (
		i_timestepping_method == "l_irk_na_erk_vd" || i_timestepping_method == "l_irk_na_erk_vd_ver0" ||
		i_timestepping_method == "l_irk_na_erk_vd_ver1"
	)
		return true;

	return false;
}



void SWE_Sphere_TS_l_irk_na_erk_vd::setup_auto()
{
	if (
		timestepping_method == "l_irk_na_erk_vd" ||
		timestepping_method == "l_irk_na_erk_vd_ver0"
	)
		setup(timestepping_order, timestepping_order2, 0);
	else if (
			timestepping_method == "l_irk_na_erk_vd_ver1"
		)
		setup(timestepping_order, timestepping_order2, 1);
	else
		SWEETError("Not implemented");
}



std::string SWE_Sphere_TS_l_irk_na_erk_vd::string_id()
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



void SWE_Sphere_TS_l_irk_na_erk_vd::run_timestep(
		SphereData_Spectral &io_phi,
		SphereData_Spectral &io_vrt,
		SphereData_Spectral &io_div,

		double i_fixed_dt,
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
			timestepping_na_erk_split_vd.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_na_erk_split_vd.run_timestep(
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
			timestepping_na_erk_split_vd.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_irk.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_na_erk_split_vd.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_irk.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_na_erk_split_vd.run_timestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
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
void SWE_Sphere_TS_l_irk_na_erk_vd::setup(
		int i_timestepping_order,	///< order of RK time stepping method
		int i_timestepping_order2,	///< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	if (i_timestepping_order2 < 0)
		i_timestepping_order2 = i_timestepping_order;

	if (i_timestepping_order != i_timestepping_order2)
		SWEETError("Orders of 1st and 2nd one must match");

	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;
	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_l_irk.setup(
			1,
			timestep_size
		);

		timestepping_na_erk_split_vd.setup(1, false, false, true, false, false);
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

		timestepping_na_erk_split_vd.setup(2, false, false, true, false, false);
	}
	else
	{
		SWEETError("Invalid timestepping order");
	}
}



SWE_Sphere_TS_l_irk_na_erk_vd::SWE_Sphere_TS_l_irk_na_erk_vd(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_l_irk(simVars, op),
		timestepping_na_erk_split_vd(simVars, op),
		version_id(0),
		timestepping_order(-1)
{

}



SWE_Sphere_TS_l_irk_na_erk_vd::~SWE_Sphere_TS_l_irk_na_erk_vd()
{
}

