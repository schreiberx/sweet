/*
 * SWE_Sphere_TS_l_erk_na_erk_uv.cpp
 *
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_erk_na_erk_uv.hpp"



bool SWE_Sphere_TS_l_erk_na_erk_uv::implements_timestepping_method(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "l_erk_na_erk_uv";
}


std::string SWE_Sphere_TS_l_erk_na_erk_uv::string_id()
{
	return "l_erk_na_erk_uv";
}


void SWE_Sphere_TS_l_erk_na_erk_uv::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		l_erk_split_uv->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_uv->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{

		na_erk_split_uv->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp
			);

		l_erk_split_uv->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_uv->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp
			);
	}
	else
	{
		SWEETError("Not yet supported!");
	}
}




/*
 * Setup
 */
void SWE_Sphere_TS_l_erk_na_erk_uv::setup(
		int i_order,	///< order of RK time stepping method for non-linear parts
		int i_order2	///< order of RK time stepping method for non-linear parts
)
{
	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	if (timestepping_order != timestepping_order2)
		SWEETError("Timestepping orders must match");

	if (l_erk_split_uv == nullptr)
	{
		l_erk_split_uv = new SWE_Sphere_TS_ln_erk_split_uv(simVars, ops);
		na_erk_split_uv = new SWE_Sphere_TS_ln_erk_split_uv(simVars, ops);
	}

	l_erk_split_uv->setup(timestepping_order, true, true, false, false, false);
	na_erk_split_uv->setup(timestepping_order, false, false, true, false, false);
}



void SWE_Sphere_TS_l_erk_na_erk_uv::setup_auto()
{
	setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2);
}



SWE_Sphere_TS_l_erk_na_erk_uv::SWE_Sphere_TS_l_erk_na_erk_uv(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_ops
)	:
		simVars(i_simVars),
		ops(i_ops)
{
}



SWE_Sphere_TS_l_erk_na_erk_uv::~SWE_Sphere_TS_l_erk_na_erk_uv()
{
	delete l_erk_split_uv;
	delete na_erk_split_uv;
}

