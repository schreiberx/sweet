/*
 * SWE_Sphere_TS_l_erk_na_erk_vd.cpp
 *
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_l_erk_na_erk_vd.hpp"



bool SWE_Sphere_TS_l_erk_na_erk_vd::implements_timestepping_method(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "l_erk_na_erk_vd";
}


std::string SWE_Sphere_TS_l_erk_na_erk_vd::string_id()
{
	return "l_erk_na_erk_vd";
}


void SWE_Sphere_TS_l_erk_na_erk_vd::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		l_erk_split_vd->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_vd->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{

		na_erk_split_vd->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp
			);

		l_erk_split_vd->run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_vd->run_timestep(
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
void SWE_Sphere_TS_l_erk_na_erk_vd::setup(
		int i_order,	///< order of RK time stepping method for non-linear parts
		int i_order2	///< order of RK time stepping method for non-linear parts
)
{
	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	if (timestepping_order != timestepping_order2)
		SWEETError("Timestepping orders must match");

	if (l_erk_split_vd == nullptr)
	{
		l_erk_split_vd = new SWE_Sphere_TS_ln_erk_split_vd(simVars, ops);
		na_erk_split_vd = new SWE_Sphere_TS_ln_erk_split_vd(simVars, ops);
	}

	l_erk_split_vd->setup(timestepping_order, true, true, false, false, false);
	na_erk_split_vd->setup(timestepping_order, false, false, true, false, false);
}



void SWE_Sphere_TS_l_erk_na_erk_vd::setup_auto()
{
	setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2);
}



SWE_Sphere_TS_l_erk_na_erk_vd::SWE_Sphere_TS_l_erk_na_erk_vd(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_ops
)	:
		simVars(i_simVars),
		ops(i_ops)
{
}



SWE_Sphere_TS_l_erk_na_erk_vd::~SWE_Sphere_TS_l_erk_na_erk_vd()
{
	delete l_erk_split_vd;
	delete na_erk_split_vd;
}

