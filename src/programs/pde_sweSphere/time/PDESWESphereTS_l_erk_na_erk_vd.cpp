/*
 * PDESWESphereTS_l_erk_na_erk_vd.cpp
 *
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_l_erk_na_erk_vd.hpp"



bool PDESWESphereTS_l_erk_na_erk_vd::implements_timestepping_method(const std::string &i_timestepping_method
									)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackDict.disc.timestepping_order;
	timestepping_order2 = shackDict.disc.timestepping_order2;
	return i_timestepping_method == "l_erk_na_erk_vd";
}


std::string PDESWESphereTS_l_erk_na_erk_vd::string_id()
{
	return "l_erk_na_erk_vd";
}


void PDESWESphereTS_l_erk_na_erk_vd::run_timestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
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
void PDESWESphereTS_l_erk_na_erk_vd::setup(
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
		l_erk_split_vd = new PDESWESphereTS_ln_erk_split_vd(shackDict, ops);
		na_erk_split_vd = new PDESWESphereTS_ln_erk_split_vd(shackDict, ops);
	}

	l_erk_split_vd->setup(timestepping_order, true, true, false, false, false);
	na_erk_split_vd->setup(timestepping_order, false, false, true, false, false);
}



void PDESWESphereTS_l_erk_na_erk_vd::setup_auto()
{
	setup(timestepping_order, shackDict.disc.timestepping_order2);
}



PDESWESphereTS_l_erk_na_erk_vd::PDESWESphereTS_l_erk_na_erk_vd(
		sweet::ShackDictionary &i_shackDict,
		sweet::SphereOperators &i_ops
)	:
		shackDict(i_shackDict),
		ops(i_ops)
{
}



PDESWESphereTS_l_erk_na_erk_vd::~PDESWESphereTS_l_erk_na_erk_vd()
{
	delete l_erk_split_vd;
	delete na_erk_split_vd;
}

