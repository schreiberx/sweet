/*
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_l_erk_na_erk_vd.hpp"


bool PDESWESphereTS_l_erk_na_erk_vd::setup_auto(
		sweet::SphereOperators *io_ops
)
{
	return setup(io_ops, shackPDESWETimeDisc->timestepping_order, shackPDESWETimeDisc->timestepping_order2);
}

bool PDESWESphereTS_l_erk_na_erk_vd::setup(
		sweet::SphereOperators *io_ops,
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
		l_erk_split_vd = new PDESWESphereTS_ln_erk_split_vd;
		na_erk_split_vd = new PDESWESphereTS_ln_erk_split_vd;
	}

	l_erk_split_vd->setup(io_ops, timestepping_order, true, true, false, false, false);
	na_erk_split_vd->setup(io_ops, timestepping_order, false, false, true, false, false);

	return true;
}




bool PDESWESphereTS_l_erk_na_erk_vd::implementsTimesteppingMethod(const std::string &i_timestepping_method
									)
{
	timestepping_method = i_timestepping_method;
	return i_timestepping_method == "l_erk_na_erk_vd";
}


std::string PDESWESphereTS_l_erk_na_erk_vd::getIDString()
{
	return "l_erk_na_erk_vd";
}


void PDESWESphereTS_l_erk_na_erk_vd::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		l_erk_split_vd->runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_vd->runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{

		na_erk_split_vd->runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp
			);

		l_erk_split_vd->runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_vd->runTimestep(
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





PDESWESphereTS_l_erk_na_erk_vd::PDESWESphereTS_l_erk_na_erk_vd()
{
}



PDESWESphereTS_l_erk_na_erk_vd::~PDESWESphereTS_l_erk_na_erk_vd()
{
	delete l_erk_split_vd;
	delete na_erk_split_vd;
}

