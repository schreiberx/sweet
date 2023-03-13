/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_irk_lc_na_erk_vd.hpp"



bool PDESWESphereTS_lg_irk_lc_na_erk_vd::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (
		i_timestepping_method == "lg_irk_lc_na_erk_vd" || i_timestepping_method == "lg_irk_lc_na_erk_vd_ver0" ||
		i_timestepping_method == "lg_irk_lc_na_erk_vd_ver1"
	)
		return true;

	return false;
}



bool PDESWESphereTS_lg_irk_lc_na_erk_vd::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (
		timestepping_method == "lg_irk_lc_na_erk_vd" ||
		timestepping_method == "lg_irk_lc_na_erk_vd_ver0"
	)
	{
		return setup(io_ops, timestepping_order, timestepping_order2, 0);
	}
	else if (
			timestepping_method == "lg_irk_lc_na_erk_vd_ver1"
		)
	{
		return setup(io_ops, timestepping_order,timestepping_order2, 1);
	}
	else
	{
		SWEETError("Not implemented");
	}

	return false;
}


bool PDESWESphereTS_lg_irk_lc_na_erk_vd::setup(
		sweet::SphereOperators *io_ops,
		int i_timestepping_order,	///< order of RK time stepping method
		int i_timestepping_order2,	///< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	ops = io_ops;

	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	if (i_timestepping_order2 < 0)
		i_timestepping_order2 = i_timestepping_order;

	if (i_timestepping_order != i_timestepping_order2)
		SWEETError("Orders of 1st and 2nd one must match");

	timestep_size = shackTimestepControl->current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_lg_irk.setup(
			ops,
			1,
			timestep_size
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_lg_irk.setup(
				ops,
				2,
				timestep_size*0.5
			);
		}
		else if (version_id == 1)
		{
			timestepping_lg_irk.setup(
				ops,
				2,
				timestep_size
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


	// Only NA part!
	timestepping_ln_erk_split_vd.setup_main(
			ops,
			i_timestepping_order2,
			false, true, true, false, false
		);

	return true;
}




std::string PDESWESphereTS_lg_irk_lc_na_erk_vd::getIDString()
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



void PDESWESphereTS_lg_irk_lc_na_erk_vd::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		if (version_id == 0)
		{
			// first order IRK for linear
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
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




PDESWESphereTS_lg_irk_lc_na_erk_vd::PDESWESphereTS_lg_irk_lc_na_erk_vd()
{

}



PDESWESphereTS_lg_irk_lc_na_erk_vd::~PDESWESphereTS_lg_irk_lc_na_erk_vd()
{
}

