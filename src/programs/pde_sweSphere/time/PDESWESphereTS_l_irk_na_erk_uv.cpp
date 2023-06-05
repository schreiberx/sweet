/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_l_irk_na_erk_uv.hpp"



bool PDESWESphereTS_l_irk_na_erk_uv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (
		timestepping_method == "l_irk_na_erk_uv" ||
		timestepping_method == "l_irk_na_erk_uv_ver0"
	)
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, shackPDESWETimeDisc->timestepping_order2, 0);
	else if (
			timestepping_method == "l_irk_na_erk_uv_ver1"
		)
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, shackPDESWETimeDisc->timestepping_order2, 1);
	else
		SWEETError("Not implemented");

	return false;
}


/*
 * Setup
 */
bool PDESWESphereTS_l_irk_na_erk_uv::setup_main(
		const sweet::SphereOperators *io_ops,
		int i_order,	///< order of RK time stepping method
		int i_order2,	///< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	ops = io_ops;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	version_id = i_version_id;

	if (timestepping_order2 < 0)
		timestepping_order2 = timestepping_order;

	if (i_order != i_order2)
		SWEETError("Orders of 1st and 2nd one must match");


	timestep_size = shackTimestepControl->current_timestepSize;

	if (timestepping_order == 1)
	{
		timestepping_l_irk.setup(
				ops,
				1,
				timestep_size
		);

		//
		// Only request 1st order time stepping methods for irk and erk
		// These 1st order methods will be combined to higher-order methods in this class
		//
		timestepping_na_erk_split_uv.setup_main(ops, 1, false, false, true, false, false);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_l_irk.setup_main(
					ops,
					2,
					timestep_size*0.5,
					shackPDESWETimeDisc->timestepping_crank_nicolson_filter,
					false
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_irk.setup_main(
					ops,
					2,
					timestep_size,
					shackPDESWETimeDisc->timestepping_crank_nicolson_filter,
					false
			);
		}
		else
		{
			SWEETError("Invalid version");
		}

		//
		// Only request 1st order time stepping methods for irk and erk
		// These 1st order methods will be combined to higher-order methods in this class
		//
		timestepping_na_erk_split_uv.setup_main(ops, 2, false, false, true, false, false);
	}
	else
	{
		SWEETError("Invalid timestepping order");
	}

	return true;
}



bool PDESWESphereTS_l_irk_na_erk_uv::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (
		i_timestepping_method == "l_irk_na_erk_uv" || i_timestepping_method == "l_irk_na_erk_uv_ver0" ||
		i_timestepping_method == "l_irk_na_erk_uv_ver1"
	)
		return true;

	return false;
}



std::string PDESWESphereTS_l_irk_na_erk_uv::getIDString()
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



void PDESWESphereTS_l_irk_na_erk_uv::runTimestep(
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
			timestepping_l_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_na_erk_split_uv.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_na_erk_split_uv.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_l_irk.runTimestep(
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
			timestepping_l_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_na_erk_split_uv.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_na_erk_split_uv.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_na_erk_split_uv.runTimestep(
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



PDESWESphereTS_l_irk_na_erk_uv::PDESWESphereTS_l_irk_na_erk_uv()	:
		version_id(0)
{

}



PDESWESphereTS_l_irk_na_erk_uv::~PDESWESphereTS_l_irk_na_erk_uv()
{
}

