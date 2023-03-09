/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_l_exp_n_erk.hpp"


bool PDESWESphereTS_l_exp_n_erk::setup(
		sweet::SphereOperators *io_ops,
		sweet::ShackExpIntegration *i_shackExpIntegration,
		const std::string &i_exp_method,
		int i_order,	///< order of RK time stepping method
		int i_order2,	///< order of RK time stepping method of non-linear parts
		double i_timestep_size,
		bool i_use_f_sphere,
		int i_version_id,
		bool i_use_rexi_sphere_solver_preallocation
)
{
	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	timestep_size = i_timestep_size;

	version_id = i_version_id;

	if (timestepping_order != timestepping_order2)
		SWEETError("Mismatch of timestepping method orders, they should be equal");


	if (timestepping_order2 == 1)
	{
		timestepping_l_rexi.setup_variant_100(
				io_ops,
				i_shackExpIntegration,
				"phi0",
				i_exp_method,
				i_timestep_size,
				i_use_f_sphere,
				false,
				timestepping_order,
				i_use_rexi_sphere_solver_preallocation
		);
	}
	else if (timestepping_order2 == 2 || timestepping_order2 == 4)
	{
		if (version_id == 0)
		{
			timestepping_l_rexi.setup_variant_100(
					io_ops,
					i_shackExpIntegration,
					"phi0",
					i_exp_method,
					i_timestep_size*0.5,
					i_use_f_sphere,
					false,
					timestepping_order,
					i_use_rexi_sphere_solver_preallocation
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_rexi.setup_variant_100(
					io_ops,
					i_shackExpIntegration,
					"phi0",
					i_exp_method,
					i_timestep_size,
					i_use_f_sphere,
					false,
					timestepping_order,
					i_use_rexi_sphere_solver_preallocation
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
	timestepping_l_erk_n_erk.setup(io_ops, 1, 1);
	return true;
}



bool PDESWESphereTS_l_exp_n_erk::setup_auto(
		sweet::SphereOperators *io_ops
)
{
	int _version_id = 0;

	if (timestepping_method == "l_exp_n_erk_ver1")
		_version_id = 1;

	return setup(
			io_ops,
			shackExpIntegration,
			shackExpIntegration->exp_method,
			timestepping_order,
			timestepping_order2,
			shackTimestepControl->current_timestep_size,
			shackPDESWESphere->sphere_use_fsphere,
			_version_id,
			shackExpIntegration->sphere_solver_preallocation
		);
}




void PDESWESphereTS_l_exp_n_erk::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order2 == 1)
	{
		if (version_id == 0)
		{
			// first order REXI for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else if (version_id == 1)
		{
			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order REXI for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETError("Invalid version id");
		}
	}
	else if (timestepping_order2 == 2 || timestepping_order2)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphereTS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);

		}
		else if (version_id == 1)
		{

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphereTS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphereTS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
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



PDESWESphereTS_l_exp_n_erk::PDESWESphereTS_l_exp_n_erk()	:
		version_id(0)
{
}



PDESWESphereTS_l_exp_n_erk::~PDESWESphereTS_l_exp_n_erk()
{
}

