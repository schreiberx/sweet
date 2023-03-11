/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_exp_lc_n_erk.hpp"



void PDESWESphereTS_lg_exp_lc_n_erk::runTimestep(
		sweet::SphereData_Spectral &io_phi,		///< prognostic variables
		sweet::SphereData_Spectral &io_vort,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		sweet::SphereData_Spectral tmp_phi = io_phi;
		sweet::SphereData_Spectral tmp_vort = io_vort;
		sweet::SphereData_Spectral tmp_div = io_div;

		// first order IRK for linear
		timestepping_lg_rexi.runTimestep(
				io_phi, io_vort, io_div,
				i_dt,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral phi_dt(io_phi.sphereDataConfig);
		sweet::SphereData_Spectral vort_dt(io_vort.sphereDataConfig);
		sweet::SphereData_Spectral div_dt(io_div.sphereDataConfig);

		// first order explicit for non-linear
		timestepping_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				tmp_phi, tmp_vort, tmp_div,
				phi_dt, vort_dt, div_dt,
				i_simulation_timestamp
			);


		io_phi += i_dt*phi_dt;
		io_vort += i_dt*vort_dt;
		io_div += i_dt*div_dt;
	}
	else if (timestepping_order == 2 || timestepping_order == 4)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_lg_rexi.runTimestep(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_n_erk,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute Euler time step updates
					io_phi, io_vort, io_div,
					i_dt,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_rexi.runTimestep(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_n_erk,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_rexi.runTimestep(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_n_erk,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
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


bool PDESWESphereTS_lg_exp_lc_n_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	int version = 0;
	if (timestepping_method == "lg_exp_lc_n_erk_ver1")
		version = 1;

	return setup_main(
			io_ops,
			shackExpIntegration,
			shackPDESWETimeDisc->timestepping_order,
			shackPDESWETimeDisc->timestepping_order2,
			shackTimestepControl->current_timestep_size,
			version
		);
}


/*
 * Setup
 */
bool PDESWESphereTS_lg_exp_lc_n_erk::setup_main(
		sweet::SphereOperators *io_ops,
		sweet::ShackExpIntegration *i_shackExpIntegration,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size,
		int i_version_id
)
{
	ops = io_ops;
	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	timestep_size = shackTimestepControl->current_timestep_size;

	if (timestepping_order2 == 1)
	{
		timestepping_lg_rexi.setup_variant_10(
				ops,
				i_shackExpIntegration,
				"phi0",
				i_timestep_size,
				false,
				true
			);
	}
	else if (timestepping_order2 == 2 || timestepping_order2 == 4)
	{
		if (version_id == 0)
		{
			timestepping_lg_rexi.setup_variant_10(
					ops,
					i_shackExpIntegration,
					"phi0",
					i_timestep_size*0.5,
					false,
					true
			);
		}
		else if (version_id == 1)
		{
			timestepping_lg_rexi.setup_variant_10(
					ops,
					i_shackExpIntegration,
					"phi0",
					i_timestep_size,
					false,
					true
			);
		}
		else
		{
			SWEETError("Invalid version id");
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
	timestepping_lg_erk_lc_n_erk.setup(ops, 1, -1);

	return true;
}



PDESWESphereTS_lg_exp_lc_n_erk::PDESWESphereTS_lg_exp_lc_n_erk()
{
}



PDESWESphereTS_lg_exp_lc_n_erk::~PDESWESphereTS_lg_exp_lc_n_erk()
{
}

