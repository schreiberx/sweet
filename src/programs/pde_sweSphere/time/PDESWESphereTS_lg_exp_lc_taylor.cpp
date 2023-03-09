/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */



#include "PDESWESphereTS_lg_exp_lc_taylor.hpp"


void PDESWESphereTS_lg_exp_lc_taylor::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
		run_timestep_lg(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
	}
	else if (timestepping_order == 2)
	{
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
		run_timestep_lg(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
	}
	else
	{
		SWEETError("This time stepping order is not yet supported!");
	}
}



void PDESWESphereTS_lg_exp_lc_taylor::run_timestep_lg(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	timestepping_lg_exp.runTimestep(
			io_phi_pert,
			io_vrt,
			io_div,
			i_fixed_dt,
			i_simulation_timestamp
		);
}



void PDESWESphereTS_lg_exp_lc_taylor::run_timestep_lc(
		sweet::SphereData_Spectral &io_phi,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	/*
	 * Strang split version
	 *
	 * Treat $L_g$ terms exponentially
	 * Treat $L_c$ terms with Taylor expansion
	 */
	sweet::SphereData_Spectral phi_Lc_pow_i_prev = io_phi;
	sweet::SphereData_Spectral vrt_Lc_pow_i_prev = io_vrt;
	sweet::SphereData_Spectral div_Lc_pow_i_prev = io_div;

	double taylor_coeff = 1.0;

	for (int i = 1; i <= shackExpIntegration->taylor_num_expansions; i++)
	{
		// Dummy output variables
		sweet::SphereData_Spectral phi_Lc_pow_i(io_phi.sphereDataConfig, 0);
		sweet::SphereData_Spectral vrt_Lc_pow_i(io_phi.sphereDataConfig, 0);
		sweet::SphereData_Spectral div_Lc_pow_i(io_phi.sphereDataConfig, 0);

		// Time tendencies
		timestepping_lc_erk.euler_timestep_update_lc_spectral_only(
				phi_Lc_pow_i_prev, vrt_Lc_pow_i_prev, div_Lc_pow_i_prev,
				phi_Lc_pow_i, vrt_Lc_pow_i, div_Lc_pow_i,
				i_simulation_timestamp
		);

		// Taylor coefficients
		taylor_coeff *= i_dt/(double)i;

		io_phi += taylor_coeff*phi_Lc_pow_i;
		io_vrt += taylor_coeff*vrt_Lc_pow_i;
		io_div += taylor_coeff*div_Lc_pow_i;

		phi_Lc_pow_i_prev.swap(phi_Lc_pow_i);
		vrt_Lc_pow_i_prev.swap(vrt_Lc_pow_i);
		div_Lc_pow_i_prev.swap(div_Lc_pow_i);
	}
}



bool PDESWESphereTS_lg_exp_lc_taylor::setup(
		sweet::SphereOperators *io_ops,
		int i_order
)
{
	if (shackExpIntegration->taylor_num_expansions <= 0)
		SWEETError("This time stepping method requires setting the number of Taylor expansions");

	timestepping_order = i_order;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETError("Only 1st and 2nd order time stepping order supported");

	timestepping_lg_exp.setup(io_ops, "phi0");

	if (shackExpIntegration->exp_method != "ss_taylor")
		SWEETError("Use --exp-method=ss_taylor to use this time stepper [TODO: This might not make sense after some internal changes]");

	return true;
}



bool PDESWESphereTS_lg_exp_lc_taylor::setup_auto(
		sweet::SphereOperators *io_ops
)
{
	return setup(io_ops, timestepping_order);
}



PDESWESphereTS_lg_exp_lc_taylor::PDESWESphereTS_lg_exp_lc_taylor()
{
}



PDESWESphereTS_lg_exp_lc_taylor::~PDESWESphereTS_lg_exp_lc_taylor()
{
}

