/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */



#include "SWE_Sphere_TS_lg_exp_lc_taylor.hpp"


void SWE_Sphere_TS_lg_exp_lc_taylor::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

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



void SWE_Sphere_TS_lg_exp_lc_taylor::run_timestep_lg(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	timestepping_lg_exp.run_timestep(
			io_phi_pert,
			io_vrt,
			io_div,
			i_fixed_dt,
			i_simulation_timestamp
		);
}



void SWE_Sphere_TS_lg_exp_lc_taylor::run_timestep_lc(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

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
	SphereData_Spectral phi_Lc_pow_i_prev = io_phi;
	SphereData_Spectral vrt_Lc_pow_i_prev = io_vrt;
	SphereData_Spectral div_Lc_pow_i_prev = io_div;

	double taylor_coeff = 1.0;

	for (int i = 1; i <= simVars.rexi.taylor_num_expansions; i++)
	{
		// Dummy output variables
		SphereData_Spectral phi_Lc_pow_i(io_phi.sphereDataConfig, 0);
		SphereData_Spectral vrt_Lc_pow_i(io_phi.sphereDataConfig, 0);
		SphereData_Spectral div_Lc_pow_i(io_phi.sphereDataConfig, 0);

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



void SWE_Sphere_TS_lg_exp_lc_taylor::setup(
		int i_order
)
{
	if (simVars.rexi.taylor_num_expansions <= 0)
		SWEETError("This time stepping method requires setting the number of Taylor expansions");

	timestepping_order = i_order;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETError("Only 1st and 2nd order time stepping order supported");

	timestepping_lg_exp.setup("phi0");

	if (simVars.rexi.exp_method != "ss_taylor")
		SWEETError("Use --exp-method=ss_taylor to use this time stepper");
}



void SWE_Sphere_TS_lg_exp_lc_taylor::setup_auto()
{
	setup(timestepping_order);
}



SWE_Sphere_TS_lg_exp_lc_taylor::SWE_Sphere_TS_lg_exp_lc_taylor(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_lg_exp(i_simVars, i_op),
		timestepping_lc_erk(i_simVars, i_op)
{
}



SWE_Sphere_TS_lg_exp_lc_taylor::~SWE_Sphere_TS_lg_exp_lc_taylor()
{
}

