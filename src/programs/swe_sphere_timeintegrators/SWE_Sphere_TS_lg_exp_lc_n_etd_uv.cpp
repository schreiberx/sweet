/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_exp_lc_n_etd_uv.hpp"


bool SWE_Sphere_TS_lg_exp_lc_n_etd_uv::implements_timestepping_method(const std::string &i_timestepping_method
									)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;
	if (	i_timestepping_method == "lg_exp_lc_n_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_na_nr_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_na_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_nr_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_etd_uv"	||
			false
	)
		return true;

	return false;
}


std::string SWE_Sphere_TS_lg_exp_lc_n_etd_uv::string_id()
{
	return "lg_exp_lc_n_etd";
}


void SWE_Sphere_TS_lg_exp_lc_n_etd_uv::setup_auto()
{
	if (simVars.sim.sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	if (	timestepping_method == "lg_exp_lc_n_etd_uv"	||
			timestepping_method == "lg_exp_lc_na_nr_etd_uv")
	{
		with_na = true;
		with_nr = true;
	}
	else if (timestepping_method == "lg_exp_lc_na_etd_uv")
	{
		with_na = true;
		with_nr = false;
	}
	else if (timestepping_method == "lg_exp_lc_nr_etd_uv")
	{
		with_na = false;
		with_nr = true;
	}
	else if (timestepping_method == "lg_exp_lc_etd_uv")
	{
		with_na = false;
		with_nr = false;
	}
	else
	{
		SWEETError("Unknown TM");
	}

	setup(
			simVars.rexi,
			timestepping_order,
			timestepping_order2,
			simVars.timecontrol.current_timestep_size
		);
}



void SWE_Sphere_TS_lg_exp_lc_n_etd_uv::print_help()
{
	std::cout << "	Exponential ETD:" << std::endl;
	std::cout << "		+ lg_exp_lc_n_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_na_nr_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_na_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_nr_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_etd_uv" << std::endl;
}




void SWE_Sphere_TS_lg_exp_lc_n_etd_uv::run_timestep(
		SphereData_Spectral &io_U_phi,	///< prognostic variables
		SphereData_Spectral &io_U_vrt,	///< prognostic variables
		SphereData_Spectral &io_U_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		SphereData_Spectral phi0_Un_phi(sphereDataConfig);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig);
		SphereData_Spectral phi0_Un_div(sphereDataConfig);

		ts_phi0_exp.run_timestep(
				io_U_phi, io_U_vrt, io_U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		SphereData_Spectral F_phi(sphereDataConfig, 0);
		SphereData_Spectral F_vrt(sphereDataConfig, 0);
		SphereData_Spectral F_div(sphereDataConfig, 0);

		ts_ln_erk_split_uv.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				F_phi, F_vrt, F_div,
				i_simulation_timestamp
		);


		if (with_na)
		{
			ts_ln_erk_split_uv.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}


		SphereData_Spectral phi1_FUn_phi(sphereDataConfig);
		SphereData_Spectral phi1_FUn_vrt(sphereDataConfig);
		SphereData_Spectral phi1_FUn_div(sphereDataConfig);

		ts_phi1_exp.run_timestep(
				F_phi, F_vrt, F_div,
				phi1_FUn_phi, phi1_FUn_vrt, phi1_FUn_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		io_U_phi = phi0_Un_phi + i_fixed_dt*phi1_FUn_phi;
		io_U_vrt = phi0_Un_vrt + i_fixed_dt*phi1_FUn_vrt;
		io_U_div = phi0_Un_div + i_fixed_dt*phi1_FUn_div;
	}
	else if (timestepping_order == 2)
	{
		const SphereData_Spectral &U_phi = io_U_phi;
		const SphereData_Spectral &U_vrt = io_U_vrt;
		const SphereData_Spectral &U_div = io_U_div;

		if (i_simulation_timestamp == 0)
		{
			/*
			 * First time step:
			 * Simply backup existing fields for multi-step parts of this algorithm.
			 */
			NU_phi_prev.spectral_set_zero();
			NU_vrt_prev.spectral_set_zero();
			NU_div_prev.spectral_set_zero();

			ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					NU_phi_prev, NU_vrt_prev, NU_div_prev,
					i_simulation_timestamp
			);

			if (with_na)
			{
				ts_ln_erk_split_uv.euler_timestep_update_na(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			if (with_nr)
			{
				ts_ln_erk_split_uv.euler_timestep_update_nr(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

		}


		/*
		 * u_{n+1} = \varphi_{0} (\Delta t L )u_{n} + \Delta t K
		 *
		 * K = \left[
		 * 		\varphi_{1}(\Delta t L)F(u_{n}) +
		 * 		\varphi_{2}(\Delta t L) (F(u_{n})-F(u_{n-1}))
		 * 	\right]
		 */

		/*
		 * phi0
		 */
		SphereData_Spectral phi0_phi(sphereDataConfig);
		SphereData_Spectral phi0_vrt(sphereDataConfig);
		SphereData_Spectral phi0_div(sphereDataConfig);

		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_phi, phi0_vrt, phi0_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi1
		 */
		SphereData_Spectral F_phi(sphereDataConfig, 0);
		SphereData_Spectral F_vrt(sphereDataConfig, 0);
		SphereData_Spectral F_div(sphereDataConfig, 0);

		ts_ln_erk_split_uv.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				F_phi, F_vrt, F_div,
				i_simulation_timestamp
		);

		if (with_na)
		{
			ts_ln_erk_split_uv.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		SphereData_Spectral phi1_phi(sphereDataConfig);
		SphereData_Spectral phi1_vrt(sphereDataConfig);
		SphereData_Spectral phi1_div(sphereDataConfig);

		ts_phi1_exp.run_timestep(
				F_phi, F_vrt, F_div,
				phi1_phi, phi1_vrt, phi1_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi2
		 */

		SphereData_Spectral phi2_phi(sphereDataConfig);
		SphereData_Spectral phi2_vrt(sphereDataConfig);
		SphereData_Spectral phi2_div(sphereDataConfig);

		ts_phi2_exp.run_timestep(
				F_phi - NU_phi_prev,
				F_vrt - NU_vrt_prev,
				F_div - NU_div_prev,

				phi2_phi,
				phi2_vrt,
				phi2_div,

				i_fixed_dt,
				i_simulation_timestamp
			);

		io_U_phi = phi0_phi + i_fixed_dt*(phi1_phi + phi2_phi);
		io_U_vrt = phi0_vrt + i_fixed_dt*(phi1_vrt + phi2_vrt);
		io_U_div = phi0_div + i_fixed_dt*(phi1_div + phi2_div);

		/*
		 * Backup nonlinear evaluations
		 */
		{
			NU_phi_prev = F_phi;
			NU_vrt_prev = F_vrt;
			NU_div_prev = F_div;
		}
	}
	else if (timestepping_order == 3)
	{
		const SphereData_Spectral &U_phi = io_U_phi;
		const SphereData_Spectral &U_vrt = io_U_vrt;
		const SphereData_Spectral &U_div = io_U_div;

		if (i_simulation_timestamp == 0)
		{
			/*
			 * First time step:
			 * Simply backup existing fields for multi-step parts of this algorithm.
			 */
			NU_phi_prev.spectral_set_zero();
			NU_vrt_prev.spectral_set_zero();
			NU_div_prev.spectral_set_zero();

			ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					NU_phi_prev, NU_vrt_prev, NU_div_prev,
					i_simulation_timestamp
			);

			if (with_na)
			{
				ts_ln_erk_split_uv.euler_timestep_update_na(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			if (with_nr)
			{
				ts_ln_erk_split_uv.euler_timestep_update_nr(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			NU_phi_prev_2 = NU_phi_prev;
			NU_vrt_prev_2 = NU_vrt_prev;
			NU_div_prev_2 = NU_div_prev;
		}


		/*
		 * Compute nonlinear tendencies
		 */
		SphereData_Spectral NU_phi(sphereDataConfig, 0);
		SphereData_Spectral NU_vrt(sphereDataConfig, 0);
		SphereData_Spectral NU_div(sphereDataConfig, 0);

		ts_ln_erk_split_uv.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				NU_phi, NU_vrt, NU_div,
				i_simulation_timestamp
		);

		if (with_na)
		{
			ts_ln_erk_split_uv.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					NU_phi, NU_vrt, NU_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					NU_phi, NU_vrt, NU_div,
					i_simulation_timestamp
			);
		}


		/*
		 * u_{n+1} = \varphi_{0} (\Delta t L )u_{n} + \Delta t K
		 *
		 * K = \left[
		 * 		\varphi_{1}(\Delta t L)F(u_{n}) +
		 * 		\varphi_{2}(\Delta t L) (F(u_{n})-F(u_{n-1}))
		 * 	\right]
		 */

		/*
		 * phi0
		 */
		SphereData_Spectral phi0_phi(sphereDataConfig);
		SphereData_Spectral phi0_vrt(sphereDataConfig);
		SphereData_Spectral phi0_div(sphereDataConfig);

		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_phi, phi0_vrt, phi0_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi1
		 */
		SphereData_Spectral phi1_phi(sphereDataConfig);
		SphereData_Spectral phi1_vrt(sphereDataConfig);
		SphereData_Spectral phi1_div(sphereDataConfig);

		ts_phi1_exp.run_timestep(
				NU_phi, NU_vrt, NU_div,
				phi1_phi, phi1_vrt, phi1_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi2
		 */
		SphereData_Spectral phi2_phi(sphereDataConfig);
		SphereData_Spectral phi2_vrt(sphereDataConfig);
		SphereData_Spectral phi2_div(sphereDataConfig);

		ts_phi2_exp.run_timestep(
				3.0/2.0*NU_phi - 2.0*NU_phi_prev + 1.0/2.0*NU_phi_prev_2,
				3.0/2.0*NU_vrt - 2.0*NU_vrt_prev + 1.0/2.0*NU_vrt_prev_2,
				3.0/2.0*NU_div - 2.0*NU_div_prev + 1.0/2.0*NU_div_prev_2,

				phi2_phi,
				phi2_vrt,
				phi2_div,

				i_fixed_dt,
				i_simulation_timestamp
			);


		/*
		 * phi3
		 */
		SphereData_Spectral phi3_phi(sphereDataConfig);
		SphereData_Spectral phi3_vrt(sphereDataConfig);
		SphereData_Spectral phi3_div(sphereDataConfig);

		ts_phi3_exp.run_timestep(
				1.0/2.0*NU_phi - NU_phi_prev + 1.0/2.0*NU_phi_prev_2,
				1.0/2.0*NU_vrt - NU_vrt_prev + 1.0/2.0*NU_vrt_prev_2,
				1.0/2.0*NU_div - NU_div_prev + 1.0/2.0*NU_div_prev_2,

				phi3_phi,
				phi3_vrt,
				phi3_div,

				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * Compute full time step
		 */
		io_U_phi = phi0_phi + i_fixed_dt*(phi1_phi + phi2_phi + phi3_phi);
		io_U_vrt = phi0_vrt + i_fixed_dt*(phi1_vrt + phi2_vrt + phi3_vrt);
		io_U_div = phi0_div + i_fixed_dt*(phi1_div + phi2_div + phi3_div);

		/*
		 * Backup nonlinear evaluations
		 */
		{
			NU_phi_prev_2 = NU_phi_prev;
			NU_vrt_prev_2 = NU_vrt_prev;
			NU_div_prev_2 = NU_div_prev;

			NU_phi_prev = NU_phi;
			NU_vrt_prev = NU_vrt;
			NU_div_prev = NU_div;
		}
	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_lg_exp_lc_n_etd_uv::setup(
		EXP_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size
)
{
	timestepping_order = i_timestepping_order;

	ts_ln_erk_split_uv.setup(i_timestepping_order, true, true, true, true, false);

	if (timestepping_order != i_timestepping_order2)
		SWEETError("Mismatch of orders, should be equal");

	if (timestepping_order == 0 || timestepping_order == 1)
	{
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup(i_rexiSimVars, "phi1", i_timestep_size, false, true, timestepping_order);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup(i_rexiSimVars, "phi1", i_timestep_size, false, true, timestepping_order);
		ts_phi2_exp.setup(i_rexiSimVars, "phi2", i_timestep_size, false, true, timestepping_order);

		NU_phi_prev.setup(ops.sphereDataConfig);
		NU_vrt_prev.setup(ops.sphereDataConfig);
		NU_div_prev.setup(ops.sphereDataConfig);
	}
	else if (timestepping_order == 3)
	{
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup(i_rexiSimVars, "phi1", i_timestep_size, false, true, timestepping_order);
		ts_phi2_exp.setup(i_rexiSimVars, "phi2", i_timestep_size, false, true, timestepping_order);
		ts_phi3_exp.setup(i_rexiSimVars, "phi3", i_timestep_size, false, true, timestepping_order);

		NU_phi_prev.setup(ops.sphereDataConfig);
		NU_vrt_prev.setup(ops.sphereDataConfig);
		NU_div_prev.setup(ops.sphereDataConfig);

		NU_phi_prev_2.setup(ops.sphereDataConfig);
		NU_vrt_prev_2.setup(ops.sphereDataConfig);
		NU_div_prev_2.setup(ops.sphereDataConfig);
	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}
}


SWE_Sphere_TS_lg_exp_lc_n_etd_uv::SWE_Sphere_TS_lg_exp_lc_n_etd_uv(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		ops(i_op),

		ts_ln_erk_split_uv(simVars, ops),

		with_na(true),
		with_nr(true),

		ts_phi0_exp(simVars, ops),
		ts_phi1_exp(simVars, ops),
		ts_phi2_exp(simVars, ops),
		ts_phi3_exp(simVars, ops)
{
}



SWE_Sphere_TS_lg_exp_lc_n_etd_uv::~SWE_Sphere_TS_lg_exp_lc_n_etd_uv()
{
}

