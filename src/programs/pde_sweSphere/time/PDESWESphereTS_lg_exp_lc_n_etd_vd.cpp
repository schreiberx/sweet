/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_exp_lc_n_etd_vd.hpp"


bool PDESWESphereTS_lg_exp_lc_n_etd_vd::implementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (	i_timestepping_method == "lg_exp_lc_n_etd_vd"	||
			i_timestepping_method == "lg_exp_lc_na_nr_etd_vd"	||
			i_timestepping_method == "lg_exp_lc_na_etd_vd"	||
			i_timestepping_method == "lg_exp_lc_nr_etd_vd"	||
			i_timestepping_method == "lg_exp_lc_etd_vd"	||
			false
	)
		return true;

	return false;
}



bool PDESWESphereTS_lg_exp_lc_n_etd_vd::setup_auto(
	const std::string &i_timestepping_method,
	sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere->sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	bool _with_na = false;
	bool _with_nr = false;

	if (	timestepping_method == "lg_exp_lc_n_etd_vd"	||
			timestepping_method == "lg_exp_lc_na_nr_etd_vd")
	{
		_with_na = true;
		_with_nr = true;
	}
	else if (timestepping_method == "lg_exp_lc_na_etd_vd")
	{
		_with_na = true;
		_with_nr = false;
	}
	else if (timestepping_method == "lg_exp_lc_nr_etd_vd")
	{
		_with_na = false;
		_with_nr = true;
	}
	else if (timestepping_method == "lg_exp_lc_etd_vd")
	{
		_with_na = false;
		_with_nr = false;
	}
	else
	{
		SWEETError("Unknown TM");
	}

	return setup_main(
			io_ops,
			shackExpIntegration,
			shackPDESWETimeDisc->timestepping_order,
			shackPDESWETimeDisc->timestepping_order2,
			shackTimestepControl->current_timestepSize,
			_with_na,
			_with_nr
		);
}



bool PDESWESphereTS_lg_exp_lc_n_etd_vd::setup_main(
		const sweet::SphereOperators *io_ops,
		sweet::ShackExpIntegration *i_shackExpIntegration,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestepSize,

		bool i_with_na,
		bool i_with_nr
)
{
	ops = io_ops;
	timestepping_order = i_timestepping_order;

	with_na = i_with_na;
	with_nr = i_with_nr;

	ts_ln_erk_split_vd.setup_main(ops, i_timestepping_order, true, true, true, true, false);

	if (timestepping_order != i_timestepping_order2)
		SWEETError("Mismatch of orders, should be equal");

	if (timestepping_order == 0 || timestepping_order == 1)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestepSize, false, true, timestepping_order);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestepSize, false, true, timestepping_order);
		ts_phi2_exp.setup_variant_50(ops, i_shackExpIntegration, "phi2", i_timestepSize, false, true, timestepping_order);

		NU_phi_prev.setup(ops->sphereDataConfig);
		NU_vrt_prev.setup(ops->sphereDataConfig);
		NU_div_prev.setup(ops->sphereDataConfig);
	}
	else if (timestepping_order == 3)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestepSize, false, true, timestepping_order);
		ts_phi2_exp.setup_variant_50(ops, i_shackExpIntegration, "phi2", i_timestepSize, false, true, timestepping_order);
		ts_phi3_exp.setup_variant_50(ops, i_shackExpIntegration, "phi3", i_timestepSize, false, true, timestepping_order);

		NU_phi_prev.setup(ops->sphereDataConfig);
		NU_vrt_prev.setup(ops->sphereDataConfig);
		NU_div_prev.setup(ops->sphereDataConfig);

		NU_phi_prev_2.setup(ops->sphereDataConfig);
		NU_vrt_prev_2.setup(ops->sphereDataConfig);
		NU_div_prev_2.setup(ops->sphereDataConfig);
	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}

	return true;
}



std::string PDESWESphereTS_lg_exp_lc_n_etd_vd::getIDString()
{
	return "lg_exp_lc_n_etd";
}



void PDESWESphereTS_lg_exp_lc_n_etd_vd::printHelp()
{
	std::cout << "	Exponential ETD:" << std::endl;
	std::cout << "		+ lg_exp_lc_n_etd_vd" << std::endl;
	std::cout << "		+ lg_exp_lc_na_nr_etd_vd" << std::endl;
	std::cout << "		+ lg_exp_lc_na_etd_vd" << std::endl;
	std::cout << "		+ lg_exp_lc_nr_etd_vd" << std::endl;
	std::cout << "		+ lg_exp_lc_etd_vd" << std::endl;
}




void PDESWESphereTS_lg_exp_lc_n_etd_vd::runTimestep(
		sweet::SphereData_Spectral &io_U_phi,
		sweet::SphereData_Spectral &io_U_vrt,
		sweet::SphereData_Spectral &io_U_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	const sweet::SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		sweet::SphereData_Spectral phi0_Un_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_div(sphereDataConfig);

		ts_phi0_exp.runTimestep(
				io_U_phi, io_U_vrt, io_U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral F_phi(sphereDataConfig, 0);
		sweet::SphereData_Spectral F_vrt(sphereDataConfig, 0);
		sweet::SphereData_Spectral F_div(sphereDataConfig, 0);

		ts_ln_erk_split_vd.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				F_phi, F_vrt, F_div,
				i_simulation_timestamp
		);


		if (with_na)
		{
			ts_ln_erk_split_vd.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_vd.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}


		sweet::SphereData_Spectral phi1_FUn_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_div(sphereDataConfig);

		ts_phi1_exp.runTimestep(
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
		const sweet::SphereData_Spectral &U_phi = io_U_phi;
		const sweet::SphereData_Spectral &U_vrt = io_U_vrt;
		const sweet::SphereData_Spectral &U_div = io_U_div;

		if (i_simulation_timestamp == 0)
		{
			/*
			 * First time step:
			 * Simply backup existing fields for multi-step parts of this algorithm.
			 */
			NU_phi_prev.spectral_set_zero();
			NU_vrt_prev.spectral_set_zero();
			NU_div_prev.spectral_set_zero();

			ts_ln_erk_split_vd.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					NU_phi_prev, NU_vrt_prev, NU_div_prev,
					i_simulation_timestamp
			);

			if (with_na)
			{
				ts_ln_erk_split_vd.euler_timestep_update_na(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			if (with_nr)
			{
				ts_ln_erk_split_vd.euler_timestep_update_nr(
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
		sweet::SphereData_Spectral phi0_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi0_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi0_div(sphereDataConfig);

		ts_phi0_exp.runTimestep(
				U_phi, U_vrt, U_div,
				phi0_phi, phi0_vrt, phi0_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi1
		 */
		sweet::SphereData_Spectral F_phi(sphereDataConfig, 0);
		sweet::SphereData_Spectral F_vrt(sphereDataConfig, 0);
		sweet::SphereData_Spectral F_div(sphereDataConfig, 0);

		ts_ln_erk_split_vd.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				F_phi, F_vrt, F_div,
				i_simulation_timestamp
		);

		if (with_na)
		{
			ts_ln_erk_split_vd.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_vd.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		sweet::SphereData_Spectral phi1_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi1_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi1_div(sphereDataConfig);

		ts_phi1_exp.runTimestep(
				F_phi, F_vrt, F_div,
				phi1_phi, phi1_vrt, phi1_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi2
		 */

		sweet::SphereData_Spectral phi2_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi2_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi2_div(sphereDataConfig);

		ts_phi2_exp.runTimestep(
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
		const sweet::SphereData_Spectral &U_phi = io_U_phi;
		const sweet::SphereData_Spectral &U_vrt = io_U_vrt;
		const sweet::SphereData_Spectral &U_div = io_U_div;

		if (i_simulation_timestamp == 0)
		{
			/*
			 * First time step:
			 * Simply backup existing fields for multi-step parts of this algorithm.
			 */
			NU_phi_prev.spectral_set_zero();
			NU_vrt_prev.spectral_set_zero();
			NU_div_prev.spectral_set_zero();

			ts_ln_erk_split_vd.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					NU_phi_prev, NU_vrt_prev, NU_div_prev,
					i_simulation_timestamp
			);

			if (with_na)
			{
				ts_ln_erk_split_vd.euler_timestep_update_na(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			if (with_nr)
			{
				ts_ln_erk_split_vd.euler_timestep_update_nr(
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
		sweet::SphereData_Spectral NU_phi(sphereDataConfig, 0);
		sweet::SphereData_Spectral NU_vrt(sphereDataConfig, 0);
		sweet::SphereData_Spectral NU_div(sphereDataConfig, 0);

		ts_ln_erk_split_vd.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				NU_phi, NU_vrt, NU_div,
				i_simulation_timestamp
		);

		if (with_na)
		{
			ts_ln_erk_split_vd.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					NU_phi, NU_vrt, NU_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_vd.euler_timestep_update_nr(
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
		sweet::SphereData_Spectral phi0_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi0_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi0_div(sphereDataConfig);

		ts_phi0_exp.runTimestep(
				U_phi, U_vrt, U_div,
				phi0_phi, phi0_vrt, phi0_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi1
		 */
		sweet::SphereData_Spectral phi1_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi1_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi1_div(sphereDataConfig);

		ts_phi1_exp.runTimestep(
				NU_phi, NU_vrt, NU_div,
				phi1_phi, phi1_vrt, phi1_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi2
		 */
		sweet::SphereData_Spectral phi2_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi2_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi2_div(sphereDataConfig);

		ts_phi2_exp.runTimestep(
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
		sweet::SphereData_Spectral phi3_phi(sphereDataConfig);
		sweet::SphereData_Spectral phi3_vrt(sphereDataConfig);
		sweet::SphereData_Spectral phi3_div(sphereDataConfig);

		ts_phi3_exp.runTimestep(
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




PDESWESphereTS_lg_exp_lc_n_etd_vd::PDESWESphereTS_lg_exp_lc_n_etd_vd()
{
}



PDESWESphereTS_lg_exp_lc_n_etd_vd::~PDESWESphereTS_lg_exp_lc_n_etd_vd()
{
}

