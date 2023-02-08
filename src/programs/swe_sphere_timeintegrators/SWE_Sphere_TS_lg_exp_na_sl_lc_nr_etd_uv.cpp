/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv.hpp"


bool SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::implements_timestepping_method(const std::string &i_timestepping_method
									)
{
#if 0
	/*
	 * Should contain _exp and _settls as well as _uv to indicate vorticity-divergence formulation
	 */
	return (
		(i_timestepping_method.find("_exp") != std::string::npos)		&&
		(i_timestepping_method.find("_etd") != std::string::npos)	&&
		(i_timestepping_method.find("_uv") != std::string::npos)		&&
		true
	);
#endif

	timestepping_method = i_timestepping_method;
	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;

	if (	i_timestepping_method == "lg_exp_na_sl_lc_nr_etd_uv"	||
			i_timestepping_method == "lg_exp_na_sl_lc_etd_uv"	||
			false
	)
		return true;

	return false;
}


std::string SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::string_id()
{
	return "lg_exp_na_sl_lc_nr_etd_uv";
}


void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::setup_auto()
{
	if (simVars.sim.sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	NLRemainderTreatment_enum nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;

	if (this->timestepping_method == "lg_exp_na_sl_lc_nr_etd_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;
	}
	else if (this->timestepping_method == "lg_exp_na_sl_lc_etd_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_IGNORE;
	}
	else
	{
		SWEETError(std::string("Timestepping method '")+this->timestepping_method+std::string("' not known"));
	}

	setup(
			simVars.rexi,
			timestepping_order,
			timestepping_order2,
			simVars.timecontrol.current_timestep_size,

			nonlinear_remainder_treatment
		);
}



void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::print_help()
{
	std::cout << "	Exponential SL ETD:" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_etd_uv" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_nr_etd_uv" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_lc_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_nr_etd_uv" << std::endl;
#if 0
	std::cout << "		+ l_exp_na_sl_etd_uv" << std::endl;
	std::cout << "		+ l_exp_na_sl_nr_etd_uv" << std::endl;
#endif
}



void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::run_timestep(
		SphereData_Spectral &io_U_phi,	///< prognostic variables
		SphereData_Spectral &io_U_vrt,	///< prognostic variables
		SphereData_Spectral &io_U_div,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;


	if (timestepping_order != 2)
	{
		SWEETError("Only 2nd order supported");
	}


	const SphereData_Spectral &U_phi = io_U_phi;
	const SphereData_Spectral &U_vrt = io_U_vrt;
	const SphereData_Spectral &U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
#if !SWEET_PARAREAL
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 */
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
#endif
	}


	/*************************************************************************************************
	 * Step 1) Compute departure points
	 *************************************************************************************************/
	SphereData_Physical U_u_lon_prev, U_v_lat_prev;
	ops.vrtdiv_to_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev);

	SphereData_Physical U_u_lon, U_v_lat;
	ops.vrtdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat);

	double dt_div_radius = i_dt / simVars.sim.sphere_radius;

	// Calculate departure points
	ScalarDataArray pos_lon_d, pos_lat_d;
	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_lon_prev, dt_div_radius*U_v_lat_prev,
			dt_div_radius*U_u_lon, dt_div_radius*U_v_lat,
			pos_lon_d, pos_lat_d		// OUTPUT
		);


	/*************************************************************************************************
	 * Step 2) Compute U_D at departure points
	 *************************************************************************************************/

	SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
	semiLagrangian.apply_sl_timeintegration_uv(
			ops,
			U_phi, U_vrt, U_div,
			pos_lon_d, pos_lat_d,
			U_phi_D, U_vrt_D, U_div_D
		);



	/*
	 * u^{n+1} = \varphi_{0} (\Delta t L )u^{n}_D + \Delta t (K + K_D)
	 *
	 * K = \left[
	 * 		\varphi_{1}(\Delta t L)F(u^{n}) +
	 * 		\varphi_{2}(\Delta t L) (F(u^{n})-F(u^{n-1}))
	 * 	\right]
	 */


	/*
	 * phi0
	 */
	SphereData_Spectral phi0_phi(sphereDataConfig);
	SphereData_Spectral phi0_vrt(sphereDataConfig);
	SphereData_Spectral phi0_div(sphereDataConfig);

	ts_phi0_exp.run_timestep(
			//U_phi, U_vrt, U_div,
			U_phi_D, U_vrt_D, U_div_D,
			phi0_phi, phi0_vrt, phi0_div,
			i_dt,
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

	if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
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
			i_dt,
			i_simulation_timestamp
		);


	/*
	 * phi2
	 */
	SphereData_Spectral F_phi_prev(sphereDataConfig, 0);
	SphereData_Spectral F_vrt_prev(sphereDataConfig, 0);
	SphereData_Spectral F_div_prev(sphereDataConfig, 0);

	ts_ln_erk_split_uv.euler_timestep_update_lc(
			U_phi_prev, U_vrt_prev, U_div_prev,
			F_phi_prev, F_vrt_prev, F_div_prev,
			i_simulation_timestamp
	);



	if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
	{
		ts_ln_erk_split_uv.euler_timestep_update_nr(
				U_phi_prev, U_vrt_prev, U_div_prev,
				F_phi_prev, F_vrt_prev, F_div_prev,
				i_simulation_timestamp
		);
	}

	SphereData_Spectral phi2_phi(sphereDataConfig);
	SphereData_Spectral phi2_vrt(sphereDataConfig);
	SphereData_Spectral phi2_div(sphereDataConfig);

	ts_phi2_exp.run_timestep(
			F_phi - F_phi_prev,
			F_vrt - F_vrt_prev,
			F_div - F_div_prev,

			phi2_phi,
			phi2_vrt,
			phi2_div,

			i_dt,
			i_simulation_timestamp
		);

	SphereData_Spectral K_phi = phi1_phi + phi2_phi;
	SphereData_Spectral K_vrt = phi1_vrt + phi2_vrt;
	SphereData_Spectral K_div = phi1_div + phi2_div;


	/*
	 * Compute K at departure points
	 */

	SphereData_Spectral K_phi_D, K_vrt_D, K_div_D;
	semiLagrangian.apply_sl_timeintegration_uv(
			ops,
			K_phi, K_vrt, K_div,
			pos_lon_d, pos_lat_d,
			K_phi_D, K_vrt_D, K_div_D
		);

	io_U_phi = phi0_phi + (i_dt*0.5)*(K_phi + K_phi_D);
	io_U_vrt = phi0_vrt + (i_dt*0.5)*(K_vrt + K_vrt_D);
	io_U_div = phi0_div + (i_dt*0.5)*(K_div + K_div_D);


	{
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::setup(
		EXP_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size,

		NLRemainderTreatment_enum i_nonlinear_remainder_treatment
)
{
	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	nonlinear_remainder_treatment = i_nonlinear_remainder_treatment;

	ts_ln_erk_split_uv.setup(i_timestepping_order, true, true, true, true, false);

	// Setup semi-lag
	semiLagrangian.setup(ops.sphereDataConfig);

	if (timestepping_order != timestepping_order2)
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
	}
	else if  (timestepping_order == 4)
	{
		SWEETError("4th order method not (yet) supported");

#if 0
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size*0.5, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup(i_rexiSimVars, "phi1", i_timestep_size*0.5, false, true, timestepping_order);
		ts_phi2_exp.setup(i_rexiSimVars, "phi2", i_timestep_size*0.5, false, true, timestepping_order);

		// phi0, but with a full time step size
		ts_ups0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true, timestepping_order);
		ts_ups1_exp.setup(i_rexiSimVars, "ups1", i_timestep_size, false, true, timestepping_order);
		ts_ups2_exp.setup(i_rexiSimVars, "ups2", i_timestep_size, false, true, timestepping_order);
		ts_ups3_exp.setup(i_rexiSimVars, "ups3", i_timestep_size, false, true, timestepping_order);
#endif
	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}
}


SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		ops(i_op),

		ts_ln_erk_split_uv(simVars, ops),

		ts_phi0_exp(simVars, ops),
		ts_phi1_exp(simVars, ops),
		ts_phi2_exp(simVars, ops),

		semiLagrangian(simVars),
		sphereSampler(semiLagrangian.sphereSampler)

#if 0
		,
		ts_ups0_exp(simVars, ops),
		ts_ups1_exp(simVars, ops),
		ts_ups2_exp(simVars, ops),
		ts_ups3_exp(simVars, ops)
#endif
{
}



SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv::~SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv()
{
}

