/*
 * SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv
 *
 * Created on: 24 Mar 2022
 *     Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#include "SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv.hpp"


bool SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::implements_timestepping_method(const std::string &i_timestepping_method)
{

	if (	i_timestepping_method == "lg_exp_na_sl_lc_nr_etdrk_uv"	||
		i_timestepping_method == "lg_exp_na_sl_lc_etdrk_uv"	||
			false
	)
		return true;

	return false;
}


std::string SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::string_id()
{
	return "lg_exp_na_sl_lc_nr_etdrk_uv";
}


void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::setup_auto()
{
	if (simVars.sim.sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	NLRemainderTreatment_enum nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;

	if (simVars.disc.timestepping_method == "lg_exp_na_sl_lc_nr_etdrk_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;
	}
	else if (simVars.disc.timestepping_method == "lg_exp_na_sl_lc_etdrk_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_IGNORE;
	}
	else
	{
		SWEETError(std::string("Timestepping method '")+simVars.disc.timestepping_method+std::string("' not known"));
	}

	setup(
			simVars.rexi,
			simVars.disc.timestepping_order,
			simVars.disc.timestepping_order2,
			simVars.timecontrol.current_timestep_size,

			nonlinear_remainder_treatment
		);
}



void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::print_help()
{
	std::cout << "	Exponential SL ETD:" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_nr_etdrk_uv" << std::endl;
}



void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::run_timestep(
		SphereData_Spectral &io_U_phi,	///< prognostic variables
		SphereData_Spectral &io_U_vrt,	///< prognostic variables
		SphereData_Spectral &io_U_div,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

	// Keep io data unchanged
	SphereData_Spectral U_phi = io_U_phi;
	SphereData_Spectral U_vrt = io_U_vrt;
	SphereData_Spectral U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
#if !SWEET_PARAREAL && !SWEET_XBRAID
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


	// N(Un) = Nonlinear term applied to Un
	SphereData_Spectral FUn_phi(sphereDataConfig, 0);
	SphereData_Spectral FUn_vrt(sphereDataConfig, 0);
	SphereData_Spectral FUn_div(sphereDataConfig, 0);

	// psi1 applied to N(Un)
	SphereData_Spectral psi1_FUn_phi(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_vrt(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_div(sphereDataConfig, 0);

	// Interpolated (.5 * psi1NUn + U) to departure points ( )_* //
	SphereData_Spectral psi1_FUn_phi_D(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_vrt_D(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_div_D(sphereDataConfig, 0);

	// All schemes use the result provided by ETD1RK
	// timestepping_order == -2 --> SL-ETD2RK-bis (double application of SL-ETD1RK; 2nd order scheme)
	if (timestepping_order == 1 || timestepping_order == 2 || timestepping_order == -2)
	{

		////////////////////
		// Compute N(U_n) //
		////////////////////
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi, U_vrt, U_div,
				FUn_phi, FUn_vrt, FUn_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi, U_vrt, U_div,
					FUn_phi, FUn_vrt, FUn_div,
					i_simulation_timestamp
			);
		}

		///////////////////////////
		// Apply psi_1 to N(U_n) //
		///////////////////////////
		ts_psi1_exp.run_timestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				i_dt,
				i_simulation_timestamp
			);


		////////////////////////////////////////////////////////
		// Add this to Un and interpolate to departure points //
		////////////////////////////////////////////////////////
		U_phi = U_phi + i_dt * psi1_FUn_phi;
		U_vrt = U_vrt + i_dt * psi1_FUn_vrt;
		U_div = U_div + i_dt * psi1_FUn_div;

		SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);

		//////////////////////////////////////////
		// Apply phi_0 to interpolated solution //
		//////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);

		ts_phi0_exp.run_timestep(
				U_phi_D, U_vrt_D, U_div_D,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;

	}
	// SL-ETD2RK-bis (2nd order)
	if (timestepping_order == -2)
	{

		////////////////////////////////////////////////////////
		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1} //
		////////////////////////////////////////////////////////
		// Save SL-ETD1RK from above
		SphereData_Spectral U_phi1(sphereDataConfig, 0);
		SphereData_Spectral U_vrt1(sphereDataConfig, 0);
		SphereData_Spectral U_div1(sphereDataConfig, 0);
		U_phi1 = U_phi;
		U_vrt1 = U_vrt;
		U_div1 = U_div;

		/////////////////////////////
		// Apply psi_1 to N(Un) //
		/////////////////////////////
		// already computed in the previous if block

		//////////////////////////////////////////////////////////////////////////
		// Add half of this to (original) U and interpolate to departure points //
		//////////////////////////////////////////////////////////////////////////
		U_phi = io_U_phi + .5 * i_dt * psi1_FUn_phi;
		U_vrt = io_U_vrt + .5 * i_dt * psi1_FUn_vrt;
		U_div = io_U_div + .5 * i_dt * psi1_FUn_div;

		/////////////////////////////////////////////////////////////
		//Interpolate (.5 * psi1NUn + U) to departure points ( )_* //
		/////////////////////////////////////////////////////////////
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				psi1_FUn_phi_D, psi1_FUn_vrt_D, psi1_FUn_div_D
			);

		////////////////////////
		// Compute N(U_{n+1}) //
		////////////////////////
		SphereData_Spectral FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral FU1_div(sphereDataConfig, 0);
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi1, U_vrt1, U_div1,
				FU1_phi, FU1_vrt, FU1_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi1, U_vrt1, U_div1,
					FU1_phi, FU1_vrt, FU1_div,
					i_simulation_timestamp
			);
		}

		//////////////////////////////
		// Apply psi1 to N(U_{n+1}) //
		//////////////////////////////
		SphereData_Spectral psi1_FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral psi1_FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral psi1_FU1_div(sphereDataConfig, 0);
		ts_psi1_exp.run_timestep(
				FU1_phi, FU1_vrt, FU1_div,
				psi1_FU1_phi, psi1_FU1_vrt, psi1_FU1_div,
				i_dt,
				i_simulation_timestamp
			);

		/////////////////////////////////////////////////////
		// Add half of this to already interpolated values //
		/////////////////////////////////////////////////////
		U_phi = psi1_FUn_phi_D + .5 * i_dt * psi1_FU1_phi;
		U_vrt = psi1_FUn_vrt_D + .5 * i_dt * psi1_FU1_vrt;
		U_div = psi1_FUn_div_D + .5 * i_dt * psi1_FU1_div;

		////////////////////////////////////////////////////////////////////////
		// Calculate phi_0 of [interpolated (.5 * psi1NUN + U) + .5 * ps1NU1] //
		////////////////////////////////////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
		);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;

	}
	else if (timestepping_order == 2)
	{

		////////////////////////////////////////////////////////
		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1} //
		////////////////////////////////////////////////////////
		// Save SL-ETD1RK from above
		SphereData_Spectral U_phi1(sphereDataConfig, 0);
		SphereData_Spectral U_vrt1(sphereDataConfig, 0);
		SphereData_Spectral U_div1(sphereDataConfig, 0);
		U_phi1 = U_phi;
		U_vrt1 = U_vrt;
		U_div1 = U_div;

		////////////////////////
		// Compute N(U_{n+1}) //
		////////////////////////
		SphereData_Spectral FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral FU1_div(sphereDataConfig, 0);
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi1, U_vrt1, U_div1,
				FU1_phi, FU1_vrt, FU1_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi1, U_vrt1, U_div1,
					FU1_phi, FU1_vrt, FU1_div,
					i_simulation_timestamp
			);
		}

		//////////////////////////////
		// Apply psi2 to N(U_{n+1}) //
		//////////////////////////////
		SphereData_Spectral psi2_FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral psi2_FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral psi2_FU1_div(sphereDataConfig, 0);
		ts_psi2_exp.run_timestep(
				FU1_phi, FU1_vrt, FU1_div,
				psi2_FU1_phi, psi2_FU1_vrt, psi2_FU1_div,
				i_dt,
				i_simulation_timestamp
		);

		//////////////////////////////
		// Apply psi2 to N(U_{n}) //
		//////////////////////////////
		SphereData_Spectral psi2_FUn_phi(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_vrt(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_div(sphereDataConfig, 0);
		ts_psi2_exp.run_timestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi2_FUn_phi, psi2_FUn_vrt, psi2_FUn_div,
				i_dt,
				i_simulation_timestamp
		);

		////////////////////////////////////////////////////
		// Interpolate psi2NUn to departure points ( )_* //
		////////////////////////////////////////////////////
		SphereData_Spectral psi2_FUn_phi_D(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_vrt_D(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_div_D(sphereDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				psi2_FUn_phi, psi2_FUn_vrt, psi2_FUn_div,
				pos_lon_d, pos_lat_d,
				psi2_FUn_phi_D, psi2_FUn_vrt_D, psi2_FUn_div_D
			);

		////////////////////////////////////////
		// Apply phi0 to psi2NU1 - (psi2NUn)_* //
		////////////////////////////////////////
		SphereData_Spectral phi0_dif2_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_dif2_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_dif2_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				psi2_FU1_phi - psi2_FUn_phi_D, psi2_FU1_vrt - psi2_FUn_vrt_D, psi2_FU1_div - psi2_FUn_div_D,
				phi0_dif2_phi, phi0_dif2_vrt, phi0_dif2_div,
				i_dt,
				i_simulation_timestamp
		);

		U_phi = U_phi1 + i_dt*phi0_dif2_phi;
		U_vrt = U_vrt1 + i_dt*phi0_dif2_vrt;
		U_div = U_div1 + i_dt*phi0_dif2_div;

	}

	// Save current time step for next step
	U_phi_prev = io_U_phi;
	U_vrt_prev = io_U_vrt;
	U_div_prev = io_U_div;

	io_U_phi = U_phi;
	io_U_vrt = U_vrt;
	io_U_div = U_div;

}


/*
 * Setup
 */
void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::setup(
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

	if (timestepping_order == 0 || timestepping_order == 1 || timestepping_order == -2)
	{
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true, timestepping_order);	/* NO Coriolis */
		ts_psi1_exp.setup(i_rexiSimVars, "psi1", i_timestep_size, false, true, timestepping_order);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true, timestepping_order);	/* NO Coriolis */
		ts_psi1_exp.setup(i_rexiSimVars, "psi1", i_timestep_size, false, true, timestepping_order);
		ts_psi2_exp.setup(i_rexiSimVars, "psi2", i_timestep_size, false, true, timestepping_order);
	}
	else if  (timestepping_order == 4)
	{
		SWEETError("4th order method not (yet) supported");

	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}
}


SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		ops(i_op),

		ts_ln_erk_split_uv(simVars, ops),

		ts_phi0_exp(simVars, ops),
		ts_phi2_exp(simVars, ops),

		ts_psi1_exp(simVars, ops),
		ts_psi2_exp(simVars, ops),

		semiLagrangian(simVars),
		sphereSampler(semiLagrangian.sphereSampler)

{
}



SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::~SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv()
{
}

