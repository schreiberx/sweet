/*
 * SWE_Sphere_TS_ln_exp_settls_uv.cpp
 *
 *  Created on: 31 Mar 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_ln_sl_exp_settls_uv.hpp"


bool SWE_Sphere_TS_ln_sl_exp_settls_uv::implements_timestepping_method(const std::string &i_timestepping_method)
{
	/*
	 * Should contain _exp and _settls as well as _uv to indicate vorticity-divergence formulation
	 */
	return (
		(i_timestepping_method.find("_exp") != std::string::npos)		&&
		(i_timestepping_method.find("_settls") != std::string::npos)	&&
		(i_timestepping_method.find("_uv") != std::string::npos)		&&
		!(i_timestepping_method.find("_only") != std::string::npos)		&&
		true
	);

	return false;
}


std::string SWE_Sphere_TS_ln_sl_exp_settls_uv::string_id()
{
	return string_id_storage;
}


void SWE_Sphere_TS_ln_sl_exp_settls_uv::run_timestep(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,				///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 2)
	{
		run_timestep_2nd_order(io_phi, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		return;
	}

	SWEETError("Only 2nd order TS method supported!");
}



void SWE_Sphere_TS_ln_sl_exp_settls_uv::interpolate_departure_point_uv(
		const SphereData_Spectral &i_phi,
		const SphereData_Spectral &i_vrt,
		const SphereData_Spectral &i_div,

		const ScalarDataArray &i_pos_lon_D,
		const ScalarDataArray &i_pos_lat_D,

		SphereData_Spectral &o_phi,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
)
{
	o_phi.setup_if_required(i_phi.sphereDataConfig);
	o_vrt.setup_if_required(i_vrt.sphereDataConfig);
	o_div.setup_if_required(i_div.sphereDataConfig);

	o_phi = sphereSampler.bicubic_scalar_ret_phys(
			i_phi.toPhys(),
			i_pos_lon_D, i_pos_lat_D,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);


	SphereData_Physical U_u_A, U_v_A;
	op.vrtdiv_to_uv(i_vrt, i_div, U_u_A, U_v_A);

	SphereData_Physical U_u_D = sphereSampler.bicubic_scalar_ret_phys(
			U_u_A,
			i_pos_lon_D, i_pos_lat_D,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical U_v_D = sphereSampler.bicubic_scalar_ret_phys(
			U_v_A,
			i_pos_lon_D, i_pos_lat_D,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

#if 1

	/*
	 * Compute vectors to arrival and departure points
	 */
	ScalarDataArray P_x_D, P_y_D, P_z_D;
	SWEETMath::latlon_to_cartesian(i_pos_lon_D, i_pos_lat_D, P_x_D, P_y_D, P_z_D);

	ScalarDataArray	&P_x_A = semiLagrangian.pos_x_A,
					&P_y_A = semiLagrangian.pos_y_A,
					&P_z_A = semiLagrangian.pos_z_A;


	/*
	 * Compute rotation angle
	 */
	ScalarDataArray rotation_angle_ =
			SWEETMath::dot_prod(
				P_x_D, P_y_D, P_z_D,
				P_x_A, P_y_A, P_z_A
			);

	// Can be slightly larger than 1, leading to NaN, hence this hack
	rotation_angle_ = SWEETMath::min(rotation_angle_, 1.0);

	ScalarDataArray rotation_angle = SWEETMath::arccos(rotation_angle_);


	/*
	 * Compute Rotation axis
	 */
	ScalarDataArray rot_x, rot_y, rot_z;
	SWEETMath::cross_prod(
			P_x_D, P_y_D, P_z_D,
			P_x_A, P_y_A, P_z_A,
			rot_x, rot_y, rot_z
		);

	SWEETMath::normalize_with_threshold(rot_x, rot_y, rot_z);


	/*
	 * Convert to Cartesian velocity space
	 */
	ScalarDataArray U_u_D_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(U_u_D);
	ScalarDataArray V_v_D_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(U_v_D);

	ScalarDataArray V_x_D, V_y_D, V_z_D;
	SWEETMath::latlon_velocity_to_cartesian_velocity(
			i_pos_lon_D,
			i_pos_lat_D,
			U_u_D_array,
			V_v_D_array,
			V_x_D,
			V_y_D,
			V_z_D
		);

	/*
	 * Rotate to velocity vector
	 */
	ScalarDataArray V_x_A, V_y_A, V_z_A;
	SWEETMath::rotate_3d_vector_normalized_rotation_axis(
			V_x_D, V_y_D, V_z_D,
			rotation_angle,
			rot_x, rot_y, rot_z,
			V_x_A, V_y_A, V_z_A
		);

#if 0
	if (V_x_A.reduce_isAnyNaNorInf())
		std::cout << "V_x_A" << std::endl;

	if (V_y_A.reduce_isAnyNaNorInf())
		std::cout << "V_y_A" << std::endl;

	if (V_z_A.reduce_isAnyNaNorInf())
		std::cout << "V_z_A" << std::endl;
#endif

	/*
	 * Return velocity in lat/lon space
	 */
	ScalarDataArray V_lon_A, V_lat_A;
	SWEETMath::cartesian_velocity_to_latlon_velocity(
			semiLagrangian.pos_lon_A,
			semiLagrangian.pos_lat_A,
			V_x_A,
			V_y_A,
			V_z_A,
			V_lon_A, V_lat_A
	);

	op.uv_to_vrtdiv(
			Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_lon_A, i_vrt.sphereDataConfig),
			Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_lat_A, i_vrt.sphereDataConfig),
			o_vrt, o_div
		);

#else

	op.uv_to_vrtdiv(U_u_D, U_v_D, o_vrt, o_div);

#endif
}




void SWE_Sphere_TS_ln_sl_exp_settls_uv::run_timestep_2nd_order(
		SphereData_Spectral &io_U_phi,
		SphereData_Spectral &io_U_vrt,
		SphereData_Spectral &io_U_div,

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp)
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
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}

	/*
	 * Step 1) Compute departure points
	 *
	 * Step 2) Compute U_D at departure points
	 *
	 * Step 3) Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - exp(dtL) N(t-dt) ]_D + N(t))
	 *
	 * Step 4) Compute update U(t+dt) = exp(dt L)(U_D + dt * N*(t+0.5dt))
	 */

	/*
	 *************************************************************************************************
	 * Step 1) Compute depature points
	 *************************************************************************************************
	 */
	SphereData_Physical U_u_lon_prev, U_v_lat_prev;
	op.vrtdiv_to_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev);

	SphereData_Physical U_u_lon, U_v_lat;
	op.vrtdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat);

	double dt_div_radius = i_dt / simVars.sim.sphere_radius;

	// Calculate departure points
	ScalarDataArray pos_lon_d, pos_lat_d;
	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_lon_prev, dt_div_radius*U_v_lat_prev,
			dt_div_radius*U_u_lon, dt_div_radius*U_v_lat,
			pos_lon_d, pos_lat_d		// OUTPUT
		);


	/*
	 *************************************************************************************************
	 * Step 2) Compute U_D at departure points
	 *************************************************************************************************
	 */

	/*
	 * Compute X_D
	 */
	SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
	interpolate_departure_point_uv(
			U_phi, U_vrt, U_div,
			pos_lon_d, pos_lat_d,
			U_phi_D, U_vrt_D, U_div_D
		);


	/*
	 *************************************************************************************************
	 * Step 3) Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - exp(dtL) N(t-dt) ]_D + N(t))
	 *************************************************************************************************
	 */
	{
		/*
		 * N*(t+0.5dt) = 1/2 ([ 2*N(t) - exp(dtL) N(t-dt) ]_D + N(t))
		 */

		/*
		 * N(t-dt)
		 */
		const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;
		SphereData_Spectral N_U_phi_prev_nr(sphereDataConfig, 0), N_U_vrt_prev_nr(sphereDataConfig, 0), N_U_div_prev_nr(sphereDataConfig, 0);

		if (nonlinear_remainder_treatment == NL_REMAINDER_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_nr(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
			);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_lc(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
				);
		}

		/*
		 * exp(dtL) N(t-dt)
		 */
		swe_sphere_ts_l_rexi->run_timestep(
			N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
			i_dt,
			i_simulation_timestamp
		);


		/*
		 * N(t)
		 */
		SphereData_Spectral N_U_phi_nr(sphereDataConfig, 0);
		SphereData_Spectral N_U_vrt_nr(sphereDataConfig, 0);
		SphereData_Spectral N_U_div_nr(sphereDataConfig, 0);

		if (nonlinear_remainder_treatment == NL_REMAINDER_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_nr(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp
				);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp
				);
		}

		/*
		 * N(t+dt)_D = [ 2*N(t) - N(t-dt) ]_D
		 */
		SphereData_Spectral N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D;
		interpolate_departure_point_uv(
				2.0 * N_U_phi_nr - N_U_phi_prev_nr,
				2.0 * N_U_vrt_nr - N_U_vrt_prev_nr,
				2.0 * N_U_div_nr - N_U_div_prev_nr,

				pos_lon_d, pos_lat_d,
				N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D
			);

		/*
		 * Compute
		 * 		N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 * and
		 * 		U_D + dt * N*(t+0.5dt)
		 */
		U_phi_D += (i_dt * 0.5) * (N_U_phi_next_D + N_U_phi_nr);
		U_vrt_D += (i_dt * 0.5) * (N_U_vrt_next_D + N_U_vrt_nr);
		U_div_D += (i_dt * 0.5) * (N_U_div_next_D + N_U_div_nr);
	}


	/*
	 *************************************************************************************************
	 * Step 4) Compute update U(t+dt) = exp(dt L)(U_D + dt * N*(t+0.5dt))
	 *************************************************************************************************
	 */
	swe_sphere_ts_l_rexi->run_timestep(
		U_phi_D, U_vrt_D, U_div_D,
		i_dt,
		i_simulation_timestamp
	);

	/*
	 * Backup previous variables for multi-step SL method
	 */
	U_phi_prev = io_U_phi;
	U_vrt_prev = io_U_vrt;
	U_div_prev = io_U_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_U_phi = U_phi_D;
	io_U_vrt = U_vrt_D;
	io_U_div = U_div_D;
}


void SWE_Sphere_TS_ln_sl_exp_settls_uv::setup_auto()

{
	SWE_Sphere_TS_ln_sl_exp_settls_uv::LinearCoriolisTreatment_enum linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_IGNORE;
	SWE_Sphere_TS_ln_sl_exp_settls_uv::NLRemainderTreatment_enum nonlinear_remainder_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::NL_REMAINDER_IGNORE;

	bool original_linear_operator_sl_treatment = true;


	// Search for Coriolis
	if (simVars.disc.timestepping_method.find("l_exp") != std::string::npos)
		linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_LINEAR;
	else if (simVars.disc.timestepping_method.find("lc_na_sl") != std::string::npos)
		linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_SEMILAGRANGIAN;
	else if (simVars.disc.timestepping_method.find("lc_") != std::string::npos)
		linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_NONLINEAR;


	// Search for nonlinear remainder term
	if (simVars.disc.timestepping_method.find("_nr_") != std::string::npos)
		nonlinear_remainder_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::NL_REMAINDER_NONLINEAR;

	if (simVars.disc.timestepping_method.find("_ver1") != std::string::npos)
		original_linear_operator_sl_treatment = false;
	else if (simVars.disc.timestepping_method.find("_ver0") != std::string::npos)
		original_linear_operator_sl_treatment = true;
	else
		original_linear_operator_sl_treatment = true;

#if 1
	string_id_storage = "";

	if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_LINEAR)
		string_id_storage += "l";
	else
		string_id_storage += "lg";

	string_id_storage += "_exp";

	if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_SEMILAGRANGIAN)
		string_id_storage += "_lc";

	string_id_storage += "_na";

	string_id_storage += "_sl";

	if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_NONLINEAR)
		string_id_storage += "_lc";

	if (nonlinear_remainder_treatment == SWE_Sphere_TS_ln_sl_exp_settls_uv::NL_REMAINDER_NONLINEAR)
		string_id_storage += "_nr";

	string_id_storage += "_settls";

	if (!original_linear_operator_sl_treatment)
		string_id_storage += "_ver1";

	std::string string_id_storage_ = string_id_storage + "_uv";

	if (simVars.disc.timestepping_method != string_id_storage_)
	{
		if (!original_linear_operator_sl_treatment)
		{
			std::cerr << "Provided time stepping method: "+simVars.disc.timestepping_method << std::endl;
			std::cerr << "Detected time stepping method: "+string_id_storage_ << std::endl;
			SWEETError("Autodetection of parts of time stepping methods failed!");
		}

		std::string string_id_storage2 = string_id_storage+"_ver0"+"_uv";
		if (simVars.disc.timestepping_method != string_id_storage2)
		{
			std::cerr << "Provided time stepping method: "+simVars.disc.timestepping_method << std::endl;
			std::cerr << "Detected time stepping method (failed): "+string_id_storage_ << std::endl;
			std::cerr << "Detected alternative time stepping method (failed): "+string_id_storage2 << std::endl;
			SWEETError("Autodetection of parts of time stepping methods failed!");
		}
	}
#endif

	setup(
			simVars.disc.timestepping_order,
			linear_coriolis_treatment,					// Coriolis treatment
			nonlinear_remainder_treatment,		// Nonlinear divergence treatment
			original_linear_operator_sl_treatment			// original SL linear operator treatment
		);
}



void SWE_Sphere_TS_ln_sl_exp_settls_uv::setup(
		int i_timestepping_order,
		LinearCoriolisTreatment_enum i_coriolis_treatment,
		NLRemainderTreatment_enum i_nonlinear_remainder_treatment,
		bool i_original_linear_operator_sl_treatment
)
{
	std::cout << " + SWE_Sphere_TS_ln_sl_exp_settls.setup() called" << std::endl;

	coriolis_treatment = i_coriolis_treatment;
	nonlinear_remainder_treatment = i_nonlinear_remainder_treatment;
	timestepping_order = i_timestepping_order;
	original_linear_operator_sl_treatment = i_original_linear_operator_sl_treatment;


	if (timestepping_order != 2)
		SWEETError("Invalid time stepping order, only 2nd order supported");

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig);

	swe_sphere_ts_ln_erk_split_uv = new SWE_Sphere_TS_ln_erk_split_uv(simVars, op);
	swe_sphere_ts_ln_erk_split_uv->setup(1, true, true, true, true, false);

	swe_sphere_ts_l_rexi = new SWE_Sphere_TS_l_exp(simVars, op);
	swe_sphere_ts_l_rexi->setup(
			simVars.rexi,
			"phi0",
			simVars.timecontrol.current_timestep_size,
			simVars.sim.sphere_use_fsphere,
			!(coriolis_treatment == CORIOLIS_LINEAR)
		);
}

SWE_Sphere_TS_ln_sl_exp_settls_uv::SWE_Sphere_TS_ln_sl_exp_settls_uv(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op,
			bool i_setup_auto
		) :
		simVars(i_simVars),
		op(i_op),
		semiLagrangian(simVars),
		sphereSampler(semiLagrangian.sphereSampler)
{
	if (i_setup_auto)
		setup_auto();
}

SWE_Sphere_TS_ln_sl_exp_settls_uv::~SWE_Sphere_TS_ln_sl_exp_settls_uv()
{
	delete swe_sphere_ts_ln_erk_split_uv;
	delete swe_sphere_ts_l_rexi;
}

