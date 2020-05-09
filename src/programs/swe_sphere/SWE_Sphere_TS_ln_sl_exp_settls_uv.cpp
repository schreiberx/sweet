/*
 * SWE_Sphere_TS_ln_exp_settls_uv.cpp
 *
 *  Created on: 31 Mar 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "SWE_Sphere_TS_ln_sl_exp_settls_uv.hpp"

#include <sweet/sphere/SphereData_DebugContainer.hpp>



void SWE_Sphere_TS_ln_sl_exp_settls_uv::run_timestep_pert(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 2)
	{
		run_timestep_2nd_order(io_phi, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		return;
	}
}



void SWE_Sphere_TS_ln_sl_exp_settls_uv::interpolate_departure_point_uv(
		const SphereData_Spectral &i_phi,
		const SphereData_Spectral &i_vrt,
		const SphereData_Spectral &i_div,

		const ScalarDataArray &i_pos_lon_d,
		const ScalarDataArray &i_pos_lat_d,

		SphereData_Spectral &o_phi,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
)
{
	o_phi = sphereSampler.bicubic_scalar_ret_phys(
			i_phi.getSphereDataPhysical(),
			i_pos_lon_d, i_pos_lat_d,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical u_tmp, v_tmp;
	op.vortdiv_to_uv(i_vrt, i_div, u_tmp, v_tmp);

	SphereData_Physical v_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
			u_tmp,
			i_pos_lon_d, i_pos_lat_d,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical u_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
			v_tmp,
			i_pos_lon_d, i_pos_lat_d,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	op.uv_to_vortdiv(u_tmp_D, v_tmp_D, o_vrt, o_div);
}



void SWE_Sphere_TS_ln_sl_exp_settls_uv::interpolate_departure_point_vd(
		const SphereData_Spectral &i_phi,
		const SphereData_Spectral &i_vrt,
		const SphereData_Spectral &i_div,

		const ScalarDataArray &i_pos_lon_d,
		const ScalarDataArray &i_pos_lat_d,

		SphereData_Spectral &o_phi,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
)
{
	o_phi = sphereSampler.bicubic_scalar_ret_phys(
			i_phi.getSphereDataPhysical(),
			i_pos_lon_d, i_pos_lat_d,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	o_vrt = sphereSampler.bicubic_scalar_ret_phys(
			i_vrt.getSphereDataPhysical(),
			i_pos_lon_d, i_pos_lat_d,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	o_div = sphereSampler.bicubic_scalar_ret_phys(
			i_div.getSphereDataPhysical(),
			i_pos_lon_d, i_pos_lat_d,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);
}



void SWE_Sphere_TS_ln_sl_exp_settls_uv::run_timestep_2nd_order(
		SphereData_Spectral &io_U_phi,	///< prognostic variables
		SphereData_Spectral &io_U_vrt,	///< prognostic variables
		SphereData_Spectral &io_U_div,	///< prognostic variables

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
	 * Step 1) Compute depature points
	 * Step 2) Compute U_D at departure points
	 * Step 3) Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - exp(dtL) N(t-dt) ]_D + N(t))
	 * Step 4) Compute update U(t+dt) = exp(dt L)(U_D + dt * N*(t+0.5dt))
	 */
	/*
	 *************************************************************************************************
	 * Step 1) Compute depature points
	 *************************************************************************************************
	 */
	SphereData_Physical U_u_lon_prev, U_v_lat_prev;
	op.vortdiv_to_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev);

	SphereData_Physical U_u_lon, U_v_lat;
	op.vortdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat);

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
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_pert_nr(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
			);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_pert_lc(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
				);
		}

		/*
		 * exp(dtL) N(t-dt)
		 */
		swe_sphere_ts_l_rexi->run_timestep_nonpert(
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
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_pert_nr(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp
			);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv->euler_timestep_update_pert_lc(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp-i_dt
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
	swe_sphere_ts_l_rexi->run_timestep_nonpert(
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
	SWE_Sphere_TS_ln_sl_exp_settls_uv::NLRemainderTreatment_enum nonlinear_divergence_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::NL_REMAINDER_IGNORE;

	bool original_linear_operator_sl_treatment = true;


	// Search for Coriolis
	if (simVars.disc.timestepping_method.find("l_irk") != std::string::npos || simVars.disc.timestepping_method.find("l_exp") != std::string::npos)
		linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_LINEAR;
	else if (simVars.disc.timestepping_method.find("lc_na_sl") != std::string::npos)
		linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_SEMILAGRANGIAN;
	else if (simVars.disc.timestepping_method.find("lc_") != std::string::npos)
		linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::CORIOLIS_NONLINEAR;


	// Search for Nonlinear divergence
	if (simVars.disc.timestepping_method.find("_nr_") != std::string::npos)
		nonlinear_divergence_treatment = SWE_Sphere_TS_ln_sl_exp_settls_uv::NL_REMAINDER_NONLINEAR;

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

	if (nonlinear_divergence_treatment == SWE_Sphere_TS_ln_sl_exp_settls_uv::NL_REMAINDER_NONLINEAR)
		string_id_storage += "_nr";

	string_id_storage += "_settls";

	if (!original_linear_operator_sl_treatment)
		string_id_storage += "_ver1";

	std::string string_id_storage_ = string_id_storage + "_uv";

	if (simVars.disc.timestepping_method != string_id_storage_)
	{
		if (!original_linear_operator_sl_treatment)
		{
			std::cerr << "Detected time stepping method: "+string_id_storage_ << std::endl;
			std::cerr << "Provided time stepping method: "+simVars.disc.timestepping_method << std::endl;
			FatalError("Autodetection of parts of time stepping methods failed!");
		}

		std::string string_id_storage2 = string_id_storage+"_ver0"+"_uv";
		if (simVars.disc.timestepping_method != string_id_storage2)
		{
			std::cerr << "Detected time stepping method: "+string_id_storage_ << std::endl;
			std::cerr << "Provided time stepping method: "+simVars.disc.timestepping_method << std::endl;
			std::cerr << "Detected alternative time stepping method: "+string_id_storage2 << std::endl;
			FatalError("Autodetection of parts of time stepping methods failed!");
		}
	}
#endif

	setup(
			simVars.disc.timestepping_order,
			linear_coriolis_treatment,					// Coriolis treatment
			nonlinear_divergence_treatment,		// Nonlinear divergence treatment
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
		FatalError("Invalid time stepping order, only 2nd order supported");

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig);

	swe_sphere_ts_ln_erk_split_uv = new SWE_Sphere_TS_ln_erk_split_uv(simVars, op);
	swe_sphere_ts_ln_erk_split_uv->setup(1, true, true, true, true);

	swe_sphere_ts_l_rexi = new SWE_Sphere_TS_l_exp(simVars, op);
	swe_sphere_ts_l_rexi->setup_new(
			simVars.rexi,
			"phi0",
			simVars.timecontrol.current_timestep_size,
			coriolis_treatment == CORIOLIS_LINEAR,
			simVars.sim.sphere_use_fsphere
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
	delete swe_sphere_ts_l_rexi;
}

