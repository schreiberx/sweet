/*
 * SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only.cpp
 *
 *  Created on: 2020-04-29
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only.hpp"



bool SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::implements_timestepping_method(const std::string &i_timestepping_method)
{
	if (i_timestepping_method == "l_irk_na_sl_nr_settls_uv_only")
		return true;

	return false;
}

std::string SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::string_id()
{
	return "l_irk_na_sl_nr_settls_uv_only";
}


void SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::setup_auto()
{
	setup(
		simVars.disc.timestepping_order
	);
}


void SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		SWEETError("TODO run_timestep_1st_order_pert");
	}
	else if (timestepping_order == 2)
	{
		run_timestep_2nd_order_pert(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
	}
	else
	{
		SWEETError("Only orders 1/2 supported (ERRID 098nd89eje)");
	}
}



void SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::interpolate_departure_point_uv(
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
	o_vrt.setup_if_required(i_phi.sphereDataConfig);
	o_div.setup_if_required(i_phi.sphereDataConfig);

	o_phi = sphereSampler.bicubic_scalar_ret_phys(
			i_phi.getSphereDataPhysical(),
			i_pos_lon_D, i_pos_lat_D,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);


	SphereData_Physical u_tmp, v_tmp;
	op.vortdiv_to_uv(i_vrt, i_div, u_tmp, v_tmp);

	SphereData_Physical u_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
			u_tmp,
			i_pos_lon_D, i_pos_lat_D,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical v_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
			v_tmp,
			i_pos_lon_D, i_pos_lat_D,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

#if 1

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
	ScalarDataArray V_lon_D = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(u_tmp_D);
	ScalarDataArray V_lat_D = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(v_tmp_D);

	ScalarDataArray V_x_D, V_y_D, V_z_D;
	SWEETMath::latlon_velocity_to_cartesian_velocity(
			i_pos_lon_D,
			i_pos_lat_D,
			V_lon_D,
			V_lat_D,
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

	op.uv_to_vortdiv(
			Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_lon_A, i_vrt.sphereDataConfig),
			Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_lat_A, i_vrt.sphereDataConfig),
			o_vrt, o_div
		);

#else

	op.uv_to_vortdiv(u_tmp_D, v_tmp_D, o_vrt, o_div);

#endif
}



void SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::run_timestep_2nd_order_pert(
	SphereData_Spectral &io_U_phi,		///< prognostic variables
	SphereData_Spectral &io_U_vrt,		///< prognostic variables
	SphereData_Spectral &io_U_div,		///< prognostic variables

	double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
	double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Sphere_TS_l_na_settls_only: Only constant time step size allowed (Please set --dt)");

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
	 * Step 1) SL
	 * Compute Lagrangian trajectories based on SETTLS.
	 * This depends on V(t-\Delta t) and V(t).
	 *
	 * See Hortal's paper for equation.
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
	 * Step 2) Midpoint rule
	 * Put everything together with midpoint rule and solve resulting Helmholtz problem
	 */

	/*
	 * Step 2a) Compute RHS
	 * R = X_D + 1/2 dt L_D + dt N*
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
	 * Compute L_D
	 */
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

	/*
	 * Method 1) First evaluate L, then sample result at departure point
	 */
	SphereData_Spectral L_U_phi(sphereDataConfig, 0), L_U_vrt(sphereDataConfig, 0), L_U_div(sphereDataConfig, 0);

	/*
	 * L_g(U): Linear gravity modes
	 */
	swe_sphere_ts_ln_erk_split_uv__l_erk_1st_order->euler_timestep_update_pert_lg(
			U_phi, U_vrt, U_div,
			L_U_phi, L_U_vrt, L_U_div,
			i_simulation_timestamp
		);

	/*
	 * L_c(U): Linear Coriolis effect
	 */
	swe_sphere_ts_ln_erk_split_uv__l_erk_1st_order->euler_timestep_update_pert_lc(
			U_phi, U_vrt, U_div,
			L_U_phi, L_U_vrt, L_U_div,
			i_simulation_timestamp
		);


	SphereData_Spectral L_U_phi_D, L_U_vrt_D, L_U_div_D;
	interpolate_departure_point_uv(
			L_U_phi, L_U_vrt, L_U_div,
			pos_lon_d, pos_lat_d,
			L_U_phi_D, L_U_vrt_D, L_U_div_D
		);

	/*
	 * Compute R = X_D + 1/2 dt L_D
	 */
	SphereData_Spectral R_phi = U_phi_D + (0.5 * i_dt) * L_U_phi_D;
	SphereData_Spectral R_vrt = U_vrt_D + (0.5 * i_dt) * L_U_vrt_D;
	SphereData_Spectral R_div = U_div_D + (0.5 * i_dt) * L_U_div_D;

	/*
	 * Nonlinear remainder term starts here
	 */
	if (1)
	{
		/*
		 * N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 *
		 * R += dt*N*(t+0.5dt)
		 */

		/*
		 * Compute
		 * [ 2*N(t) - N(t-dt) ]_D
		 */

		/*
		 * N(t-dt)
		 */
		SphereData_Spectral N_U_phi_prev_nr(sphereDataConfig, 0);
		SphereData_Spectral N_U_vrt_prev_nr(sphereDataConfig, 0);
		SphereData_Spectral N_U_div_prev_nr(sphereDataConfig, 0);

		swe_sphere_ts_ln_erk_split_uv__l_erk_1st_order->euler_timestep_update_pert_nr(
				U_phi_prev, U_vrt_prev, U_div_prev,
				N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
				i_simulation_timestamp-i_dt
		);

		/*
		 * N(t)
		 */
		SphereData_Spectral N_U_phi_nr(sphereDataConfig, 0);
		SphereData_Spectral N_U_vrt_nr(sphereDataConfig, 0);
		SphereData_Spectral N_U_div_nr(sphereDataConfig, 0);

		swe_sphere_ts_ln_erk_split_uv__l_erk_1st_order->euler_timestep_update_pert_nr(
				U_phi, U_vrt, U_div,
				N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
				i_simulation_timestamp
		);


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
		 * Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 * and add to R terms
		 */
		R_phi += (i_dt * 0.5) * (N_U_phi_next_D + N_U_phi_nr);
		R_vrt += (i_dt * 0.5) * (N_U_vrt_next_D + N_U_vrt_nr);
		R_div += (i_dt * 0.5) * (N_U_div_next_D + N_U_div_nr);
	}

	/*
	 * Step 2b) Solve Helmholtz problem
	 * X - 1/2 dt LX = R
	 */
	swe_sphere_ts_l_irk->run_timestep_pert(
			R_phi, R_vrt, R_div,
			0.5 * i_dt,
			i_simulation_timestamp
		);

	/*
	 * Backup previous variables for multi-step SL method
	 */
	U_phi_prev = U_phi;
	U_vrt_prev = U_vrt;
	U_div_prev = U_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_U_phi = R_phi;
	io_U_vrt = R_vrt;
	io_U_div = R_div;
}




void SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::setup(
		int i_timestepping_order
)
{
	timestepping_order = i_timestepping_order;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETError("Invalid time stepping order, must be 1 or 2");

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig);

	// Initialize with 1st order
	swe_sphere_ts_ln_erk_split_uv__l_erk_1st_order = new SWE_Sphere_TS_ln_erk_split_uv(simVars, op);
	swe_sphere_ts_ln_erk_split_uv__l_erk_1st_order->setup(1, true, true, false, false);

	// Initialize with 1st order and half time step size
	swe_sphere_ts_l_irk = new SWE_Sphere_TS_l_irk(simVars, op);
	swe_sphere_ts_l_irk->setup(1, 0.5 * simVars.timecontrol.current_timestep_size, 0);
}



SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only(
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



SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::~SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only()
{
	delete swe_sphere_ts_ln_erk_split_uv__l_erk_1st_order;
	delete swe_sphere_ts_l_irk;
}

