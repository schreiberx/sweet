/*
 * SWE_Sphere_TS_l_na_settls_only.cpp
 *
 *  Created on: 2020-04-22
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2020-04-22: Based on ln_settls version, but just exactly this version: l implicit, na sl
 */

#include "SWE_Sphere_TS_l_irk_na_sl_settls_only.hpp"
#include <sweet/sphere/SphereData_DebugContainer.hpp>



void SWE_Sphere_TS_l_irk_na_sl_settls_only::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		FatalError("TODO run_timestep_1st_order_pert");
	}
	else if (timestepping_order == 2)
	{
		run_timestep_2nd_order_pert(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
	}
	else
	{
		FatalError("Only orders 1/2 supported (ERRID 098nd89eje)");
	}
}



void SWE_Sphere_TS_l_irk_na_sl_settls_only::run_timestep_2nd_order_pert(
	SphereData_Spectral &io_U_phi_pert,	///< prognostic variables
	SphereData_Spectral &io_U_vrt,		///< prognostic variables
	SphereData_Spectral &io_U_div,		///< prognostic variables

	double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
	double i_simulation_timestamp
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi_pert.sphereDataConfig;
//	double dt_div_radius = i_dt/simVars.sim.sphere_radius;

	if (i_dt <= 0)
		FatalError("SWE_Sphere_TS_l_na_settls_only: Only constant time step size allowed (Please set --dt)");

	const SphereData_Spectral &U_phi_pert = io_U_phi_pert;
	const SphereData_Spectral &U_vrt = io_U_vrt;
	const SphereData_Spectral &U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 * TODO: Maybe we should replace this with some 1st order SL time step or some subcycling
		 */
		U_phi_pert_prev = U_phi_pert;
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

	SphereData_Physical U_u_lon_prev(sphereDataConfig);
	SphereData_Physical U_v_lat_prev(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev, false);

	SphereData_Physical U_u_lon(sphereDataConfig);
	SphereData_Physical U_v_lat(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat, false);


	// Initialize departure points with arrival points for iterative corrections
	ScalarDataArray pos_lon_d(sphereDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray pos_lat_d(sphereDataConfig->physical_array_data_number_of_elements);


	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			U_u_lon_prev, U_v_lat_prev,
			U_u_lon, U_v_lat,

			i_dt,
			//i_simulation_timestamp,
			// doesn't matter
			-98123098120938102,
			simVars.sim.sphere_radius,
			nullptr,

			pos_lon_d, pos_lat_d,		// OUTPUT

			op,

			//simVars.disc.timestepping_order,
			// 2nd order
			2,

			//simVars.disc.semi_lagrangian_max_iterations,
			// 2 iterations
			2,

			//simVars.disc.semi_lagrangian_convergence_threshold,
			// no convergence threshold, 2 iterations only
			-1,

			//simVars.disc.semi_lagrangian_approximate_sphere_geometry,
			// no approximation of sphere geometry
			false,

			//simVars.disc.semi_lagrangian_interpolation_limiter
			// no limiters
			false
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
	SphereData_Spectral U_phi_pert_D(sphereDataConfig);
	SphereData_Spectral U_vrt_D(sphereDataConfig);
	SphereData_Spectral U_div_D(sphereDataConfig);

	SphereData_Physical U_phi_pert_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(
			U_phi_pert.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			U_phi_pert_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points
		);
	U_phi_pert_D = U_phi_pert_D_phys;

#if 0
	/*
	 * DO NOT USE VRT/DIV!!!
	 * THIS DOESN'T WORK!!!
	 */

	SphereData_Physical U_vort_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(U_vrt.getSphereDataPhysical(), pos_lon_d, pos_lat_d, U_vort_D_phys, false, simVars.disc.semi_lagrangian_interpolation_limiter);
	U_vrt_D = U_vort_D_phys;

	SphereData_Physical U_div_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(U_div.getSphereDataPhysical(), pos_lon_d, pos_lat_d, U_div_D_phys, false, simVars.disc.semi_lagrangian_interpolation_limiter);
	U_div_D = U_div_D_phys;

#else

	SphereData_Physical U_u_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(
			U_u_lon,
			pos_lon_d, pos_lat_d,
			U_u_D_phys,
			true,
			simVars.disc.semi_lagrangian_interpolation_limiter,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points
		);

	SphereData_Physical U_v_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(
			U_v_lat,
			pos_lon_d, pos_lat_d,
			U_v_D_phys,
			true,
			simVars.disc.semi_lagrangian_interpolation_limiter,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points
		);


	op.uv_to_vortdiv(U_u_D_phys, U_v_D_phys, U_vrt_D, U_div_D, false);

#endif

	/*
	 * Compute L_D
	 */

	SphereData_Spectral L_phi_pert_D(sphereDataConfig);
	SphereData_Spectral L_vort_D(sphereDataConfig);
	SphereData_Spectral L_div_D(sphereDataConfig);

	/*
	 * Method 1) First evaluate L, then sample result at departure point
	 */
	SphereData_Spectral L_phi_pert(sphereDataConfig);
	SphereData_Spectral L_vort(sphereDataConfig);
	SphereData_Spectral L_div(sphereDataConfig);

	// evaluate
	swe_sphere_ts_l_erk->euler_timestep_update(U_phi_pert, U_vrt, U_div, L_phi_pert, L_vort, L_div, i_simulation_timestamp);

	// sample phi
	SphereData_Physical L_phi_pert_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(L_phi_pert.getSphereDataPhysical(), pos_lon_d, pos_lat_d, L_phi_pert_D_phys, false, simVars.disc.semi_lagrangian_interpolation_limiter, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);
	L_phi_pert_D.loadSphereDataPhysical(L_phi_pert_D_phys);


	// sample u/v
	SphereData_Physical U_L_u_phys(sphereDataConfig);
	SphereData_Physical U_L_v_phys(sphereDataConfig);
	op.vortdiv_to_uv(L_vort, L_div, U_L_u_phys, U_L_v_phys, false);

	SphereData_Physical U_L_u_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(U_L_u_phys, pos_lon_d, pos_lat_d, U_L_u_D_phys, true, simVars.disc.semi_lagrangian_interpolation_limiter, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);

	SphereData_Physical U_L_v_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar(U_L_v_phys, pos_lon_d, pos_lat_d, U_L_v_D_phys, true, simVars.disc.semi_lagrangian_interpolation_limiter, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);

	op.uv_to_vortdiv(U_L_u_D_phys, U_L_v_D_phys, L_vort_D, L_div_D, false);


	/*
	 * Compute R = X_D + 1/2 dt L_D
	 */

	SphereData_Spectral R_phi_pert(sphereDataConfig);
	SphereData_Spectral R_vrt(sphereDataConfig);
	SphereData_Spectral R_div(sphereDataConfig);

	R_phi_pert = U_phi_pert_D + (0.5 * i_dt) * L_phi_pert_D;
	R_vrt = U_vrt_D + (0.5 * i_dt) * L_vort_D;
	R_div = U_div_D + (0.5 * i_dt) * L_div_D;


	/*
	 * Step 2b) Solve Helmholtz problem
	 * X - 1/2 dt LX = R
	 */
	swe_sphere_ts_l_irk->run_timestep_nonpert(
			R_phi_pert, R_vrt, R_div,
			0.5 * i_dt,
			i_simulation_timestamp
		);

	/*
	 * Backup previous variables for multi-step SL method
	 */
	U_phi_pert_prev = U_phi_pert;
	U_vrt_prev = U_vrt;
	U_div_prev = U_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_U_phi_pert = R_phi_pert;
	io_U_vrt = R_vrt;
	io_U_div = R_div;
}





void SWE_Sphere_TS_l_irk_na_sl_settls_only::setup(
		int i_timestepping_order
)
{
	std::cout << " + SWE_Sphere_TS_l_na_settls_only.setup() called" << std::endl;

	timestepping_order = i_timestepping_order;

	if (timestepping_order != 1 && timestepping_order != 2)
		FatalError("Invalid time stepping order, must be 1 or 2");

	// Setup sampler for future interpolations
	sphereSampler.setup(op.sphereDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig, simVars);

	// Initialize with 1st order
	swe_sphere_ts_l_erk = new SWE_Sphere_TS_l_erk(simVars, op);
	swe_sphere_ts_l_erk->setup(1);

	// Initialize with 1st order and half time step size
	swe_sphere_ts_l_irk = new SWE_Sphere_TS_l_irk(simVars, op);
	swe_sphere_ts_l_irk->setup(1, 0.5 * simVars.timecontrol.current_timestep_size, 0);
}



SWE_Sphere_TS_l_irk_na_sl_settls_only::SWE_Sphere_TS_l_irk_na_sl_settls_only(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		) :
		simVars(i_simVars), op(i_op),

		/*
		 * Initialize these fields as well to avoid any setup overheads during the 1st time step
		 */
		U_phi_pert_prev(i_op.sphereDataConfig),
		U_vrt_prev(i_op.sphereDataConfig),
		U_div_prev(i_op.sphereDataConfig)
{
}



SWE_Sphere_TS_l_irk_na_sl_settls_only::~SWE_Sphere_TS_l_irk_na_sl_settls_only()
{
	delete swe_sphere_ts_l_erk;
	delete swe_sphere_ts_l_irk;
}

