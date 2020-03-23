/*
 * SWE_Sphere_TS_ln_settls.cpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2019-10-24: Partly based on plane version
 */

#include "SWE_Sphere_TS_ln_settls.hpp"


/**
 * SETTLS implementation (see Hortal 2002)
 *
 * Solve SWE with Crank-Nicolson implicit time stepping
 * (spectral formulation for Helmholtz equation) with semi-Lagrangian
 * SL-SI-SP
 *
 * U_t = L U(0)
 *
 * Fully implicit version:
 *
 * (U(tau) - U(0)) / tau = 0.5*(L U(tau) + L U(0))
 *
 * <=> U(tau) - U(0) =  tau * 0.5*(L U(tau) + L U(0))
 *
 * <=> U(tau) - 0.5* L tau U(tau) = U(0) + tau * 0.5*L U(0)
 *
 * <=> (1 - 0.5 L tau) U(tau) = (1 + tau*0.5*L) U(0)
 *
 * <=> (2/tau - L) U(tau) = (2/tau + L) U(0)
 *
 * <=> U(tau) = (2/tau - L)^{-1} (2/tau + L) U(0)
 *
 * Semi-implicit has Coriolis term as totally explicit
 *
 * Semi-Lagrangian:
 *   U(tau) is on arrival points
 *   U(0) is on departure points
 *
 * Nonlinear term is added following Hortal (2002)
 * http://onlinelibrary.wiley.com/doi/10.1002/qj.200212858314/pdf
 */

void SWE_Sphere_TS_ln_settls::run_timestep(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Sphere_TS_ln_settls: Only constant time step size allowed (Please set --dt)");

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 * TODO: Maybe we should replace this with some 1st order SL time step or some subcycling
		 */
		phi_prev = io_phi;
		vort_prev = io_vort;
		div_prev = io_div;
	}

	// Output variables
	SphereData_Spectral phi(io_phi.sphereDataConfig);
	SphereData_Spectral vort(io_phi.sphereDataConfig);
	SphereData_Spectral div(io_phi.sphereDataConfig);

	/*
	 * Step 1) SL
	 * Compute Lagrangian trajectories based on SETTLS.
	 * This depends on V(t-\Delta t) and V(t).
	 * See Hortal for equation.
	 */

	// Departure points and arrival points
	ScalarDataArray pos_lon_d = pos_lon_a;
	ScalarDataArray pos_lat_d = pos_lat_a;

	SphereData_Physical u_lon_prev(io_phi.sphereDataConfig);
	SphereData_Physical v_lat_prev(io_phi.sphereDataConfig);
	op.vortdiv_to_uv(vort_prev, div_prev, u_lon_prev, v_lat_prev);

	SphereData_Physical u_lon(io_phi.sphereDataConfig);
	SphereData_Physical v_lat(io_phi.sphereDataConfig);
	op.vortdiv_to_uv(io_vort, io_div, u_lon, v_lat);

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			u_lon_prev,	v_lat_prev,
			u_lon, v_lat,
			pos_lon_a,	pos_lat_a,

			i_dt,
			simVars.sim.sphere_radius,

			pos_lon_d,	pos_lat_d,		// OUTPUT

			simVars.disc.timestepping_order,
			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold,

			simVars.disc.semi_lagrangian_approximate_sphere_geometry,
			simVars.disc.semi_lagrangian_interpolation_limiter
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
	SphereData_Physical phi_D_phys(io_phi.sphereDataConfig);
	SphereData_Physical vort_D_phys(io_phi.sphereDataConfig);
	SphereData_Physical div_D_phys(io_phi.sphereDataConfig);

	sphereSampler.bicubic_scalar(
			io_phi.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			phi_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	sphereSampler.bicubic_scalar(
			io_vort.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			vort_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	sphereSampler.bicubic_scalar(
			io_div.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			div_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Spectral phi_D(io_phi.sphereDataConfig);
	SphereData_Spectral vort_D(io_phi.sphereDataConfig);
	SphereData_Spectral div_D(io_phi.sphereDataConfig);

	phi_D.loadSphereDataPhysical(phi_D_phys);
	vort_D.loadSphereDataPhysical(vort_D_phys);
	div_D.loadSphereDataPhysical(div_D_phys);


	/*
	 * Compute L_D
	 */

	SphereData_Spectral L_phi_D(io_phi.sphereDataConfig);
	SphereData_Spectral L_vort_D(io_phi.sphereDataConfig);
	SphereData_Spectral L_div_D(io_phi.sphereDataConfig);

	if (original_linear_operator_sl_treatment)
	{
		/*
		 * Method 1) First evaluate L, then sample result at departure point
		 */
		SphereData_Spectral L_phi(io_phi.sphereDataConfig);
		SphereData_Spectral L_vort(io_phi.sphereDataConfig);
		SphereData_Spectral L_div(io_phi.sphereDataConfig);

		if (coriolis_treatment == CORIOLIS_LINEAR)
		{
			swe_sphere_ts_l_erk->euler_timestep_update(
					io_phi, io_vort, io_div,
					L_phi, L_vort, L_div,
					i_simulation_timestamp
			);
		}
		else
		{
			swe_sphere_ts_lg_erk->euler_timestep_update(
					io_phi, io_vort, io_div,
					L_phi, L_vort, L_div,
					i_simulation_timestamp
			);
		}

		SphereData_Physical L_phi_D_phys(io_phi.sphereDataConfig);
		SphereData_Physical L_vort_D_phys(io_phi.sphereDataConfig);
		SphereData_Physical L_div_D_phys(io_phi.sphereDataConfig);

		sphereSampler.bicubic_scalar(
				L_phi.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				L_phi_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		sphereSampler.bicubic_scalar(
				L_vort.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				L_vort_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		sphereSampler.bicubic_scalar(
				L_div.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				L_div_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		L_phi_D.loadSphereDataPhysical(L_phi_D_phys);
		L_vort_D.loadSphereDataPhysical(L_vort_D_phys);
		L_div_D.loadSphereDataPhysical(L_div_D_phys);
	}
	else
	{
		/*
		 * Method 2) First get variables on departure points, then evaluate L
		 */

		SphereData_Physical io_phi_D_phys(io_phi.sphereDataConfig);
		SphereData_Physical io_vort_D_phys(io_phi.sphereDataConfig);
		SphereData_Physical io_div_D_phys(io_phi.sphereDataConfig);

		sphereSampler.bicubic_scalar(
				io_phi.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				io_phi_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		sphereSampler.bicubic_scalar(
				io_vort.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				io_vort_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		sphereSampler.bicubic_scalar(
				io_div.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				io_div_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		SphereData_Spectral io_phi_D(io_phi.sphereDataConfig);
		SphereData_Spectral io_vort_D(io_phi.sphereDataConfig);
		SphereData_Spectral io_div_D(io_phi.sphereDataConfig);

		io_phi_D.loadSphereDataPhysical(io_phi_D_phys);
		io_vort_D.loadSphereDataPhysical(io_vort_D_phys);
		io_div_D.loadSphereDataPhysical(io_div_D_phys);

		if (coriolis_treatment == CORIOLIS_LINEAR)
		{
			swe_sphere_ts_l_erk->euler_timestep_update(
					io_phi_D, io_vort_D, io_div_D,
					L_phi_D, L_vort_D, L_div_D,
					i_simulation_timestamp
			);
		}
		else
		{
			swe_sphere_ts_lg_erk->euler_timestep_update(
					io_phi_D, io_vort_D, io_div_D,
					L_phi_D, L_vort_D, L_div_D,
					i_simulation_timestamp
			);
		}
	}


	/*
	 * Compute R = X_D + 1/2 dt L_D + dt N*
	 */

	SphereData_Spectral R_phi(io_phi.sphereDataConfig);
	SphereData_Spectral R_vort(io_phi.sphereDataConfig);
	SphereData_Spectral R_div(io_phi.sphereDataConfig);

	R_phi = phi_D + 0.5*i_dt*L_phi_D;
	R_vort = vort_D + 0.5*i_dt*L_vort_D;
	R_div = div_D + 0.5*i_dt*L_div_D;


	/**
	 * Now we care about dt*N*
	 */
	if (include_nonlinear_divergence)
	{
		/**
		 * Compute non-linear term N*
		 *
		 * Extrapolate non-linear terms with the equation
		 * N*(t+0.5dt) = 1/2 [ 2*N(t) - N(t-dt) ]_D + N^t
		 *
		 * Here we assume that
		 * 		N = -Phi' * Div (U)
		 * with non-constant U!
		 * 
		 * The divergence is already 
		 */

		double gh = simVars.sim.gravitation*simVars.sim.h0;



		SphereData_Spectral N_phi_t(io_phi.sphereDataConfig);
		SphereData_Spectral N_phi_prev(io_phi.sphereDataConfig);

		/**
		 * Compute N = - Phi' * Div (U)
		 */
		
		// Compute N(t)
		N_phi_t.loadSphereDataPhysical(
			(-(io_phi-gh)).getSphereDataPhysical()*io_div.getSphereDataPhysical()
		);

		// Compute N(t-dt)
		N_phi_prev.loadSphereDataPhysical(
			(-(phi_prev-gh)).getSphereDataPhysical()*div_prev.getSphereDataPhysical()
		);

		// [ 2*N(t) - N(t-dt) ]_D
		SphereData_Physical N_D_phys(io_phi.sphereDataConfig);
		sphereSampler.bicubic_scalar(
				(2.0*N_phi_t-N_phi_prev).getSphereDataPhysical(),	// field to sample
				pos_lon_d, pos_lat_d,
				N_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		// N_D spectral
		SphereData_Spectral N_D(N_D_phys);

		// Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N^t)
		SphereData_Spectral N_star = 0.5*(N_D + N_phi_t);

		// add nonlinear divergence
		R_phi += i_dt*N_star;
	}

	if (coriolis_treatment == CORIOLIS_NONLINEAR)
	{
		/**
		 * Apply Coriolis Effect in physical VELOCITY space
		 * 
		 * This is require to stay compatible to all other things
		 */
		SphereData_Spectral f_vort_t(io_phi.sphereDataConfig);
		SphereData_Spectral f_vort_t_prev(io_phi.sphereDataConfig);

		SphereData_Spectral f_div_t(io_phi.sphereDataConfig);
		SphereData_Spectral f_div_t_prev(io_phi.sphereDataConfig);

		/**
		 * Calculate f_vort_t and f_div_t caused by Coriolis effect
		 */
		SphereData_Physical tmpg1 = u_lon*fg;
		SphereData_Physical tmpg2 = v_lat*fg;
		if (simVars.misc.sphere_use_robert_functions)
			op.robert_uv_to_vortdiv(tmpg1, tmpg2, f_vort_t, f_div_t);
		else
			op.uv_to_vortdiv(tmpg1, tmpg2, f_vort_t, f_div_t);
		f_vort_t *= -1.0;


		/**
		 * Calculate f_vort_t_prev and f_div_t_prev caused by Coriolis effect
		 */
		tmpg1 = u_lon_prev*fg;
		tmpg2 = v_lat_prev*fg;
		if (simVars.misc.sphere_use_robert_functions)
			op.robert_uv_to_vortdiv(tmpg1, tmpg2, f_vort_t_prev, f_div_t_prev);
		else
			op.uv_to_vortdiv(tmpg1, tmpg2, f_vort_t_prev, f_div_t_prev);
		f_vort_t_prev *= -1.0;


		/**
		 * Handle vorticity
		 */
		// [ 2*N(t) - N(t-dt) ]_D
		SphereData_Physical f_vort_t_D_phys(io_phi.sphereDataConfig);
		sphereSampler.bicubic_scalar(
				(2.0*f_vort_t-f_vort_t_prev).getSphereDataPhysical(),	// field to sample
				pos_lon_d, pos_lat_d,
				f_vort_t_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		// vort_D physical
		SphereData_Spectral f_vort_t_D(f_vort_t_D_phys);

		// Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N^t)
		SphereData_Spectral f_vort_t_star = 0.5*(f_vort_t_D + f_vort_t);

		// add nonlinear vorticity
		R_vort += i_dt*f_vort_t_star;


		/**
		 * Handle divergence
		 */
		// [ 2*N(t) - N(t-dt) ]_D
		SphereData_Physical f_div_t_D_phys(io_phi.sphereDataConfig);
		sphereSampler.bicubic_scalar(
				(2.0*f_div_t-f_div_t_prev).getSphereDataPhysical(),	// field to sample
				pos_lon_d, pos_lat_d,
				f_div_t_D_phys,
				false,
				simVars.disc.semi_lagrangian_interpolation_limiter
			);

		// div_D physical
		SphereData_Spectral f_div_t_D(f_div_t_D_phys);

		// Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N^t)
		SphereData_Spectral f_div_t_star = 0.5*(f_div_t_D + f_div_t);

		// add nonlinear divergence
		R_div += i_dt*f_div_t_star;
	}

	/*
	 * Step 2b) Solve Helmholtz problem
	 * X - 1/2 dt LX = R
	 */
	if (coriolis_treatment == CORIOLIS_LINEAR)
	{
		swe_sphere_ts_l_irk->run_timestep(
				R_phi,
				R_vort,
				R_div,
				0.5*i_dt,
				i_simulation_timestamp
		);
	}
	else
	{
		swe_sphere_ts_lg_irk->run_timestep(
				R_phi,
				R_vort,
				R_div,
				0.5*i_dt,
				i_simulation_timestamp
		);
	}

	/*
	 * Backup previous variables for multi-step SL method
	 */
	phi_prev = io_phi;
	vort_prev = io_vort;
	div_prev = io_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_phi = R_phi;
	io_vort = R_vort;
	io_div = R_div;
}



/*
 * Setup
 */
void SWE_Sphere_TS_ln_settls::setup(
		bool i_include_nonlinear_divergence,
		bool i_original_linear_operator_sl_tretment,
		const std::string &i_coriolis_treatment		// "linear", "nonlinear", "semi-lagrangian"
)
{
	include_nonlinear_divergence = i_include_nonlinear_divergence;
	original_linear_operator_sl_treatment = i_original_linear_operator_sl_tretment;
	
	if (i_coriolis_treatment == "linear")
		coriolis_treatment = CORIOLIS_LINEAR;
	else if (i_coriolis_treatment == "nonlinear")
		coriolis_treatment = CORIOLIS_NONLINEAR;
	else if (i_coriolis_treatment == "semi-lagrangian")
		coriolis_treatment = CORIOLIS_SEMILAGRANGIAN;
	else
		FatalError(std::string("Only coriolis methods 'linear', 'nonlinear', 'semi-lagrangian' supported, but not")+i_coriolis_treatment);

	if (simVars.disc.space_grid_use_c_staggering)
		FatalError("SWE_Sphere_TS_ln_settls: Staggering not supported for l_cn_na_sl_nd_settls");

	// Setup sampler for future interpolations
	sphereSampler.setup(op.sphereDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig);


	pos_lon_a.setup(op.sphereDataConfig->physical_array_data_number_of_elements);
	pos_lat_a.setup(op.sphereDataConfig->physical_array_data_number_of_elements);

	// setup some test sampling points
	// we use 2 arrays - one for each sampling position

	pos_lon_a.update_lambda_array_indices(
		[&](int idx, double &io_data)
		{
			int i = idx % op.sphereDataConfig->physical_num_lon;
			//int j = idx / sphereDataConfig->physical_data_size[0];

			io_data = 2.0*M_PI*(double)i/(double)op.sphereDataConfig->physical_num_lon;
			assert(io_data >= 0);
			assert(io_data < 2.0*M_PI);
		}
	);

	pos_lat_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
		{
			//int i = idx % sphereDataConfig->physical_data_size[0];
			int j = idx / op.sphereDataConfig->physical_num_lon;

			io_data = op.sphereDataConfig->lat[j];

			assert(io_data >= -M_PI*0.5);
			assert(io_data <= M_PI*0.5);
		}
	);

	if (coriolis_treatment == CORIOLIS_LINEAR)
	{
		// initialize with 1st order
		swe_sphere_ts_l_erk = new SWE_Sphere_TS_l_erk(simVars, op);
		swe_sphere_ts_l_erk->setup(1);

		// initialize with 1st order and half time step size
		swe_sphere_ts_l_irk = new SWE_Sphere_TS_l_irk(simVars, op);
		swe_sphere_ts_l_irk->setup(1, 0.5*simVars.timecontrol.current_timestep_size, 0);

		swe_sphere_ts_lg_erk = nullptr;
		swe_sphere_ts_lg_irk = nullptr;
	}
	else
	{
		// initialize with 1st order
		swe_sphere_ts_lg_erk = new SWE_Sphere_TS_lg_erk(simVars, op);
		swe_sphere_ts_lg_erk->setup(1);

		// initialize with 1st order and half time step size
		swe_sphere_ts_lg_irk = new SWE_Sphere_TS_lg_irk(simVars, op);
		swe_sphere_ts_lg_irk->setup(1, 0.5*simVars.timecontrol.current_timestep_size, 0);

		swe_sphere_ts_l_erk = nullptr;
		swe_sphere_ts_l_irk = nullptr;
	}

	fg.setup(op.sphereDataConfig);
	if (simVars.sim.sphere_use_fsphere)
	{
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = simVars.sim.sphere_fsphere_f0;
			}
		);
	}
	else
	{
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu*2.0*simVars.sim.sphere_rotating_coriolis_omega;
			}
		);
	}
}



SWE_Sphere_TS_ln_settls::SWE_Sphere_TS_ln_settls(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		phi_prev(i_op.sphereDataConfig),
		vort_prev(i_op.sphereDataConfig),
		div_prev(i_op.sphereDataConfig),

		pos_lon_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		pos_lat_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.sphereDataConfig->physical_array_data_number_of_elements)
{
}



SWE_Sphere_TS_ln_settls::~SWE_Sphere_TS_ln_settls()
{
	delete swe_sphere_ts_l_erk;
	delete swe_sphere_ts_lg_erk;
	delete swe_sphere_ts_l_irk;
	delete swe_sphere_ts_lg_irk;
}

