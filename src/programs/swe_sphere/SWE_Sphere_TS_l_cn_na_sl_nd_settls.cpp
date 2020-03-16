/*
 * SWE_Sphere_TS_l_cn_na_sl_nd_settls.cpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2019-10-24: Partly based on plane version
 */

#include "SWE_Sphere_TS_l_cn_na_sl_nd_settls.hpp"


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
 *
 */

void SWE_Sphere_TS_l_cn_na_sl_nd_settls::run_timestep(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Sphere_TS_l_cn_na_sl_nd_settls: Only constant time step size allowed (Please set --dt)");

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

			simVars.disc.semi_lagrangian_approximate_sphere_geometry
	);

	#if 0
		/**
		 * Use this to debug for valid lat/lon coordinates
		 */
		pos_lon_d.update_lambda_array_indices(
			[](int, double &v)
			{
				if (! ((v >= 0) && (v <= 2.0*M_PI)))
				{
					std::cout << "lon: " << v << std::endl;
					FatalError("LON");
				}
			}
		);


		pos_lat_d.update_lambda_array_indices(
			[](int, double &v)
			{
				if (! ((v >= -0.5*M_PI) && (v <= 0.5*M_PI)))
				{
					std::cout << "lat: " << v << std::endl;
					FatalError("LAT");
				}
			}
		);
	#endif

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
			false
		);

	sphereSampler.bicubic_scalar(
			io_vort.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			vort_D_phys,
			false
		);

	sphereSampler.bicubic_scalar(
			io_div.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			div_D_phys,
			false
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

	if (original_linear_operator_sl_tretment)
	{
		/*
		 * Method 1) First evaluate L, then sample result at departure point
		 */
		SphereData_Spectral L_phi(io_phi.sphereDataConfig);
		SphereData_Spectral L_vort(io_phi.sphereDataConfig);
		SphereData_Spectral L_div(io_phi.sphereDataConfig);

		swe_sphere_ts_l_erk.euler_timestep_update(
				io_phi, io_vort, io_div,
				L_phi, L_vort, L_div,
				i_simulation_timestamp
		);

		SphereData_Physical L_phi_D_phys(io_phi.sphereDataConfig);
		SphereData_Physical L_vort_D_phys(io_phi.sphereDataConfig);
		SphereData_Physical L_div_D_phys(io_phi.sphereDataConfig);

		sphereSampler.bicubic_scalar(
				L_phi.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				L_phi_D_phys,
				false
			);

		sphereSampler.bicubic_scalar(
				L_vort.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				L_vort_D_phys,
				false
			);

		sphereSampler.bicubic_scalar(
				L_div.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				L_div_D_phys,
				false
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
				false
			);

		sphereSampler.bicubic_scalar(
				io_vort.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				io_vort_D_phys,
				false
			);

		sphereSampler.bicubic_scalar(
				io_div.getSphereDataPhysical(),
				pos_lon_d, pos_lat_d,
				io_div_D_phys,
				false
			);

		SphereData_Spectral io_phi_D(io_phi.sphereDataConfig);
		SphereData_Spectral io_vort_D(io_phi.sphereDataConfig);
		SphereData_Spectral io_div_D(io_phi.sphereDataConfig);

		io_phi_D.loadSphereDataPhysical(io_phi_D_phys);
		io_vort_D.loadSphereDataPhysical(io_vort_D_phys);
		io_div_D.loadSphereDataPhysical(io_div_D_phys);

		swe_sphere_ts_l_erk.euler_timestep_update(
				io_phi_D, io_vort_D, io_div_D,
				L_phi_D, L_vort_D, L_div_D,
				i_simulation_timestamp
		);
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
		/*
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


		/**
		 * Compute N = - Phi' * Div (U)
		 */
		
		// Compute N(t)
		SphereData_Spectral N_t(io_phi.sphereDataConfig);
		N_t.loadSphereDataPhysical(
			-(io_phi-gh).getSphereDataPhysical()*io_div.getSphereDataPhysical()
		);

		// Compute N(t-dt)
		SphereData_Spectral N_t_sub_dt(io_phi.sphereDataConfig);
		N_t_sub_dt.loadSphereDataPhysical(
			-(phi_prev-gh).getSphereDataPhysical()*div_prev.getSphereDataPhysical()
		);

		// [ 2*N(t) - N(t-dt) ]_D
		SphereData_Physical N_D_phys(io_phi.sphereDataConfig);
		sphereSampler.bicubic_scalar(
				(N_t*2.0-N_t_sub_dt).getSphereDataPhysical(),	// field to sample
				pos_lon_d, pos_lat_d,
				N_D_phys,
				false	// no velocity sampling
			);

		SphereData_Spectral N_D(io_phi.sphereDataConfig);
		N_D.loadSphereDataPhysical(N_D_phys);

		// Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N^t)
		SphereData_Spectral N_star = 0.5*(N_D + N_t);

		// add nonlinear divergence
		R_phi += i_dt*N_star;
	}


	/*
	 * Step 2b) Solve Helmholtz problem
	 * X - 1/2 dt LX = R
	 */

	swe_sphere_ts_l_irk.run_timestep(
			R_phi,
			R_vort,
			R_div,
			0.5*i_dt,
			i_simulation_timestamp
	);


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
void SWE_Sphere_TS_l_cn_na_sl_nd_settls::setup(
		bool i_include_nonlinear_divergence,
		bool i_original_linear_operator_sl_tretment
)
{
	include_nonlinear_divergence = i_include_nonlinear_divergence;
	original_linear_operator_sl_tretment = i_original_linear_operator_sl_tretment;

	if (simVars.disc.space_grid_use_c_staggering)
		FatalError("SWE_Sphere_TS_l_cn_na_sl_nd_settls: Staggering not supported for l_cn_na_sl_nd_settls");

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

	swe_sphere_ts_l_erk.setup(1);
	swe_sphere_ts_l_irk.setup(1, 0.5*simVars.timecontrol.current_timestep_size, 0);
}



SWE_Sphere_TS_l_cn_na_sl_nd_settls::SWE_Sphere_TS_l_cn_na_sl_nd_settls(
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
		posy_d(i_op.sphereDataConfig->physical_array_data_number_of_elements),

		swe_sphere_ts_l_erk(i_simVars, i_op),
		swe_sphere_ts_l_irk(i_simVars, i_op)
{
}



SWE_Sphere_TS_l_cn_na_sl_nd_settls::~SWE_Sphere_TS_l_cn_na_sl_nd_settls()
{
}

