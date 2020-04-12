/*
 * SWE_Sphere_TS_ln_settls.cpp
 *
 *  Created on: 31 Mar 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "SWE_Sphere_TS_ln_sl_exp_settls.hpp"
#include <sweet/sphere/SphereData_DebugContainer.hpp>



void SWE_Sphere_TS_ln_sl_exp_settls::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation*simVars.sim.h0;
	io_phi_pert += gh0;
	run_timestep_nonpert(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
	io_phi_pert -= gh0;
}




void SWE_Sphere_TS_ln_sl_exp_settls::run_timestep_1st_order(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	FatalError("1st order formulation not available!");
}



void SWE_Sphere_TS_ln_sl_exp_settls::run_timestep_2nd_order(
		SphereData_Spectral &io_U_phi,	///< prognostic variables
		SphereData_Spectral &io_U_vort,	///< prognostic variables
		SphereData_Spectral &io_U_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;
	double gh = simVars.sim.gravitation * simVars.sim.h0;
	double dt_inv_radius = i_dt/simVars.sim.sphere_radius;

	if (i_dt <= 0)
		FatalError("SWE_Sphere_TS_ln_settls: Only constant time step size allowed (Please set --dt)");

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 * TODO: Maybe we should replace this with some 1st order SL time step or some subcycling
		 */
		U_phi_prev = io_U_phi;
		U_vort_prev = io_U_vort;
		U_div_prev = io_U_div;
	}

	/*
	 *************************************************************************************************
	 * 1. Step: Semi-Lagrangian
	 * GENERATES: pos_lon_d, pos_lat_d
	 *************************************************************************************************
	 *
	 * Compute Lagrangian trajectories based on SETTLS.
	 * This depends on V(t-\Delta t) and V(t).
	 *
	 * See
	 * Hortal, M. (2002).
	 * The development and testing of a new two-time-level semi-Lagrangian scheme (SETTLS) in the ECMWF forecast model.
	 * Q. J. R. Meteorol. Soc., 2, 1671–1687.
	 *
	 */

	SphereData_Physical U_u_prev(sphereDataConfig);
	SphereData_Physical U_v_prev(sphereDataConfig);
	op.vortdiv_to_uv(U_vort_prev, U_div_prev, U_u_prev, U_v_prev, false);

	SphereData_Physical U_u(sphereDataConfig);
	SphereData_Physical U_v(sphereDataConfig);
	op.vortdiv_to_uv(io_U_vort, io_U_div, U_u, U_v, false);

	// Departure points and arrival points
	ScalarDataArray pos_lon_d = pos_lon_a;
	ScalarDataArray pos_lat_d = pos_lat_a;

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			dt_inv_radius*U_u_prev, dt_inv_radius*U_v_prev,
			dt_inv_radius*U_u, dt_inv_radius*U_v,

			pos_lon_a, pos_lat_a,
			pos_lon_d, pos_lat_d,		// OUTPUT

			simVars.disc.timestepping_order, simVars.disc.semi_lagrangian_max_iterations, simVars.disc.semi_lagrangian_convergence_threshold,
			simVars.disc.semi_lagrangian_approximate_sphere_geometry, simVars.disc.semi_lagrangian_interpolation_limiter
	);

	/*
	 *************************************************************************************************
	 * Step 2) Compute state variables at departure points
	 * GENERATES: phi_D, vort_D, div_D
	 *************************************************************************************************
	 */

	/*
	 * Compute U_D
	 */
	SphereData_Spectral U_phi_D;
	SphereData_Spectral U_vort_D;
	SphereData_Spectral U_div_D;

	if (nonlinear_advection_treatment == NL_ADV_SEMILAGRANGIAN)
	{
		SphereData_Physical U_phi_D_phys(sphereDataConfig);
		sphereSampler.bicubic_scalar(io_U_phi.getSphereDataPhysical(), pos_lon_d, pos_lat_d, U_phi_D_phys, false, simVars.disc.semi_lagrangian_interpolation_limiter);
		U_phi_D = U_phi_D_phys;

		SphereData_Physical U_vort_D_phys(sphereDataConfig);
		sphereSampler.bicubic_scalar(io_U_vort.getSphereDataPhysical(), pos_lon_d, pos_lat_d, U_vort_D_phys, false, simVars.disc.semi_lagrangian_interpolation_limiter);
		U_vort_D = U_vort_D_phys;

		SphereData_Physical U_div_D_phys(sphereDataConfig);
		sphereSampler.bicubic_scalar(io_U_div.getSphereDataPhysical(), pos_lon_d, pos_lat_d, U_div_D_phys, false, simVars.disc.semi_lagrangian_interpolation_limiter);
		U_div_D = U_div_D_phys;
	}
	else
	{
		U_phi_D = io_U_phi;
		U_vort_D = io_U_vort;
		U_div_D = io_U_div;
	}


	/*
	 *************************************************************************************************
	 * Step 3) Compute non-linearities
	 * UPDATES: phi_N_t, vort_N_t, div_N_t
	 *************************************************************************************************
	 */


	if (nonlinear_divergence_treatment == NL_DIV_NONLINEAR || linear_coriolis_treatment == CORIOLIS_NONLINEAR)
	{
		/*
		 * Compute nonlinearities at N^n and N^{n-1}
		 */

		SphereData_Spectral N_phi_n(sphereDataConfig);
		SphereData_Spectral N_phi_n_prev(sphereDataConfig);

		if (nonlinear_divergence_treatment == NL_DIV_NONLINEAR)
		{
			// Compute N(t)
			N_phi_n.loadSphereDataPhysical((-(io_U_phi - gh)).getSphereDataPhysical() * io_U_div.getSphereDataPhysical());

			// Compute N(t-dt)
			N_phi_n_prev.loadSphereDataPhysical((-(U_phi_prev - gh)).getSphereDataPhysical() * U_div_prev.getSphereDataPhysical());
		}
		else
		{
			N_phi_n.spectral_set_zero();
			N_phi_n_prev.spectral_set_zero();
		}


		SphereData_Spectral N_vort_n(sphereDataConfig);
		SphereData_Spectral N_div_n(sphereDataConfig);

		SphereData_Spectral N_vort_n_prev(sphereDataConfig);
		SphereData_Spectral N_div_n_prev(sphereDataConfig);

		if (linear_coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			SphereData_Physical ufg = U_u * op.fg;
			SphereData_Physical vfg = U_v * op.fg;

			op.uv_to_vortdiv(ufg, vfg, N_div_n, N_vort_n, false);
			N_vort_n *= -1.0;

			// u_lon and u_lat already available in physical space
			SphereData_Physical ufg_prev = U_u_prev * op.fg;
			SphereData_Physical vfg_prev = U_v_prev * op.fg;

			SphereData_Spectral f_vort_t_prev(sphereDataConfig);
			SphereData_Spectral f_div_t_prev(sphereDataConfig);
			op.uv_to_vortdiv(ufg_prev, vfg_prev, N_div_n_prev, N_vort_n_prev, false);
			N_vort_n_prev *= -1.0;
		}
		else
		{
			N_vort_n.spectral_set_zero();
			N_div_n.spectral_set_zero();

			N_vort_n_prev.spectral_set_zero();
			N_div_n_prev.spectral_set_zero();
		}


		/*
		 * Compute
		 * E = exp(dtL)N^{n-1}
		 */
		swe_sphere_ts_l_rexi->run_timestep_nonpert(
			N_phi_n_prev, N_vort_n_prev, N_div_n_prev,
			i_dt,
			i_simulation_timestamp
		);


		/*
		 * Interpolate
		 * 	2*N^n - E
		 * at departure points
		 */

		/*
		 * Compute nonlinearities at N^n and N^{n-1}
		 */
		SphereData_Spectral N_phi_D;
		if (nonlinear_divergence_treatment == NL_DIV_NONLINEAR)
		{
			SphereData_Physical N_phi_D_phys(sphereDataConfig);
			sphereSampler.bicubic_scalar(
					(2.0*N_phi_n - N_phi_n_prev).getSphereDataPhysical(),
					pos_lon_d,
					pos_lat_d,
					N_phi_D_phys,
					false,
					simVars.disc.semi_lagrangian_interpolation_limiter
				);
			N_phi_D = N_phi_D_phys;
		}
		else
		{
			N_phi_D.setup(sphereDataConfig);
			N_phi_D.spectral_set_zero();
		}


		SphereData_Spectral N_vort_D;
		SphereData_Spectral N_div_D;
		if (linear_coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			SphereData_Physical N_vort_D_phys(sphereDataConfig);
			sphereSampler.bicubic_scalar(
					(2.0*N_vort_n - N_vort_n_prev).getSphereDataPhysical(),
					pos_lon_d,
					pos_lat_d,
					N_vort_D_phys,
					false,
					simVars.disc.semi_lagrangian_interpolation_limiter
				);
			N_vort_D = N_vort_D_phys;

			SphereData_Physical N_div_D_phys(sphereDataConfig);
			sphereSampler.bicubic_scalar(
					(2.0*N_div_n - N_div_n_prev).getSphereDataPhysical(),
					pos_lon_d,
					pos_lat_d,
					N_div_D_phys,
					false,
					simVars.disc.semi_lagrangian_interpolation_limiter
				);
			N_div_D = N_div_D_phys;

		}
		else
		{
			N_vort_D.setup(sphereDataConfig);
			N_vort_D.spectral_set_zero();

			N_div_D.setup(sphereDataConfig);
			N_div_D.spectral_set_zero();
		}

		U_phi_D += (i_dt*0.5)*(N_phi_D + N_phi_n);
		U_vort_D += (i_dt*0.5)*(N_vort_D + N_vort_n);
		U_div_D += (i_dt*0.5)*(N_div_D + N_div_n);
	}


	/*
	 * Compute
	 * E = exp(dtL)(K)
	 */
	swe_sphere_ts_l_rexi->run_timestep_nonpert(
		U_phi_D, U_vort_D, U_div_D,
		i_dt,
		i_simulation_timestamp
	);


	/*
	 * Backup previous variables for multi-step SL method
	 */
	U_phi_prev = io_U_phi;
	U_vort_prev = io_U_vort;
	U_div_prev = io_U_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_U_phi = U_phi_D;
	io_U_vort = U_vort_D;
	io_U_div = U_div_D;
}



void SWE_Sphere_TS_ln_sl_exp_settls::run_timestep_nonpert(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp)
{
	if (timestepping_order == 1)
	{
		run_timestep_1st_order(io_phi, io_vort, io_div, i_dt, i_simulation_timestamp);
	}
	else if (timestepping_order == 2)
	{
		run_timestep_2nd_order(io_phi, io_vort, io_div, i_dt, i_simulation_timestamp);
	}
	else
	{
		FatalError("Only orders 1/2 supported (ERRID 098nd89eje)");
	}
}


void SWE_Sphere_TS_ln_sl_exp_settls::setup(
		int i_timestepping_order,
		LinearGravityTreatment_enum i_linear_treatment,
		LinearCoriolisTreatment_enum i_coriolis_treatment,
		NLAdvectionTreatment_enum i_nonlinear_advection_treatment,
		NLDivergenceTreatment_enum i_nonlinear_divergence_treatment,
		bool i_original_linear_operator_sl_treatment
)
{
	std::cout << " + SWE_Sphere_TS_ln_sl_exp_settls.setup() called" << std::endl;

	linear_treatment = i_linear_treatment;
	linear_coriolis_treatment = i_coriolis_treatment;
	nonlinear_advection_treatment = i_nonlinear_advection_treatment;
	nonlinear_divergence_treatment = i_nonlinear_divergence_treatment;
	timestepping_order = i_timestepping_order;
	original_linear_operator_sl_treatment = i_original_linear_operator_sl_treatment;

	if (nonlinear_advection_treatment == NL_ADV_IGNORE)
	{
		std::cerr << "Using no treatment of Semi-Lagrangian advection would result in instabilities" << std::endl;
		std::cerr << "Benchmark: geostropic_balance, timestepping method: l_irk_settls, happens after about 3 days with checkerboard pattern" << std::endl;
		FatalError("Doesn't make much sense");
	}

	if (timestepping_order != 2)
		FatalError("Invalid time stepping order, only 2nd order supported");

	// Setup sampler for future interpolations
	sphereSampler.setup(op.sphereDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig, simVars);

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

	swe_sphere_ts_l_rexi = new SWE_Sphere_TS_l_rexi(simVars, op);
	swe_sphere_ts_l_rexi->setup_new(
			simVars.rexi,
			"phi0",
			simVars.timecontrol.current_timestep_size,
			linear_coriolis_treatment == CORIOLIS_LINEAR,
			simVars.sim.sphere_use_fsphere
		);
}

SWE_Sphere_TS_ln_sl_exp_settls::SWE_Sphere_TS_ln_sl_exp_settls(SimulationVariables &i_simVars, SphereOperators_SphereData &i_op) :
		simVars(i_simVars), op(i_op),

		U_phi_prev(i_op.sphereDataConfig),
		U_vort_prev(i_op.sphereDataConfig),
		U_div_prev(i_op.sphereDataConfig),
		pos_lon_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		pos_lat_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		posx_d(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.sphereDataConfig->physical_array_data_number_of_elements)
{
}

SWE_Sphere_TS_ln_sl_exp_settls::~SWE_Sphere_TS_ln_sl_exp_settls()
{
	delete swe_sphere_ts_l_rexi;
}

