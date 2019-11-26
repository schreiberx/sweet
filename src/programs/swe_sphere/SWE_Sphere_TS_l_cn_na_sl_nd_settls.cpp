/*
 * SWE_Sphere_TS_l_cn_na_sl_nd_settls.cpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2019-10-24: Based on sphere version
 */

#include "SWE_Sphere_TS_l_cn_na_sl_nd_settls.hpp"


/**
 * Solve  SWE with Crank-Nicolson implicit time stepping
 *  (spectral formulation for Helmholtz eq) with semi-Lagrangian
 *   SL-SI-SP
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
		 * First time step
		 */
		phi_prev = io_phi;
		vort_prev = io_vort;
		div_prev = io_div;
	}

	// Out vars
	SphereData_Spectral phi(io_phi.sphereDataConfig);
	SphereData_Spectral vort(io_phi.sphereDataConfig);
	SphereData_Spectral div(io_phi.sphereDataConfig);

	// Departure points and arrival points
	ScalarDataArray posx_d = posx_a;
	ScalarDataArray posy_d = posy_a;

	// Parameters
//	double h_bar = simVars.sim.h0;
//	double g = simVars.sim.gravitation;

	SphereData_Physical u_prev(io_phi.sphereDataConfig);
	SphereData_Physical v_prev(io_phi.sphereDataConfig);
	op.vortdiv_to_uv(vort_prev, div_prev, u_prev, v_prev);

	SphereData_Physical u(io_phi.sphereDataConfig);
	SphereData_Physical v(io_phi.sphereDataConfig);
	op.vortdiv_to_uv(io_vort, io_div, u, v);

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev,	v_prev,
			u, v,
			posx_a,	posy_a,
			i_dt,
			simVars.sim.sphere_radius,
			posx_d,	posy_d,

			simVars.disc.timestepping_order,
			simVars.disc.semi_lagrangian_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold

	);



//	swe_sphere_ts_l_cn.run_timestep(io_phi, io_vort, io_div, i_dt, i_simulation_timestamp);

#if 0
	// Calculate the RHS
	SphereData rhs_u = alpha * io_u + f0 * io_v    - g * op.diff_c_x(io_h);
	SphereData rhs_v =  - f0 * io_u + alpha * io_v - g * op.diff_c_y(io_h);
	SphereData rhs_h = alpha * io_h  - h_bar * div;

	// all the RHS are to be evaluated at the departure points
	rhs_u=sampler2D.bicubic_scalar(rhs_u, posx_d, posy_d, -0.5, -0.5);
	rhs_v=sampler2D.bicubic_scalar(rhs_v, posx_d, posy_d, -0.5, -0.5);
	rhs_h=sampler2D.bicubic_scalar(rhs_h, posx_d, posy_d, -0.5, -0.5);

	//Get data in spectral space
	rhs_u.request_data_spectral();
	rhs_v.request_data_spectral();
	rhs_h.request_data_spectral();

	// Calculate nonlinear term at half timestep and add to RHS of phi eq.

	if (!use_only_linear_divergence) //full nonlinear case
	{
		SphereData hdiv = 2.0 * io_h * div - phi_prev * div_prev;
		SphereData nonlin(io_h.sphereDataConfig);
		if(simVars.misc.use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_SPHERE_SPECTRAL_SPACE
			FatalError("Implicit diffusion only supported with spectral space activated");
#else
			hdiv = op.implicit_diffusion(hdiv, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
#endif
		}
		nonlin = 0.5 * io_h * div + 0.5 * sampler2D.bicubic_scalar(hdiv, posx_d, posy_d, -0.5, -0.5);
		rhs_h = rhs_h - 2.0*nonlin;
		rhs_h.request_data_spectral();
	}


	// Build Helmholtz eq.
	SphereData rhs_div = op.diff_c_x(rhs_u)+op.diff_c_y(rhs_v);
	SphereData rhs_vort = op.diff_c_x(rhs_v)-op.diff_c_y(rhs_u);
	SphereData rhs     = kappa* rhs_h / alpha - h_bar * rhs_div - f0 * h_bar * rhs_vort / alpha;

	// Helmholtz solver
	helmholtz_spectral_solver(kappa, g*h_bar, rhs, h);

	// Update vort and div
	vort = (1/kappa)*
			( alpha *rhs_u + f0 * rhs_v
					- g * alpha * op.diff_c_x(phi)
					- g * f0 * op.diff_c_y(phi))
					;

	div = (1/kappa)*
			( alpha *rhs_v - f0 * rhs_u
					+ g * f0 * op.diff_c_x(phi)
					- g * alpha * op.diff_c_y(phi))
					;


	// Set time (n) as time (n-1)
	phi_prev = io_h;
	vort_prev = io_u;
	div_prev = io_v;

	// output data
	io_h = phi;
	io_u = u;
	io_v = div;
#endif
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_cn_na_sl_nd_settls::setup(
		bool i_use_only_linear_divergence
)
{
	use_only_linear_divergence = i_use_only_linear_divergence;

	if (simVars.disc.space_grid_use_c_staggering)
		FatalError("SWE_Sphere_TS_l_cn_na_sl_nd_settls: Staggering not supported for l_cn_na_sl_nd_settls");

	// Setup sampler for future interpolations
	sampler2D.setup(op.sphereDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig);


	posx_a.setup(op.sphereDataConfig->physical_array_data_number_of_elements);
	posy_a.setup(op.sphereDataConfig->physical_array_data_number_of_elements);

	// setup some test sampling points
	// we use 2 arrays - one for each sampling position

	posx_a.update_lambda_array_indices(
		[&](int idx, double &io_data)
		{
			int i = idx % op.sphereDataConfig->physical_num_lon;
			//int j = idx / sphereDataConfig->physical_data_size[0];

			io_data = 2.0*M_PI*(double)i/(double)op.sphereDataConfig->physical_num_lon;
			assert(io_data >= 0);
			assert(io_data < 2.0*M_PI);
		}
	);
	posy_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
		{
			//int i = idx % sphereDataConfig->physical_data_size[0];
			int j = idx / op.sphereDataConfig->physical_num_lon;

			io_data = op.sphereDataConfig->lat[j];

			assert(io_data >= -M_PI*0.5);
			assert(io_data <= M_PI*0.5);
		}
	);

	swe_sphere_ts_l_cn.setup(0.5, simVars.timecontrol.current_timestep_size, simVars.rexi.use_sphere_extended_modes);
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

		posx_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.sphereDataConfig->physical_array_data_number_of_elements),

		swe_sphere_ts_l_cn(i_simVars, i_op)
{
}



SWE_Sphere_TS_l_cn_na_sl_nd_settls::~SWE_Sphere_TS_l_cn_na_sl_nd_settls()
{
}

