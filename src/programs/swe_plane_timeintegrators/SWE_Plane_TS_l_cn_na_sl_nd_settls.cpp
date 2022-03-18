/*
 * SWE_Plane_TS_l_cn_na_sl_nd_settls.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_cn_na_sl_nd_settls.hpp"

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
void SWE_Plane_TS_l_cn_na_sl_nd_settls::run_timestep(
		PlaneData_Spectral &io_h,	///< prognostic variables
		PlaneData_Spectral &io_u,	///< prognostic variables
		PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_cn_na_sl_nd_settls: Only constant time step size allowed (Please set --dt)");

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step
		 */
		h_prev = io_h;
		u_prev = io_u;
		v_prev = io_v;
	}

	// Out vars
	PlaneData_Spectral h(io_h.planeDataConfig);
	PlaneData_Spectral u(io_h.planeDataConfig);
	PlaneData_Spectral v(io_h.planeDataConfig);

	// Departure points and arrival points
	ScalarDataArray posx_d = posx_a;
	ScalarDataArray posy_d = posy_a;

	// Parameters
	double h_bar = simVars.sim.h0;
	double g = simVars.sim.gravitation;
	double f0 = simVars.sim.plane_rotating_f0;
	double dt = i_dt;
	double alpha = 2.0/dt;
	double kappa = alpha*alpha;
	double kappa_bar = kappa;
	kappa += f0*f0;
	kappa_bar -= f0*f0;

	Staggering staggering;
	assert(staggering.staggering_type == 'a');

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev.toPhys(),	v_prev.toPhys(),
			io_u.toPhys(),		io_v.toPhys(),
			posx_a,	posy_a,
			dt,
			posx_d,	posy_d,
			simVars.sim.plane_domain_size,
			&staggering,
			simVars.disc.timestepping_order,

			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold

	);


	// Calculate Divergence and vorticity spectrally
	PlaneData_Spectral div = op.diff_c_x(io_u) + op.diff_c_y(io_v);

	// This could be pre-stored
	PlaneData_Spectral div_prev = op.diff_c_x(u_prev) + op.diff_c_y(v_prev);

	/**
	 * Calculate the RHS
	 * 
	 * The terms related to alpha include the current solution.
	 * 
	 * In this implementation, the original formulation is rearranged to
	 * 
	 *    U + 1/2 dt L U + dt N(U)
	 * 
	 *    = 1/2 * dt (2.0/dt U + L U + 2.0 * N(U))
	 */
	PlaneData_Spectral rhs_u = alpha * io_u + f0 * io_v    - g * op.diff_c_x(io_h);
	PlaneData_Spectral rhs_v =  - f0 * io_u + alpha * io_v - g * op.diff_c_y(io_h);
	PlaneData_Spectral rhs_h = alpha * io_h - h_bar * div;

	// All the RHS are to be evaluated at the departure points
	PlaneData_Physical rhs_u_phys = rhs_u.toPhys();
	PlaneData_Physical rhs_v_phys = rhs_v.toPhys();
	PlaneData_Physical rhs_h_phys = rhs_h.toPhys();
	rhs_u = sampler2D.bicubic_scalar(rhs_u_phys, posx_d, posy_d, -0.5, -0.5);
	rhs_v = sampler2D.bicubic_scalar(rhs_v_phys, posx_d, posy_d, -0.5, -0.5);
	rhs_h = sampler2D.bicubic_scalar(rhs_h_phys, posx_d, posy_d, -0.5, -0.5);

///	// Get data in spectral space
///	rhs_u.request_data_spectral();
///	rhs_v.request_data_spectral();
///	rhs_h.request_data_spectral();

	// Calculate nonlinear term at half timestep and add to RHS of h eq.

	if (!use_only_linear_divergence) //full nonlinear case
	{
		// Extrapolation
		PlaneData_Spectral hdiv = 2.0 * io_h * div - h_prev * div_prev;
		PlaneData_Spectral nonlin(io_h.planeDataConfig);
		if(simVars.misc.use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
			SWEETError("Implicit diffusion only supported with spectral space activated");
#else
			// Add diffusion (stabilisation)
			hdiv = op.implicit_diffusion(hdiv, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
#endif
		}
		// Average
	        PlaneData_Physical hdiv_phys = hdiv.toPhys();
		nonlin = 0.5*(io_h*div) + 0.5*sampler2D.bicubic_scalar(hdiv_phys, posx_d, posy_d, -0.5, -0.5);

		// Add to RHS h (TODO (2020-03-16): No clue why there's a -2.0)
		rhs_h = rhs_h - 2.0*nonlin;
///		rhs_h.request_data_spectral();
	}


	// Build Helmholtz eq.
	PlaneData_Spectral rhs_div = op.diff_c_x(rhs_u)+op.diff_c_y(rhs_v);
	PlaneData_Spectral rhs_vort = op.diff_c_x(rhs_v)-op.diff_c_y(rhs_u);
	PlaneData_Spectral rhs     = kappa* rhs_h / alpha - h_bar * rhs_div - f0 * h_bar * rhs_vort / alpha;

	// Helmholtz solver
	helmholtz_spectral_solver(kappa, g*h_bar, rhs, h);

	// Update u and v
	u = (1/kappa)*
			( alpha *rhs_u + f0 * rhs_v
					- g * alpha * op.diff_c_x(h)
					- g * f0 * op.diff_c_y(h))
					;

	v = (1/kappa)*
			( alpha *rhs_v - f0 * rhs_u
					+ g * f0 * op.diff_c_x(h)
					- g * alpha * op.diff_c_y(h))
					;


	// Set time (n) as time (n-1)
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	// output data
	io_h = h;
	io_u = u;
	io_v = v;
}



/*
 * Setup
 */
void SWE_Plane_TS_l_cn_na_sl_nd_settls::setup(
		bool i_use_only_linear_divergence
)
{
	use_only_linear_divergence = i_use_only_linear_divergence;

	if (simVars.disc.space_grid_use_c_staggering)
		SWEETError("SWE_Plane_TS_l_cn_na_sl_nd_settls: Staggering not supported for l_cn_na_sl_nd_settls");

	// Setup sampler for future interpolations
	sampler2D.setup(simVars.sim.plane_domain_size, op.planeDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(simVars.sim.plane_domain_size, op.planeDataConfig);


	PlaneData_Physical tmp_x(op.planeDataConfig);
	tmp_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)i)*simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0];
			},
			false
	);
	PlaneData_Physical tmp_y(op.planeDataConfig);
	tmp_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)j)*simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1];
			},
			false
	);

	//Initialize arrival points with h position
	ScalarDataArray pos_x = Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_x);
	ScalarDataArray pos_y = Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_y);

	double cell_size_x = simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0];
	double cell_size_y = simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1];

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*cell_size_x;
	posy_a = pos_y+0.5*cell_size_y;


}


SWE_Plane_TS_l_cn_na_sl_nd_settls::SWE_Plane_TS_l_cn_na_sl_nd_settls(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		h_prev(i_op.planeDataConfig),
		u_prev(i_op.planeDataConfig),
		v_prev(i_op.planeDataConfig),

		posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements)
{
}



SWE_Plane_TS_l_cn_na_sl_nd_settls::~SWE_Plane_TS_l_cn_na_sl_nd_settls()
{
}

