/*
 * SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_lg_rexi_lc_erk_na_sl_nd_settls.hpp"




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
void SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{

	if (i_fixed_dt <= 0)
		FatalError("SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time - i_simulation_timestamp;

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
	PlaneData h(io_h.planeDataConfig);
	PlaneData u(io_h.planeDataConfig);
	PlaneData v(io_h.planeDataConfig);

	// Departure points and arrival points
	ScalarDataArray posx_d = posx_a;
	ScalarDataArray posy_d = posy_a;

	// Parameters
	double h_bar = simVars.sim.h0;
	double g = simVars.sim.gravitation;
	double f0 = simVars.sim.f0;
	double dt = i_fixed_dt;
	double alpha = 2.0/dt;
	double kappa = alpha*alpha;
	double kappa_bar = kappa;
	kappa += f0*f0;
	kappa_bar -= f0*f0;

	if (with_nonlinear > 0)
	{
		Staggering staggering;
		assert(staggering.staggering_type == 'a');

		// Calculate departure points
		semiLagrangian.semi_lag_departure_points_settls(
				u_prev,	v_prev,
				io_u,		io_v,
				posx_a,	posy_a,
				dt,
				posx_d,	posy_d,
				staggering
		);
	}


	// Calculate Divergence and vorticity spectrally
	PlaneData div = op.diff_c_x(io_u) + op.diff_c_y(io_v);

	// this could be pre-stored
	PlaneData div_prev = op.diff_c_x(u_prev) + op.diff_c_y(v_prev);

	// Calculate the RHS
	PlaneData rhs_u = alpha * io_u + f0 * io_v    - g * op.diff_c_x(io_h);
	PlaneData rhs_v =  - f0 * io_u + alpha * io_v - g * op.diff_c_y(io_h);
	PlaneData rhs_h = alpha * io_h  - h_bar * div;

	if (with_nonlinear > 0)
	{
		// all the RHS are to be evaluated at the departure points
		rhs_u=sampler2D.bicubic_scalar(rhs_u, posx_d, posy_d, -0.5, -0.5);
		rhs_v=sampler2D.bicubic_scalar(rhs_v, posx_d, posy_d, -0.5, -0.5);
		rhs_h=sampler2D.bicubic_scalar(rhs_h, posx_d, posy_d, -0.5, -0.5);

		//Get data in spectral space
		rhs_u.request_data_spectral();
		rhs_v.request_data_spectral();
		rhs_h.request_data_spectral();
	}

	// Calculate nonlinear term at half timestep and add to RHS of h eq.
	if (with_nonlinear > 0)
	{
		// Calculate nonlinear term interpolated to departure points
		// h*div is calculate in cartesian space (pseudo-spectrally)
		//div.aliasing_zero_high_modes();
		//div_prev.aliasing_zero_high_modes();
		PlaneData hdiv = 2.0 * io_h * div - h_prev * div_prev;
		//hdiv.aliasing_zero_high_modes();
		//std::cout<<offcent<<std::endl;
		PlaneData nonlin = 0.5 * io_h * div + 0.5 * sampler2D.bicubic_scalar(hdiv, posx_d, posy_d, -0.5, -0.5);
		//add diffusion
		//nonlin.printSpectrumEnergy_y();
		//nonlin.printSpectrumIndex();

		//nonlin.aliasing_zero_high_modes();
		//nonlin.printSpectrumEnergy_y();
		//nonlin.printSpectrumIndex();
		//nonlin=diff(nonlin);
		//nonlin=op.implicit_diffusion(nonlin,i_simVars.sim.viscosity,i_simVars.sim.viscosity_order );
		//nonlin.printSpectrumIndex();
		//nonlin.aliasing_zero_high_modes();
		//nonlin.printSpectrumIndex();

		//std::cout << "blocked: "  << std::endl;
		//nonlin.printSpectrumEnergy();
		//std::cout << "Nonlinear error: " << nonlin.reduce_maxAbs() << std::endl;
		//std::cout << "Div: " << div.reduce_maxAbs() << std::endl;
		//nonlin=0;
		rhs_h = rhs_h - 2.0*nonlin;
		rhs_h.request_data_spectral();	/// why is there a request_data_spectral()?
	}


	// Build Helmholtz eq.
	PlaneData rhs_div = op.diff_c_x(rhs_u)+op.diff_c_y(rhs_v);
	PlaneData rhs_vort = op.diff_c_x(rhs_v)-op.diff_c_y(rhs_u);
	PlaneData rhs     = kappa* rhs_h / alpha - h_bar * rhs_div - f0 * h_bar * rhs_vort / alpha;

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

	o_dt = i_fixed_dt;
}



/*
 * Setup
 */
void SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk::setup(
		REXI_SimulationVariables &i_rexi,

		int i_with_nonlinear
)
{
	ts_l_rexi.setup(i_rexi);

	with_nonlinear = i_with_nonlinear;

	// Setup sampler for future interpolations
	sampler2D.setup(simVars.sim.domain_size, op.planeDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(simVars.sim.domain_size, op.planeDataConfig);


	PlaneData tmp_x(op.planeDataConfig);
	tmp_x.physical_update_lambda_array_indices(
		[&](int i, int j, double &io_data)
		{
			io_data = ((double)i)*simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
		}
	);
	PlaneData tmp_y(op.planeDataConfig);
	tmp_y.physical_update_lambda_array_indices(
		[&](int i, int j, double &io_data)
		{
			io_data = ((double)j)*simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];
		}
	);

	// Initialize arrival points with h position
	ScalarDataArray pos_x = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_x);
	ScalarDataArray pos_y = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_y);

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*simVars.disc.cell_size[0];
	posy_a = pos_y+0.5*simVars.disc.cell_size[1];


}


SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk::SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		ts_l_rexi(i_simVars, i_op),

		h_prev(i_op.planeDataConfig),
		u_prev(i_op.planeDataConfig),
		v_prev(i_op.planeDataConfig),

		posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements)
{
}



SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk::~SWE_Plane_TS_lg_rexi_lc_erk_nt_sl_nd_erk()
{
}

