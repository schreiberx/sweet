/*
 * rexi_swe.hpp
 *
 *  Created on: 24 Jul 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <schreiberx@gmail.com>
 */
#include "SWE_Plane_REXI.hpp"

#include <sweet/sweetmath.hpp>


#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>


#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>



/**
 * Solve  SWE with Crank-Nicolson implicit time stepping
 *  (spectral formulation for Helmholtz eq) with semi-Lagrangian
 *   SL-SI-SP
 *
 * U_t = L U(0)
 *
 * MS@PP: Are you sure that this is fully-implicit?
 * Fully-implicit would mean that all terms on RHS are dependent on U(tau)
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
bool SWE_Plane_REXI::run_timestep_cn_sl_ts(
	PlaneData &io_h,  ///< Current fields
	PlaneData &io_u,
	PlaneData &io_v,

	PlaneData &io_h_prev,	///< past fields
	PlaneData &io_u_prev,
	PlaneData &io_v_prev,

	ScalarDataArray &i_posx_a, //Arrival point positions in x and y (this is basically the grid)
	ScalarDataArray &i_posy_a,

	double i_timestep_size,	///< timestep size
	int i_param_nonlinear, ///< degree of nonlinearity (0-linear, 1-full nonlinear, 2-only nonlinear adv)

	const SimulationVariables &i_simVars, ///< Parameters for simulation

	PlaneOperators &op,     ///< Operator class
	PlaneDataSampler &sampler2D, ///< Interpolation class
	PlaneDataSemiLagrangian &semiLagrangian  ///< Semi-Lag class
)
{
#if 0
	// Out vars
	PlaneData h(io_h.planeDataConfig);
	PlaneData u(io_h.planeDataConfig);
	PlaneData v(io_h.planeDataConfig);

	// Departure points and arrival points

	ScalarDataArray posx_d = i_posx_a;
	ScalarDataArray posy_d = i_posy_a;

	// Parameters
	double h_bar = i_simVars.sim.h0;
	double g = i_simVars.sim.gravitation;
	double f0 = i_simVars.sim.f0;
	double dt = i_timestep_size;
	double alpha = 2.0/dt;
	double kappa = alpha*alpha;
	double kappa_bar = kappa;
	kappa += f0*f0;
	kappa_bar -= f0*f0;

	if (i_param_nonlinear > 0)
	{
		Staggering staggering;
		assert(staggering.staggering_type == 'a');

		// Calculate departure points
		semiLagrangian.semi_lag_departure_points_settls(
				io_u_prev,	io_v_prev,
				io_u,		io_v,
				i_posx_a,	i_posy_a,
				dt,
				posx_d,	posy_d,
				staggering
		);

	}

	// Calculate Divergence and vorticity spectrally
	PlaneData div = op.diff_c_x(io_u) + op.diff_c_y(io_v);

	// this could be pre-stored
	PlaneData div_prev = op.diff_c_x(io_u_prev) + op.diff_c_y(io_v_prev);

	// Calculate the RHS
	PlaneData rhs_u = alpha * io_u + f0 * io_v    - g * op.diff_c_x(io_h);
	PlaneData rhs_v =  - f0 * io_u + alpha * io_v - g * op.diff_c_y(io_h);
	PlaneData rhs_h = alpha * io_h  - h_bar * div;

	if (i_param_nonlinear > 0)
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
	if (i_param_nonlinear == 1)
	{
		// Calculate nonlinear term interpolated to departure points
		// h*div is calculate in cartesian space (pseudo-spectrally)
		//div.aliasing_zero_high_modes();
		//div_prev.aliasing_zero_high_modes();
		PlaneData hdiv = 2.0 * io_h * div - io_h_prev * div_prev;
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
	helmholtz_spectral_solver(kappa, g*h_bar, rhs, h, op);

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


	//Set time (n) as time (n-1)
	io_h_prev = io_h;
	io_u_prev = io_u;
	io_v_prev = io_v;

	//output data
	io_h = h;
	io_u = u;
	io_v = v;
#endif
	return true;
}


