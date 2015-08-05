/*
 * rexi_swe.hpp
 *
 *  Created on: 24 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_PROGRAMS_REXISWE_HPP_
#define SRC_PROGRAMS_REXISWE_HPP_

#if SWEET_USE_SPECTRAL_SPACE != 1
	#error	"Spectral space required for solvers"
#endif

#include <rexi/REXI.hpp>
#include <sweet/DataArray.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/SimulationParameters.hpp>
#include <complex>
#include <ostream>
#include <fstream>
#include <vector>

/**
 * This class implements the REXI (rational approximation of exponential integrator) solver for the SWE,
 * see High-order time-parallel approximation of evolution operators, T. Haut et. al.
 */
class RexiSWE
{
	typedef std::complex<double> complex;

	double tau;
	double h;
	int M;
	double f;

	REXI rexi;

public:
	RexiSWE()
	{
	}

	/**
	 * setup the REXI
	 */
public:
	void setup(
			double i_tau,	///< time step size
			double i_h,
			int i_M,		///< number of sampling points
			double i_f		///< Coriolis force
	)
	{
		M = i_M;
		h = i_h;
		tau = i_tau;
		f = i_f;

		rexi.setup(h, M);
	}



	/**
	 * Solve the REXI of U(t) = exp(L*t)
	 *
	 * See also
	 * 	apply_rational_func_expL_wave_prop
	 * in Terrys code
	 *
	 * See also doc/rexi/understanding_rexi.pdf
	 */
	void run_timestep(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		Operators2D &op,
		const SimulationParameters &i_parameters
	)
	{
		if (i_parameters.sim_g != 1.0 || i_parameters.setup_h0 != 1.0)
		{
			std::cerr << "Only non-dimensional formulation supported yet, use g=1 and h0=1" << std::endl;
			exit(1);
		}

		DataArray<2> o_h(io_h.resolution);
		DataArray<2> o_u(io_h.resolution);
		DataArray<2> o_v(io_h.resolution);

		o_h.setAll(0);
		o_u.setAll(0);
		o_v.setAll(0);

		Complex2DArrayFFT op_diff_c_x(io_h.resolution), op_diff_c_y(io_h.resolution);
		op_diff_c_x.op_setup_diff_x(i_parameters.sim_domain_size);
		op_diff_c_y.op_setup_diff_y(i_parameters.sim_domain_size);

		Complex2DArrayFFT op_diff2_c_x(io_h.resolution), op_diff2_c_y(io_h.resolution);
		op_diff2_c_x.op_setup_diff2_x(i_parameters.sim_domain_size);
		op_diff2_c_y.op_setup_diff2_y(i_parameters.sim_domain_size);

#define TEST_IMAG_ZERO	1

#if TEST_IMAG_ZERO
		DataArray<2> o_h_im(o_h.resolution);
		DataArray<2> o_u_im(o_u.resolution);
		DataArray<2> o_v_im(o_v.resolution);

		o_h_im.setAll(0);
		o_u_im.setAll(0);
		o_v_im.setAll(0);
#endif

		for (int m = 0; m < M*2+1; m++)
		{
			/*
			 * (D^2 - K) h1 = k/a h0 - f/a D x v0 + D.v0
			 */

			Complex2DArrayFFT h0(o_h.resolution);	h0.setAll(0,0);		h0.loadRealFromDataArray(io_h);
			Complex2DArrayFFT u0(o_u.resolution);	u0.setAll(0,0);		u0.loadRealFromDataArray(io_u);
			Complex2DArrayFFT v0(o_v.resolution);	v0.setAll(0,0);		v0.loadRealFromDataArray(io_v);

			h0 = h0.toSpec();
			u0 = u0.toSpec();
			v0 = v0.toSpec();

			// load alpha (a) and scale by inverse of tau
			complex alpha = rexi.alpha[m]/tau;

			// load kappa (k)
			complex kappa = alpha*alpha + f*f;

			Complex2DArrayFFT rhs =
					kappa/alpha * h0
					-(f/alpha) * (op_diff_c_x(v0) - op_diff_c_y(u0))
					+(op_diff_c_x(u0) + op_diff_c_y(v0))
				;

			Complex2DArrayFFT lhs = (op_diff2_c_x + op_diff2_c_y).subScalar_Spec(kappa);

			Complex2DArrayFFT h1 = rhs.spec_div_element_wise(lhs);

			// v1 = A * (D h(1) - v(0))
			Complex2DArrayFFT uh = op_diff_c_x(h1) - u0;
			Complex2DArrayFFT vh = op_diff_c_y(h1) - v0;

			Complex2DArrayFFT u1 = 1.0/kappa * (alpha * uh     - f * vh);
			Complex2DArrayFFT v1 = 1.0/kappa * (    f * uh + alpha * vh);

			// TO CARTESIAN and scale with TAU
			h1 = h1.toCart()*tau;
			u1 = u1.toCart()*tau;
			v1 = v1.toCart()*tau;

			o_h += (h1*rexi.beta_re[m]).getRealWithDataArray();
			o_u += (u1*rexi.beta_re[m]).getRealWithDataArray();
			o_v += (v1*rexi.beta_re[m]).getRealWithDataArray();

#if TEST_IMAG_ZERO
			o_h_im += (h1*rexi.beta_im[m]).getRealWithDataArray();
			o_u_im += (u1*rexi.beta_im[m]).getRealWithDataArray();
			o_v_im += (v1*rexi.beta_im[m]).getRealWithDataArray();
#endif
		}

		io_h = o_h;
		io_u = o_u;
		io_v = o_v;

		std::cout << "o_h : " << o_h.reduce_norm1_quad() << std::endl;
		std::cout << "o_u : " << o_u.reduce_norm1_quad() << std::endl;
		std::cout << "o_v : " << o_v.reduce_norm1_quad() << std::endl;

#if TEST_IMAG_ZERO
		o_h_im /= tau;
		o_u_im /= tau;
		o_v_im /= tau;

		std::cout << "o_h_im : " << o_h_im.reduce_norm1_quad() << std::endl;
		std::cout << "o_u_im : " << o_u_im.reduce_norm1_quad() << std::endl;
		std::cout << "o_v_im : " << o_v_im.reduce_norm1_quad() << std::endl;
#endif
	}


	~RexiSWE()
	{
	}
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
