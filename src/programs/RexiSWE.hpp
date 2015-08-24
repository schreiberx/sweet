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
#include <sweet/Complex2DArrayFFT.hpp>
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

	Complex2DArrayFFT op_diff_c_x, op_diff_c_y;
	Complex2DArrayFFT op_diff2_c_x, op_diff2_c_y;

	Complex2DArrayFFT eta0;
	Complex2DArrayFFT u0;
	Complex2DArrayFFT v0;


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
			double i_h,		///< sampling size
			int i_M,		///< number of sampling points
			double i_f,		///< Coriolis force
			std::size_t *i_resolution,		///< resolution of domain
			const double *i_domain_size		///< size of domain
	)
	{
		M = i_M;
		h = i_h;
		tau = i_tau;
		f = i_f;

//		std::cout << "REXI setup: M=" << M << ", h=" << h << ", tau=" << tau << ", f=" << f << std::endl;

		rexi.setup(h, M);

		if (op_diff_c_x.data == nullptr)
		{
			op_diff_c_x.setup(i_resolution);
			op_diff_c_x.op_setup_diff_x(i_domain_size);
			op_diff_c_y.setup(i_resolution);
			op_diff_c_y.op_setup_diff_y(i_domain_size);

			op_diff2_c_x.setup(i_resolution);
			op_diff2_c_x.op_setup_diff2_x(i_domain_size);
			op_diff2_c_y.setup(i_resolution);
			op_diff2_c_y.op_setup_diff2_y(i_domain_size);

			eta0.setup(i_resolution);
			u0.setup(i_resolution);
			v0.setup(i_resolution);
		}
	}



	/**
	 * Solve the REXI of U(t) = exp(L*t)
	 *
	 * See
	 * 		doc/rexi/understanding_rexi.pdf
	 * for further information
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
//			exit(1);
		}

		double eta_bar = i_parameters.setup_h0;
		double g = i_parameters.sim_g;

		eta0.loadRealFromDataArray(io_h);
		u0.loadRealFromDataArray(io_u);
		v0.loadRealFromDataArray(io_v);

		// convert to spectral space
		// scale with inverse of tau
		eta0 = eta0.toSpec()*std::complex<double>(1.0/tau, 0);
		u0 = u0.toSpec()*std::complex<double>(1.0/tau, 0);
		v0 = v0.toSpec()*std::complex<double>(1.0/tau, 0);

		io_h.setAll(0);
		io_u.setAll(0);
		io_v.setAll(0);

#define USE_REXI_HALF	0

		// TODO: compute only half of it
		std::size_t N = rexi.alpha.size();

#if USE_REXI_HALF
		N = N >> 1;
#endif
//		for (std::size_t n = 0; n < N; n++)
//			std::cout << " + alpha " << n << ": " << rexi.alpha[n] << std::endl;

//		for (std::size_t n = 0; n < N; n++)
//			std::cout << " + beta " << n << ": " << rexi.beta_re[n] << std::endl;

		for (std::size_t n = 0; n < N; n++)
		{
#if USE_REXI_HALF
			if (n == N-1)
			{
				// This is the last (middle) pole we process.
				// First, we have to merge the other already computed data
				// with the equation
			}
#endif

			// load alpha (a) and scale by inverse of tau
			complex alpha = rexi.alpha[n]/tau;

			// load kappa (k)
			complex kappa = alpha*alpha + f*f;

			// compute
			// 		kappa - g * eta_bar * D2
			// NOTE!!! We add kappa in cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
			// This is *NOT* straightforward and different to adding a constant for computations.
			// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.
			Complex2DArrayFFT lhs = (-g*eta_bar*(op_diff2_c_x + op_diff2_c_y)).addScalar_Cart(kappa);

			Complex2DArrayFFT rhs =
					kappa/alpha * eta0
					- eta_bar*(op_diff_c_x(u0) + op_diff_c_y(v0))
					- (f*eta_bar/alpha) * (op_diff_c_x(v0) - op_diff_c_y(u0))
				;

			Complex2DArrayFFT eta = rhs.spec_div_element_wise(lhs);

			// v1 = A * (D h(1) - v(0))
			Complex2DArrayFFT uh = u0 - g*op_diff_c_x(eta);
			Complex2DArrayFFT vh = v0 - g*op_diff_c_y(eta);

			Complex2DArrayFFT u1 = 1.0/kappa * (alpha * uh     + f * vh);
			Complex2DArrayFFT v1 = 1.0/kappa * (   -f * uh + alpha * vh);

			// TO CARTESIAN
			Complex2DArrayFFT eta1_cart = eta.toCart();
			Complex2DArrayFFT u1_cart = u1.toCart();
			Complex2DArrayFFT v1_cart = v1.toCart();

			io_h += (eta1_cart*(rexi.beta_re[n])).getRealWithDataArray();
			io_u += (u1_cart*(rexi.beta_re[n])).getRealWithDataArray();
			io_v += (v1_cart*(rexi.beta_re[n])).getRealWithDataArray();
		}
	}


	~RexiSWE()
	{
	}
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
