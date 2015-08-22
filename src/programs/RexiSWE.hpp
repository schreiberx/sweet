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

		std::cout << "REXI setup: M=" << M << ", h=" << h << ", tau=" << tau << ", f=" << f << std::endl;

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

		std::cout << "REXI SWE START" << std::endl;
		std::cout << "  io_h : " << io_h.reduce_norm1_quad()/(double)(i_parameters.res[0]*i_parameters.res[1]) << std::endl;
		std::cout << "  io_u : " << io_u.reduce_norm1_quad()/(double)(i_parameters.res[0]*i_parameters.res[1]) << std::endl;
		std::cout << "  io_v : " << io_v.reduce_norm1_quad()/(double)(i_parameters.res[0]*i_parameters.res[1]) << std::endl;

		io_h.setAll(0);
		io_u.setAll(0);
		io_v.setAll(0);

		std::size_t N = rexi.alpha.size();
		for (std::size_t n = 0; n < N; n++)
		{
			// load alpha (a) and scale by inverse of tau
			complex alpha = rexi.alpha[n]/tau;

#if 0
			/*
			 * VALIDATION CHECK FOR SOLVER (based on Pedros suggestion)
			 *
			 * we test for a valid spectral solver for the side constraints
			 * 		g=1, h0=1, u=0, v=0, f=0
			 *
			 * As a result, we assume
			 * 		\eta := exp(i*(k1*x+k2*y))
			 *
			 * To test the D^2 operator, we can compute
			 * 		D^2 \eta = -(k1^2 + k2^2) \eta
			 *
			 * Then, we can evaluate the following expression for \eta:
			 * 		(\alpha^2 - D^2) \eta
			 * 			= (\alpha^2 + (k1^2 + k2^2)) \eta
			 * 			= (\alpha^2 + |k|^2) \eta
			 *
			 * Given the side constraints, the solver for \eta is given by
			 * 		(\alpha^2 - D^2) \eta = \alpha^2 / alpha \eta_0
			 *
			 * Given the above specialization on \eta, we can solve for \eta_0:
			 *		\eta_0 = (\alpha^2 + |k|^2) / \alpha^2 * alpha * \eta
			 *		       = (1.0 + |k|^2 / \alpha^2) * alpha * \eta
			 *		       = (alpha + |k|^2 / \alpha) * \eta
			 */
			{
				complex alpha = rexi.alpha[n];

				assert(g == 1);
				assert(eta_bar == 1);

				Complex2DArrayFFT test_eta0(i_parameters.res);
				Complex2DArrayFFT test_eta(i_parameters.res);
				Complex2DArrayFFT test_diff1(i_parameters.res);
				Complex2DArrayFFT test_diff2(i_parameters.res);

				Complex2DArrayFFT test_step1(i_parameters.res);
				Complex2DArrayFFT test_step2(i_parameters.res);
				Complex2DArrayFFT test_step3(i_parameters.res);

				double k1 = 2.0*M_PIl*1.0;
				double k2 = 2.0*M_PIl*2.0;

//				alpha = {1,0};
				for (std::size_t j = 0; j < i_parameters.res[1]; j++)
				{
					for (std::size_t i = 0; i < i_parameters.res[0]; i++)
					{
						double x = ((double)i+0.5)/(double)i_parameters.res[0];
						double y = ((double)j+0.5)/(double)i_parameters.res[1];

						double exponent = x*k1 + y*k2;
						std::complex<double> eta_val(cos(exponent), sin(exponent));

						// test function
						// \eta := exp(i*(k1*x+k2*y))
						test_eta.set(j, i, eta_val);

						// \eta_0 = (\alpha*\alpha + |k|^2) / (\alpha*\alpha) * \alpha * \eta
						std::complex<double> eta_fac = (alpha*alpha + (k1*k1 + k2*k2))/(alpha*alpha) * alpha;
						test_eta0.set(j, i, eta_val * eta_fac);

						// \eta_x + \eta_y = (i k1 + i k2) \eta
						std::complex<double> test_d1_fac(0, (k1 + k2));
						test_diff1.set(j, i, eta_val * test_d1_fac);

						// D^2 \eta = -(k1^2 + k2^2) \eta
						std::complex<double> test_d2_fac(-(k1*k1 + k2*k2), 0);
						test_diff2.set(j, i, eta_val * test_d2_fac);

						std::complex<double> test_step1_fac(alpha*alpha + (k1*k1 + k2*k2));
//						std::complex<double> test_step1_fac((k1*k1 + k2*k2));
						test_step1.set(j, i, eta_val*test_step1_fac);

						std::complex<double> test_step2_fac(alpha*alpha);
						test_step2.set(j, i, eta_val*test_step2_fac);

						std::complex<double> test_step3_fac(alpha*alpha + (k1*k1 + k2*k2));
						test_step3.set(j, i, eta_val*test_step3_fac);
					}
				}

				{
					// Test SOLVER, step 1
//					Complex2DArrayFFT lhs = ((-1.0*(op_diff2_c_x + op_diff2_c_y)).addScalar_Spec(alpha*alpha))(test_eta.toSpec());
					Complex2DArrayFFT lhs = ((-1.0*(op_diff2_c_x + op_diff2_c_y)))(test_eta.toSpec());
					lhs = lhs + test_eta.toSpec()*(alpha*alpha);
					double error_0 = (lhs.toCart()-test_step1).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 0 step1: " << error_0 << std::endl;
				}


				{
					// Test SOLVER, step 2
					Complex2DArrayFFT lhs = test_eta.toSpec()*(alpha*alpha);
					double error_0 = (lhs.toCart()-test_step2).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 0 step2: " << error_0 << std::endl;
				}

				{
					// Test SOLVER, step 3
					// NOTE!!!! Here we use addScalar_Cart, since \eta has to be scaled by alpha.
					// If we would use addScalar_Spec, only the zero mode would be scaled
					Complex2DArrayFFT lhs = ((-1.0*(op_diff2_c_x + op_diff2_c_y)).addScalar_Cart(alpha*alpha))(test_eta.toSpec());
//					Complex2DArrayFFT lhs = ((-1.0*(op_diff2_c_x + op_diff2_c_y)))(test_eta.toSpec())
//											 + test_eta.toSpec()*(alpha*alpha);

					double error_0 = (lhs.toCart()-test_step1).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 0 step3: " << error_0 << std::endl;
				}


//				{
					// NOTE!!!!
					// Here we use addScalar_Cart, since \eta has to be scaled by alpha.
					// If we would use addScalar_Spec, only the zero mode would be scaled
					Complex2DArrayFFT lhs = (-1.0*(op_diff2_c_x + op_diff2_c_y)).addScalar_Cart(alpha*alpha);
					Complex2DArrayFFT rhs =	test_eta0.toSpec()*alpha;
					Complex2DArrayFFT eta = rhs.spec_div_element_wise(lhs).toCart();

					// Test SOLVER
					double error_1 = (eta - test_eta).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 1: " << error_1 << std::endl;
//				}

				{
					// Test DIFF
					double error_2 = ((op_diff_c_x + op_diff_c_y)(test_eta.toSpec()).toCart()-test_diff1).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 2: " << error_2 << std::endl;
				}

				{
					// Test DIFF2
					double error_3 = ((op_diff2_c_x + op_diff2_c_y)(test_eta.toSpec()).toCart()-test_diff2).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 3: " << error_3 << std::endl;
				}

				{
					// Test inverse of DIFF2
					double error_3b = (test_diff2.toSpec().spec_div_element_wise(op_diff2_c_x + op_diff2_c_y).toCart() - test_eta).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 3b: " << error_3b << std::endl;
				}

				{
					std::complex<double> test_value(0.7, -2.0);
					Complex2DArrayFFT test_a(i_parameters.res);

					// test add scalar spec
					test_a.setAll(test_value);
					double error_4 = (test_a.toSpec().addScalar_Spec(-test_value)).toCart().reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 4: " << error_4 << std::endl;
				}

				{
					std::complex<double> test_value(0.7, -2.0);
					Complex2DArrayFFT test_a(i_parameters.res);

					// test add scalar cart
					test_a.setAll(test_value);
					double error_5 = (test_a.addScalar_Cart(-test_value)).reduce_norm2_quad()/(i_parameters.res[0]*i_parameters.res[1]);
					std::cout << "ERROR 5: " << error_5 << std::endl;
				}

				io_h = test_eta.getRealWithDataArray();
				io_u = op_diff_c_x(test_eta).getRealWithDataArray();
				io_v = eta.getRealWithDataArray();
				continue;
			}
#endif

#if 0
			// NOTE!!!!
			// Here we use addScalar_Cart, since \eta has to be scaled by alpha.
			// If we would use addScalar_Spec, only the zero mode would be scaled
			Complex2DArrayFFT lhs = (-1.0*(op_diff2_c_x + op_diff2_c_y)).addScalar_Cart(alpha*alpha);
			Complex2DArrayFFT rhs =	eta0*alpha;
			Complex2DArrayFFT eta = rhs.spec_div_element_wise(lhs);

			// convert to real space
			// scaling done by 1.0/(res[0]*res[1])
			eta = eta*rexi.beta_re[n];
			Complex2DArrayFFT eta1_cart = eta.toCart();

//			eta1_cart = eta1_cart;//*(i_parameters.res[0]*i_parameters.res[1]);
			io_h += (eta1_cart).getRealWithDataArray();

#else

			// load kappa (k)
			complex kappa = alpha*alpha + f*f;

			// compute
			// 		kappa - g * eta_bar * D2
			Complex2DArrayFFT lhs = (-g*eta_bar*(op_diff2_c_x + op_diff2_c_y)).addScalar_Cart(kappa);

			Complex2DArrayFFT rhs =
					kappa/alpha * eta0
					- eta_bar*(op_diff_c_x(u0) + op_diff_c_y(v0))
					- (f*eta_bar/alpha) * (op_diff_c_x(v0) - op_diff_c_y(u0))
				;

			Complex2DArrayFFT eta = rhs.spec_div_element_wise(lhs);

//			eta1 = eta1*(i_parameters.res[0]*i_parameters.res[1]);
//			eta1 = eta1*(-1.0);


			// CHECKED, this is OK
//			double error = (lhs*eta1-rhs).reduce_norm2_quad()/(double)(i_parameters.res[0]*i_parameters.res[1]);
//			std::cout << error << std::endl;

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
#endif
		}

		std::cout << "REXI SWE END" << std::endl;
		std::cout << "  o_h : " << io_h.reduce_norm1_quad()/(double)(i_parameters.res[0]*i_parameters.res[1]) << std::endl;
		std::cout << "  o_u : " << io_u.reduce_norm1_quad()/(double)(i_parameters.res[0]*i_parameters.res[1]) << std::endl;
		std::cout << "  o_v : " << io_v.reduce_norm1_quad()/(double)(i_parameters.res[0]*i_parameters.res[1]) << std::endl;
	}


	~RexiSWE()
	{
	}
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
