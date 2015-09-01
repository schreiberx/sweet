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
#include <sweet/SimulationVariables.hpp>
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


	Complex2DArrayFFT op_diff_c_x, op_diff_c_y;
	Complex2DArrayFFT op_diff2_c_x, op_diff2_c_y;

	Complex2DArrayFFT eta0;
	Complex2DArrayFFT u0;
	Complex2DArrayFFT v0;

public:
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
			double i_h,		///< sampling size
			int i_M,		///< number of sampling points
			int i_L,		///< number of sampling points for Gaussian approx
			double i_f,		///< Coriolis force
			std::size_t *i_resolution,		///< resolution of domain
			const double *i_domain_size,		///< size of domain
			bool i_rexi_half = true	///< use half-pole reduction
	)
	{
		M = i_M;
		h = i_h;
		tau = i_tau;
		f = i_f;

//		std::cout << "REXI setup: M=" << M << ", h=" << h << ", tau=" << tau << ", f=" << f << std::endl;

		rexi.setup(h, M, i_L, i_rexi_half);


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
	 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
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
		const SimulationVariables &i_parameters,
		bool i_use_half_reduction = false			///< reduce the REXI computations to its half
	)
	{
		double eta_bar = i_parameters.setup.h0;
		double g = i_parameters.sim.g;

		eta0.loadRealFromDataArray(io_h);
		u0.loadRealFromDataArray(io_u);
		v0.loadRealFromDataArray(io_v);

		// convert to spectral space
		// scale with inverse of tau
		eta0 = eta0.toSpec()*(1.0/tau);
		u0 = u0.toSpec()*(1.0/tau);
		v0 = v0.toSpec()*(1.0/tau);

		io_h.set_all(0);
		io_u.set_all(0);
		io_v.set_all(0);


		// TODO: compute only half of it
		std::size_t N = rexi.alpha.size();

		for (std::size_t n = 0; n < N; n++)
		{
			// load alpha (a) and scale by inverse of tau
			// we flip the sign to account for the -L used in exp(\tau (-L))
			complex alpha = -rexi.alpha[n]/tau;
			complex beta = -rexi.beta_re[n];

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

			Complex2DArrayFFT uh = u0 - g*op_diff_c_x(eta);
			Complex2DArrayFFT vh = v0 - g*op_diff_c_y(eta);

			Complex2DArrayFFT u1 = 1.0/kappa * (alpha * uh     + f * vh);
			Complex2DArrayFFT v1 = 1.0/kappa * (   -f * uh + alpha * vh);

			// TO CARTESIAN
			Complex2DArrayFFT eta1_cart = eta.toCart();
			Complex2DArrayFFT u1_cart = u1.toCart();
			Complex2DArrayFFT v1_cart = v1.toCart();

			io_h += (eta1_cart*beta).getRealWithDataArray();
			io_u += (u1_cart*beta).getRealWithDataArray();
			io_v += (v1_cart*beta).getRealWithDataArray();
		}
	}


	complex conj(complex v)	const
	{
		return complex(v.real(), -v.imag());
	}


	/**
	 * This method computes the analytical solution based on the given initial values.
	 *
	 * See Embid/Madja/1996
	 * Terry/Beth/2014, page 16 and
	 * the maple worksheet in doc/swe_solution_for_L/swe_solution_for_L.mw
	 */
	void run_timestep_direct_solution(
			DataArray<2> &io_h,
			DataArray<2> &io_u,
			DataArray<2> &io_v,

			double i_timestep_size,

			Operators2D &op,
			const SimulationVariables &i_parameters
	)
	{
		double eta_bar = i_parameters.setup.h0;
		double g = i_parameters.sim.g;
		double f = i_parameters.sim.f;
		complex I(0.0,1.0);

		Complex2DArrayFFT i_h(io_h.resolution);
		Complex2DArrayFFT i_u(io_h.resolution);
		Complex2DArrayFFT i_v(io_h.resolution);

		Complex2DArrayFFT o_h(io_h.resolution);
		Complex2DArrayFFT o_u(io_h.resolution);
		Complex2DArrayFFT o_v(io_h.resolution);

		i_h.loadRealFromDataArray(io_h);
		i_h = i_h.toSpec();

		i_u.loadRealFromDataArray(io_u);
		i_u = i_u.toSpec();

		i_v.loadRealFromDataArray(io_v);
		i_v = i_v.toSpec();

		double s0 = i_parameters.sim.domain_size[0];
		double s1 = i_parameters.sim.domain_size[1];

		for (std::size_t ik1 = 0; ik1 < i_h.resolution[1]; ik1++)
		{
			for (std::size_t ik0 = 0; ik0 < i_h.resolution[0]; ik0++)
			{
				if (ik0 == i_h.resolution[0]/2 || ik1 == i_h.resolution[1]/2)
				{
					o_h.set(ik1, ik0, 0, 0);
					o_u.set(ik1, ik0, 0, 0);
					o_v.set(ik1, ik0, 0, 0);
				}

				complex U_hat[3];
				U_hat[0] = i_h.get(ik1, ik0);
				U_hat[1] = i_u.get(ik1, ik0);
				U_hat[2] = i_v.get(ik1, ik0);

				double k0, k1;
				if (ik0 < i_h.resolution[0]/2)
					k0 = (double)ik0;
				else
					k0 = (double)((int)ik0-(int)i_h.resolution[0]);

				if (ik1 < i_h.resolution[1]/2)
					k1 = (double)ik1;
				else
					k1 = (double)((int)ik1-(int)i_h.resolution[1]);

				/*
				 * dimensionful formulation
				 * see doc/swe_solution_for_L
				 */

				double H0 = eta_bar;

				//////////////////////////////////////
				// GENERATED CODE START
				//////////////////////////////////////
				std::complex<double> eigenvalues[3];
				std::complex<double> eigenvectors[3][3];

				if (k0 == 0 && k1 == 0)
				{
//					complex wg = std::sqrt((complex)f*f*s0*s0*s1*s1);

					eigenvalues[0] = 0.0;
					eigenvalues[1] = -1.0*f;
					eigenvalues[2] = f;

					eigenvectors[0][0] = 1.00000000000000;
					eigenvectors[0][1] = 0.0;
					eigenvectors[0][2] = 0.0;
					eigenvectors[1][0] = 0.0;
					eigenvectors[1][1] = -1.0*I;
					eigenvectors[1][2] = 1.00000000000000;
					eigenvectors[2][0] = 0.0;
					eigenvectors[2][1] = I;
					eigenvectors[2][2] = 1.00000000000000;
				}
				else if (k0 == 0)
				{
//					complex wg = std::sqrt((complex)s0*s0*(f*f*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1));

					eigenvalues[0] = 0.0;
					eigenvalues[1] = -1.0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
					eigenvalues[2] = -1.0*I*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k1*k1 - 1.0*f*f*s1*s1);

					eigenvectors[0][0] = (1.0/2.0)*I*1.0/M_PI*f*1.0/g*1.0/k1*s1;
					eigenvectors[0][1] = 1.00000000000000;
					eigenvectors[0][2] = 0.0;
					eigenvectors[1][0] = -2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
					eigenvectors[1][1] = -1.0*I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
					eigenvectors[1][2] = 1.00000000000000;
					eigenvectors[2][0] = 2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
					eigenvectors[2][1] = I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
					eigenvectors[2][2] = 1.00000000000000;
				}
				else if (k1 == 0)
				{
//					complex wg = std::sqrt((complex)s1*s1*(f*f*s0*s0 + 4.0*M_PI*M_PI*g*g*k0*k0));

					eigenvalues[0] = 0.0;
					eigenvalues[1] = -1.0*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
					eigenvalues[2] = -1.0*I*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);

					eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
					eigenvectors[0][1] = 0.0;
					eigenvectors[0][2] = 1.00000000000000;
					eigenvectors[1][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
					eigenvectors[1][1] = -1.0*I*1.0/f*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
					eigenvectors[1][2] = 1.00000000000000;
					eigenvectors[2][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
					eigenvectors[2][1] = 1.0/f*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);
					eigenvectors[2][2] = 1.00000000000000;
				}
				else
				{
//					complex K2 = M_PI*M_PI*k0*k0 + M_PI*M_PI*k1*k1;
					complex w = std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);

//					complex wg = std::sqrt((complex)f*f*s0*s0*s1*s1 + 4.0*M_PI*M_PI*g*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1*s0*s0);
					eigenvalues[0] = 0.0;
					eigenvalues[1] = -1.0*1.0/s0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);
					eigenvalues[2] = -1.0*I*1.0/s0*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 - 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 - 1.0*f*f*s0*s0*s1*s1);

					eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
					eigenvectors[0][1] = -1.0*1.0/k0*k1*s0*1.0/s1;
					eigenvectors[0][2] = 1.00000000000000;
					eigenvectors[1][0] = 2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 + f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
					eigenvectors[1][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 - 1.0*I*f*w);
					eigenvectors[1][2] = 1.00000000000000;
					eigenvectors[2][0] = -2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 - 1.0*f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
					eigenvectors[2][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 + I*f*w);
					eigenvectors[2][2] = 1.00000000000000;
				}




				//////////////////////////////////////
				// GENERATED CODE END
				//////////////////////////////////////


				if (f == 0)
				{
					/*
					 * override if f == 0, see ./sympy_L_spec_decomposition.py executed with LNr=4
					 */
					if (k0 != 0 || k1 != 0)
					{
						double K2 = K2;

						eigenvalues[0] = 0.0;
						eigenvalues[1] = -2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);
						eigenvalues[2] = 2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);

						eigenvectors[0][0] = 0.0;
						eigenvectors[0][1] = -1.0*k1/sqrt(k0*k0 + k1*k1);
						eigenvectors[0][2] = k0/sqrt(k0*k0 + k1*k1);
						eigenvectors[1][0] = -1.0*sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
						eigenvectors[1][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
						eigenvectors[1][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
						eigenvectors[2][0] = sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
						eigenvectors[2][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
						eigenvectors[2][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					}
					else
					{

						eigenvalues[0] = 0.0;
						eigenvalues[1] = 0.0;
						eigenvalues[2] = 0.0;

						eigenvectors[0][0] = 1.00000000000000;
						eigenvectors[0][1] = 0.0;
						eigenvectors[0][2] = 0.0;
						eigenvectors[1][0] = 0.0;
						eigenvectors[1][1] = 1.00000000000000;
						eigenvectors[1][2] = 0.0;
						eigenvectors[2][0] = 0.0;
						eigenvectors[2][1] = 0.0;
						eigenvectors[2][2] = 1.00000000000000;
					}
				}


				/*
				 * Compute inverse of Eigenvectors.
				 * This generalizes to the case that the Eigenvectors are not orthonormal.
				 */
				std::complex<double> eigenvectors_inv[3][3];

				eigenvectors_inv[0][0] =  (eigenvectors[1][1]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][1]);
				eigenvectors_inv[0][1] = -(eigenvectors[0][1]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][1]);
				eigenvectors_inv[0][2] =  (eigenvectors[0][1]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][1]);

				eigenvectors_inv[1][0] = -(eigenvectors[1][0]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][0]);
				eigenvectors_inv[1][1] =  (eigenvectors[0][0]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][0]);
				eigenvectors_inv[1][2] = -(eigenvectors[0][0]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][0]);

				eigenvectors_inv[2][0] =  (eigenvectors[1][0]*eigenvectors[2][1] - eigenvectors[1][1]*eigenvectors[2][0]);
				eigenvectors_inv[2][1] = -(eigenvectors[0][0]*eigenvectors[2][1] - eigenvectors[0][1]*eigenvectors[2][0]);
				eigenvectors_inv[2][2] =  (eigenvectors[0][0]*eigenvectors[1][1] - eigenvectors[0][1]*eigenvectors[1][0]);

			    complex s = eigenvectors[0][0]*eigenvectors_inv[0][0] + eigenvectors[0][1]*eigenvectors_inv[1][0] + eigenvectors[0][2]*eigenvectors_inv[2][0];

			    for (int j = 0; j < 3; j++)
				    for (int i = 0; i < 3; i++)
				    	eigenvectors_inv[j][i] /= s;


			    // check
			    for (int j = 0; j < 3; j++)
			    {
				    for (int i = 0; i < 3; i++)
				    {
				    	if (
				    			std::isnan(eigenvectors[j][i].real()) || std::isinf(eigenvectors[j][i].real())	||
				    			std::isnan(eigenvectors[j][i].imag()) || std::isinf(eigenvectors[j][i].imag())
				    	)
				    	{
				    		std::cerr << "Invalid number in Eigenvector " << j << " detected: " << eigenvectors[j][0] << ", " << eigenvectors[j][1] << ", " << eigenvectors[j][2] << std::endl;
				    	}

				    	if (
				    			std::isnan(eigenvectors_inv[j][i].real()) || std::isinf(eigenvectors_inv[j][i].real())	||
				    			std::isnan(eigenvectors_inv[j][i].imag()) || std::isinf(eigenvectors_inv[j][i].imag())
				    	)
				    	{
				    		std::cerr << "Invalid number in inverse of Eigenvector " << j << " detected: " << eigenvectors_inv[j][0] << ", " << eigenvectors_inv[j][1] << ", " << eigenvectors_inv[j][2] << std::endl;
				    	}
				    }
			    }

				/*
				 * Solve based on previously computed data.
				 * Note, that this data can be also precomputed and reused every time.
				 */
				complex UEV0_sp[3];
				for (int k = 0; k < 3; k++)
				{
					UEV0_sp[k] = {0, 0};
					for (int j = 0; j < 3; j++)
						UEV0_sp[k] += eigenvectors_inv[j][k] * U_hat[j];
				}

				complex omega[3];
				omega[0] = std::exp(-I*eigenvalues[0]*i_timestep_size);
				omega[1] = std::exp(-I*eigenvalues[1]*i_timestep_size);
				omega[2] = std::exp(-I*eigenvalues[2]*i_timestep_size);

				complex U_hat_sp[3];
				for (int k = 0; k < 3; k++)
				{
					U_hat_sp[k] = {0, 0};
					for (int j = 0; j < 3; j++)
						U_hat_sp[k] += eigenvectors[j][k] * omega[j] * UEV0_sp[j];
				}

				o_h.set(ik1, ik0, U_hat_sp[0]);
				o_u.set(ik1, ik0, U_hat_sp[1]);
				o_v.set(ik1, ik0, U_hat_sp[2]);
			}
		}

		io_h = o_h.toCart().getRealWithDataArray();
		io_u = o_u.toCart().getRealWithDataArray();
		io_v = o_v.toCart().getRealWithDataArray();
	}

	~RexiSWE()
	{
	}
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
