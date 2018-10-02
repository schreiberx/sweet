/*
 * REXI_CI.hpp
 *
 *  Created on: 18 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_REXI_CI_HPP_
#define SRC_INCLUDE_REXI_CI_HPP_

#include <sweet/FatalError.hpp>
#include <libmath/DQStuff.hpp>
#include <rexi/REXIFunctions.hpp>
#include <vector>


template <
#if SWEET_QUADMATH
	typename T =  __float128
#else
	typename T =  double
#endif
>
class REXI_CI
{
	typedef std::complex<T> TComplex;

	std::vector<TComplex> alpha_eval;
	std::vector<TComplex> beta_eval;

	const T pi = DQStuff::fromString<T>("3.14159265358979323846264338327950288");
	const T pi2 = pi*DQStuff::fromString<T>("2.0");
	const T pi4 = pi*DQStuff::fromString<T>("4.0");
	const std::complex<T> I = std::complex<T>(0, 1);

public:
	std::vector<std::complex<double>> alpha;
	std::vector<std::complex<double>> beta;

	REXI_CI()
	{
	}


	void setup_shifted_circle(
			const std::string &i_function_name,
			int N,
			T max_real_evalue,
			T max_imag_evalue,
			T i_ci_gaussian_filter_scale_a,
			T i_ci_gaussian_filter_dt_norm,
			T i_ci_gaussian_filter_exp_N,
			T i_ci_gaussian_filter_dt
	)
	{
		/*
		 * See CI documentation
		 *
		 * Here we compute a circle which is given by the boundary points.
		 * This avoids large positive EValues leading to numerical cancellation effects
		 */
		T x0 = max_real_evalue;
		T xm = max_imag_evalue;
		T r = (x0*x0 + xm*xm)/(2.0*x0);
		T center = max_real_evalue-r;

		setup(
				i_function_name,
				N,
				"circle",
				r*2.0,
				r*2.0,
				center,
				i_ci_gaussian_filter_scale_a,
				i_ci_gaussian_filter_dt_norm,
				i_ci_gaussian_filter_exp_N,
				i_ci_gaussian_filter_dt
			);
	}


	void setup(
			const std::string &i_function_name,
			int N,
			const std::string &i_primitive_name,
			T i_size_real,
			T i_size_imag,
			T i_mu,
			T i_ci_gaussian_filter_scale_a,
			T i_ci_gaussian_filter_dt_norm,
			T i_ci_gaussian_filter_exp_N,
			T i_ci_gaussian_filter_dt
	)
	{
#if !SWEET_QUADMATH
		FatalError("Don't use this without quad precision support to generate the coefficients!");
#endif

		alpha_eval.resize(N);
		beta_eval.resize(N);

		REXIFunctions<T> fun(i_function_name);

		if (i_primitive_name == "circle")
		{
			if (i_size_real != i_size_imag)
				FatalError("size along real and imaginary axis must be the same");

			T r = i_size_real*0.5;

			for (int j = 0; j < N; j++)
			{
//				T theta_j = (T)pi2*((T)j+(T)0.5)/(T)N;

				// avoid points directoy on axes
				//T theta_j = (T)pi2*((T)j+(T)0.5)/(T)N;

				// allow points directoy on axes
				T theta_j = (T)pi2*((T)j)/(T)N;

				// sampling position of support point
				TComplex pos = r*DQStuff::exp(I*theta_j);

				// shifted position
				TComplex gamma_j = pos + i_mu;

				/*
				 * Handle filter, see
				 * doc/rexi/filter
				 *
				 * We use a Gaussian bump to filter out the fast modes
				 */
				T filter_value = 1.0;
				if (i_ci_gaussian_filter_dt_norm != 0 && i_ci_gaussian_filter_scale_a != 0)
				{
//					T dist = DQStuff::sqrt(gamma_j.real()*gamma_j.real() + gamma_j.imag()*gamma_j.imag());
					T dist = DQStuff::abs(gamma_j.imag());

//					std::cout << (double)gamma_j.real() << std::endl;
//					std::cout << (double)gamma_j.imag() << std::endl;

					filter_value = DQStuff::exp(
										-DQStuff::pow(
												i_ci_gaussian_filter_dt_norm/(i_ci_gaussian_filter_scale_a*i_ci_gaussian_filter_dt)*dist,
												i_ci_gaussian_filter_exp_N
										)
									);
					std::cout << (double)dist << "\t" << (double)filter_value << std::endl;
				}

				beta_eval[j] = filter_value*fun.eval(gamma_j)*pos;
				alpha_eval[j] = gamma_j;

				beta_eval[j] /= (T)N;
			}
		}
		else if (i_primitive_name == "rectangle")
		{
//			FatalError("Not yet working :-(");

			T SRe = i_size_real;
			T SIm = i_size_imag;

			int NRe = 0.5 * SRe * N / (SRe + SIm);
			int NIm = 0.5 * SIm * N / (SRe + SIm);

#if 0
			/*
			 * Make N even numbered
			 */
			if (NRe & 1)
				NRe--;

			if (NIm & 1)
				NIm--;
#endif

			std::cout << "Size: " << (double)i_size_real << " x " << (double)i_size_imag << std::endl;
			std::cout << "Points: " << NRe << " x " << NIm << std::endl;
			int total_N = (2*NRe + 2*NIm);
			std::cout << "Points (total): " << total_N << std::endl;

			if (total_N > N)
				FatalError("Total number of points exceeds the total number of requested points");

			// Sampling distances
			T DRe = SRe/NRe;
			T DIm = SIm/NIm;

			// Sampling points along each rectangle boundary
			auto PRe = [&](int n) -> T { return (T)0.5*(DRe-SRe) + n*DRe; };
			auto PIm = [&](int n) -> T { return (T)0.5*(DIm-SIm) + n*DIm; };

			int j = 0;

			// BOTTOM
			for (int n = 0; n < NRe; n++)
			{
				TComplex z = -I*(T)0.5*SIm + PRe(n) + i_mu;

				beta_eval[j] = fun.eval(z);
				beta_eval[j] *= SRe/NRe;
				beta_eval[j] /= I*pi2;

				alpha_eval[j] = z;
				j++;
			}

			// RIGHT
			for (int n = 0; n < NIm; n++)
			{
				TComplex z = (T)0.5*SRe + I*PIm(n) + i_mu;

				beta_eval[j] = fun.eval(z);
				beta_eval[j] *= I;
				beta_eval[j] *= SIm/NIm;
				beta_eval[j] /= I*pi2;

				alpha_eval[j] = z;
				j++;
			}

			// TOP
			for (int n = 0; n < NRe; n++)
			{
				TComplex z = I*(T)0.5*SIm - PRe(n) + i_mu;

				beta_eval[j] = fun.eval(z);
				beta_eval[j] *= -(T)1.0;
				beta_eval[j] *= SRe/NRe;
				beta_eval[j] /= I*pi2;

				alpha_eval[j] = z;
				j++;
			}


			// LEFT
			for (int n = 0; n < NIm; n++)
			{
				TComplex z = -(T)0.5*SRe - I*PIm(n) + i_mu;

				beta_eval[j] = fun.eval(z);
				beta_eval[j] *= -I;
				beta_eval[j] *= SIm/NIm;
				beta_eval[j] /= I*pi2;

				alpha_eval[j] = z;
				j++;
			}

			N = j;
		}
		else
		{
			FatalError("This primitive is not known");
		}


		alpha.resize(N);
		beta.resize(N);

		for (int j = 0; j < N; j++)
		{
			alpha[j] = {(double)alpha_eval[j].real(), (double)alpha_eval[j].imag()};
			alpha[j] = -alpha[j];

			beta[j] = {(double)beta_eval[j].real(), (double)beta_eval[j].imag()};
			beta[j] = -beta[j];
//			std::cout << j << ":\t" << alpha[j] << std::endl;
		}
	}
};


#endif /* SRC_INCLUDE_REXI_CI_HPP_ */
