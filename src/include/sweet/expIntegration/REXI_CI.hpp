/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_REXI_CI_HPP_
#define SRC_INCLUDE_REXI_CI_HPP_

#include <vector>
#include <sweet/core/SWEETError.hpp>
#include <sweet/expIntegration/ExpFunction.hpp>
#include <sweet/libmath/DQStuff.hpp>


namespace sweet
{

template <
#if SWEET_QUADMATH
	typename T =  __float128,
#else
	typename T =  double,
#endif
	typename TStorage = double	///< storage precision of coefficients - use quad precision per default
>
class REXI_CI
{
	typedef std::complex<T> TComplex;

	std::vector<TComplex> alpha_eval;
	std::vector<TComplex> beta_eval;

	T pi, pi2;
	std::complex<T> I;

public:
	std::vector<std::complex<TStorage>> alpha;
	std::vector<std::complex<TStorage>> beta;

	REXI_CI()
	{
#if SWEET_QUADMATH
		if (sizeof(T) == 16)
		{
			// T == __float128
			pi = DQStuff::fromString<T>("3.14159265358979323846264338327950288");
			pi2 = pi*DQStuff::fromString<T>("2.0");
			I = std::complex<T>(0, 1);
		}
		else
#endif
		if (sizeof(T) == 8)
		{
			// T == double
			pi = M_PI;
			pi2 = M_PI*2.0;
			I = std::complex<T>(0, 1);

		}
		else
		{
			SWEETError("Type T not supported");
		}
	}


	void setup_shifted_circle(
			const std::string &i_function_name,
			int N,
			T max_real_evalue,	/// maximum real value which must be part of the circle
			T inc_imag_evalue	/// imaginary value which should be included in circle
	)
	{
		/*
		 * See CI documentation
		 *
		 * Here we compute a circle which is given by the boundary points.
		 * This avoids large positive EValues leading to numerical cancellation effects
		 */
		if (max_real_evalue > inc_imag_evalue)
			max_real_evalue = inc_imag_evalue;

		T x0 = max_real_evalue;
		T xm = inc_imag_evalue;
		T r = (x0*x0 + xm*xm)/(2.0*x0);
		T center = max_real_evalue-r;

		setup(
				i_function_name,
				N,
				"circle",
				r*2.0,
				r*2.0,
				center
			);
	}


	void setup(
			const std::string &i_function_name,
			int N,
			const std::string &i_primitive_name,
			T i_size_real,
			T i_size_imag,
			T i_mu
	)
	{
		alpha_eval.resize(N);
		beta_eval.resize(N);

		ExpFunction<T> fun;
		fun.setup(i_function_name);

		if (i_primitive_name == "circle")
		{
			if (i_size_real != i_size_imag)
				SWEETError("size along real and imaginary axis must be the same");

			T r = i_size_real*0.5;

			for (int j = 0; j < N; j++)
			{
				/*
				 * TODO:
				 * Figure out why the half-shifted version is a must-do
				 * for the upsN functions (forgot the reason for it)
				 */
				T theta_j;
				if (i_function_name.substr(0, 3) == "ups" || 1)
				{
					// avoid points directly on axes for upsN functions
					theta_j = (T)pi2*((T)j+(T)0.5)/(T)N;
				}
				else
				{
					// allow points directly on axes
					theta_j = (T)pi2*((T)j)/(T)N;
				}

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

				beta_eval[j] = filter_value*fun.eval(gamma_j)*pos;
				alpha_eval[j] = gamma_j;

				beta_eval[j] /= (T)N;
			}
		}
		else if (i_primitive_name == "rectangle")
		{
//			SWEETError("Not yet working :-(");

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
				SWEETError("Total number of points exceeds the total number of requested points");

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
			SWEETError("This primitive is not known");
		}


		alpha.resize(N);
		beta.resize(N);

		for (int j = 0; j < N; j++)
		{
			alpha[j] = {(TStorage)alpha_eval[j].real(), (TStorage)alpha_eval[j].imag()};
			//alpha[j] = -alpha[j];

			beta[j] = {(TStorage)beta_eval[j].real(), (TStorage)beta_eval[j].imag()};
			beta[j] = -beta[j];
		}
	}
};

}

#endif
