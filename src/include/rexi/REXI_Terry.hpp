/*
 * REXI.hpp
 *
 *  Created on: 3 Aug 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_INCLUDE_REXI_REXI_HPP_
#define SRC_INCLUDE_REXI_REXI_HPP_



#include <complex>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <rexi/REXI_Terry_ExponentialApproximation.hpp>
#include <rexi/REXI_Terry_FunApproximation.hpp>
#include <rexi/REXI_Terry_GaussianApproximation.hpp>



/**
 * This class provides the abstraction layer for various forms
 * of REXI.
 * (This is now the extended form beyond e^{ix}.)
 *
 * After calling the setup, the rational approximation can be calculated
 * in the following way:
 *
 * \sum_i Re( (rexi.alpha[i]*I + L)^{-1} * rexi.beta_re[i] )
 *
 * where 'rexi' is an instantiation of the REXI class.
 *
 * The different supported phi functions are as follows
 * (see also Terry et al. paper)
 *
 * 	phi0(x) = e^{ix}
 * 	phi1(x) = (e^{ix} - 1)/(ix)
 * 	phi2(x) = (e^{ix} - ix - 1)/(ix)^2
 *
 * All these functions play an important role in the context
 * of exponential integrators.
 */
template <
	typename T = __float128,	///< evaluation accuracy of coefficients
	typename TStorage = double	///< storage precision of coefficients - use quad precision per default
>
class REXI_Terry
{
public:
	typedef std::complex<T> TComplex;
	typedef std::complex<TStorage> TStorageComplex;

//	REXIFunctions<T> rexiFunctions;

public:
	std::vector<TStorageComplex> alpha;
	std::vector<TStorageComplex> beta_re;


	std::vector<TComplex> alpha_eval;
	std::vector<TComplex> beta_re_eval;



public:
	REXI_Terry()
	{
	}



public:
	REXI_Terry(
			const std::string &i_function_id,
			TStorage i_h,	///< sampling width
			int i_M,		///< approximation area
			int i_L = 0,	///< L, see Gaussian approximation
			bool i_reduce_to_half = true,
			bool i_normalization = true
	)
	{
		setup(i_function_id, i_h, i_M, i_L, i_reduce_to_half, i_normalization);
	}



	static
	std::complex<T> conj(
			const std::complex<T> &i
	)
	{
		return std::complex<T>(i.real(), -i.imag());
	}



public:
	void setup(
		const std::string &i_function_name,
		T i_h,		///< sampling width
		int i_M,				///< approximation area
		int i_L = 0,			///< L value for Gaussian approximation, use 0 for autodetection
		bool i_reduce_to_half = true,	///< reduce the number of poles to half
		bool i_normalization = false
	)
	{
		REXI_Terry_GaussianApproximation<T> ga(i_L);

		int L = ga.L;
		int N = i_M+ga.L;
		int M = i_M;

		alpha_eval.resize(2*N+1);
		beta_re_eval.resize(2*N+1);

		/// temporary storage vector for generalization
		/// over phi functions
		std::vector<TComplex> b;

		if (i_h < 0)
			FatalError("Specify the sampling distance width of REXI (parameter --rexi-h)");

		if (i_function_name == "phi0")
		{
#if 1
			REXI_Terry_ExponentialApproximation<T> exp_approx(i_h, i_M);
			b = exp_approx.b;
#else
			REXI_Terry_FunApproximation<T,T> phia(i_function_name, i_h, i_M);
			b = phia.b;
#endif


#if 0
			typedef double T;

			REXIFunctions<T> rexiFunctions;
			rexiFunctions.setup(i_function_name);

			double d = i_h*i_M;
			for (double x = -d; x <= d+1e-10; x += 1.0)
			{
				std::complex<double> approx = exp_approx.approx(x);
				std::complex<T> tmp = rexiFunctions.eval_returnComplex(std::complex<double>(0, x));

				std::complex<double> analyt(tmp.real(), tmp.imag());
				std::cout << x << ": " << approx << "	" << analyt << std::endl;
			}
#endif
		}
		else
		{
			REXI_Terry_FunApproximation<T> phia(i_function_name, i_h, i_M);
			b = phia.b;

#if 0
			typedef double T;
			//typedef __float128 T;

			REXIFunctions<T> rexiFunctions;
			rexiFunctions.setup(i_function_name);

			double d = i_h*i_M;
			for (double x = -d; x <= d+1e-10; x += 1.0)
			{
				std::complex<double> approx = phia.approx(x);
				std::complex<T> tmp = rexiFunctions.eval_returnComplex(std::complex<double>(0, x));

				std::complex<double> analyt(tmp.real(), tmp.imag());
				std::cout << x << ": " << approx << "	" << analyt << std::endl;
			}
#endif

		}

//		rexiFunctions.setup(i_function_name);


#if 1
		for (int n = 0; n < 2*N+1; n++)
		{
			alpha_eval[n] = {0,0};
			beta_re_eval[n] = {0,0};
		}

		TComplex cmu = ga.mu;
		for (int l = -L; l < L+1; l++)
		{
			for (int m = -M; m < M+1; m++)
			{
				int n = l+m;
				alpha_eval[n+N] = i_h*(cmu + TComplex(0, n));

				beta_re_eval[n+N] += b[m+M].real()*i_h*ga.a[l+L];
			}
		}

#else
		for (int n = -N; n < N+1; n++)
		{
			alpha_eval[n+N] = i_h*(ga.mu + TComplex(0, n));

			int L1 = std::max(-L, n-M);
			int L2 = std::min(L, n+M);

			beta_re_eval[n+N] = 0;
			for (int k = L1; k < L2; k++)
			{
				assert(k+L >= 0);
				assert(k+L < 2*L+1);
				assert(n-k+M >= 0);
				assert(n-k+M < 2*M+1);

				beta_re_eval[n+N] += ga.a[k+L]*ea.b[n-k+M].real();
			}

			beta_re_eval[n+N] *= i_h;
		}
#endif


		if (i_reduce_to_half)
		{
#if 0
			/**
			 * reduce the computational amount to its half,
			 * see understanding REXI in the documentation folder
			 */
			alpha_eval.resize(N+1);
			beta_re_eval.resize(N+1);

			// N+1 contains the pole and we don't rescale this one by 2 but all the other ones
			for (int i = 0; i < N; i++)
				beta_re_eval[i] *= 2.0;

#else
			/*
			 * This is slightly more accurate
			 */
			for (int i = 0; i < N; i++)
			{
//				alpha_eval[i] = (alpha_eval[i] + alpha_eval[N*2-i])*complexEvaluation(0.5);
				beta_re_eval[i] += conj(beta_re_eval[N*2-i]);
			}
			alpha_eval.resize(N+1);
			beta_re_eval.resize(N+1);
#endif
		}


		if (i_normalization)
		{
			if (i_function_name == "phi0")
			{
				{
					TStorageComplex sum = 0;

					std::cout << "Normalize for 0-dispersion modes:" << std::endl;

					for (std::size_t n = 0; n < alpha_eval.size(); n++)
					{
						TStorageComplex b(beta_re_eval[n].real(), beta_re_eval[n].imag());
						TStorageComplex a(alpha_eval[n].real(), alpha_eval[n].imag());
						TStorageComplex val = b/a;
						sum += val;
					}

					T normalization = sum.real();
					std::cout << "REXI sum for geostrophic modes with double precision: " << (double)((TStorage)normalization) << std::endl;
					std::cout << "REXI Error with coefficients used with double precision: " << (double)((TStorage)1.0-normalization) << std::endl;
				}

				/*
				 * Only available for e^{ix}
				 */
				/*
				 * Apply normalization to beta coefficients
				 * This improves problems of over/undershooting for geostrophic modes
				 */
				{
					TComplex sum = 0;
					for (std::size_t n = 0; n < alpha_eval.size(); n++)
						sum += beta_re_eval[n]/alpha_eval[n];

					T normalization = sum.real();
					std::cout << "REXI Error: " << (double)(TStorage)((T)1.0-normalization) << std::endl;
					std::cout << "Using REXI normalization: " << (double)((TStorage)1.0/normalization) << std::endl;

					for (std::size_t n = 0; n < beta_re_eval.size(); n++)
						beta_re_eval[n] /= normalization;
				}
			}
		}

		alpha.resize(alpha_eval.size());
		beta_re.resize(beta_re_eval.size());


		for (std::size_t n = 0; n < alpha.size(); n++)
		{
			alpha[n] = alpha_eval[n];
			beta_re[n] = beta_re_eval[n];
		}
	}


	void output()
	{
		int N = alpha.size();
		std::cout << "N: " << N << std::endl;

//		std::cout << "Alpha:" << std::endl;
		for (int i = 0; i < N; i++)
			std::cout << "alpha[" << i << "] = " << alpha[i] << std::endl;

//		std::cout << "Beta:" << std::endl;
		for (int i = 0; i < N; i++)
			std::cout << "beta_re[" << i << "] = " << beta_re[i] << std::endl;
	}


#if 0
	/**
	 * \return \f$ cos(x) + i*sin(x) \f$
	 */
	TStorageComplex eval_returnComplex(
			TStorage i_x	///< sampling position
	)
	{
		TComplex ix = {0, (T)i_x};
		TComplex retval = rexiFunctions.eval(ix);
		TStorageComplex retdata = {(TStorage)retval.real(), (TStorage)retval.imag()};
		return retdata;
	}
#endif


	/**
	 * \return \f$ Re(cos(x) + i*sin(x)) = cos(x) \f$
	 */
	std::complex<TStorage> approx_returnComplex(
			TStorage i_x
	)
	{
		std::complex<TStorage> sum = 0;

		std::size_t S = alpha.size();

		for (std::size_t n = 0; n < S; n++)
			sum += (DQStuff::convertComplex<TStorage>(beta_re[n]) / (TStorageComplex(0, i_x) + DQStuff::convertComplex<TStorage>(alpha[n])));

		return sum;
	}


	/**
	 * \return \f$ Re(cos(x) + i*sin(x)) = cos(x) \f$
	 */
	T approx_returnReal(
			T i_x
	)
	{
		T sum = 0;

		std::size_t S = alpha.size();

		for (std::size_t n = 0; n < S; n++)
			sum += (DQStuff::convertComplex<T>(beta_re[n]) / (TComplex(0, i_x) + DQStuff::convertComplex<T>(alpha[n]))).real();

		return sum;
	}



	/**
	 * \return \f$ Im(cos(x) + i*sin(x)) = sin(x) \f$
	 *
	 * we simply use a phase shift of M_PI and use the returnReal variant
	 */
	std::complex<T> approx_returnImag(
			T i_x		///< sampling position
	)
	{
		std::size_t S = alpha.size();

		std::complex<T> sum = 0;
		for (std::size_t n = 0; n < S; n++)
			sum += (DQStuff::convertComplex<T>(beta_re[n]) / (TComplex(0, i_x) + DQStuff::convertComplex<T>(alpha[n])));
		return sum;
	}
};


#endif
