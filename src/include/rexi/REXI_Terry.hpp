/*
 * REXI.hpp
 *
 *  Created on: 3 Aug 2015
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
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
 * \sum_i Re( (rexi.alpha[i]*I + L)^{-1} * rexi.beta[i] )
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
#if SWEET_QUADMATH
	typename T = __float128,	///< evaluation accuracy of coefficients
#else
	typename T = double,
#endif
	typename TStorage = double	///< storage precision of coefficients - use quad precision per default
>
class REXI_Terry
{
public:
	typedef std::complex<T> TComplex;
	typedef std::complex<TStorage> TStorageComplex;

public:
	std::vector<TStorageComplex> alpha;
	std::vector<TStorageComplex> beta;

	std::vector<TComplex> alpha_eval;
	std::vector<TComplex> beta_eval;

private:
	std::vector<TComplex> alpha_reim_eval;
	std::vector<TComplex> beta_re_eval;
	std::vector<TComplex> beta_im_eval;



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

		alpha_reim_eval.resize(2*N+1);
		beta_re_eval.resize(2*N+1);
		beta_im_eval.resize(2*N+1);

		/// temporary storage vector for generalization
		/// over phi functions
		std::vector<TComplex> b;

		if (i_h < 0)
			SWEETError("Specify the sampling distance width of REXI (parameter --rexi-h)");

		if (i_function_name == "phi0")
		{
			REXI_Terry_ExponentialApproximation<T> exp_approx(i_h, i_M);
			b = exp_approx.b;
		}
		else
		{
			REXI_Terry_FunApproximation<T> phia(i_function_name, i_h, i_M);
			b = phia.b;
		}

//		rexiFunctions.setup(i_function_name);


#if 1
		for (int n = 0; n < 2*N+1; n++)
		{
			alpha_reim_eval[n] = {0,0};
			beta_re_eval[n] = {0,0};
			beta_im_eval[n] = {0,0};
		}

		TComplex cmu = ga.mu;
		for (int l = -L; l < L+1; l++)
		{
			for (int m = -M; m < M+1; m++)
			{
				int n = l+m;
				alpha_reim_eval[n+N] = i_h*(cmu + TComplex(0, n));

				beta_re_eval[n+N] += b[m+M].real()*i_h*ga.a[l+L];
				beta_im_eval[n+N] += b[m+M].imag()*i_h*ga.a[l+L];
			}
		}

#else
		for (int n = -N; n < N+1; n++)
		{
			alpha_reim_eval[n+N] = i_h*(ga.mu + TComplex(0, n));

			int L1 = std::max(-L, n-M);
			int L2 = std::min(L, n+M);

			beta_re_eval[n+N] = 0;
			beta_im_eval[n+N] = 0;
			for (int k = L1; k < L2; k++)
			{
				assert(k+L >= 0);
				assert(k+L < 2*L+1);
				assert(n-k+M >= 0);
				assert(n-k+M < 2*M+1);

				beta_re_eval[n+N] += ga.a[k+L]*ea.b[n-k+M].real();
				beta_im_eval[n+N] += ga.a[k+L]*ea.b[n-k+M].imag();
			}

			beta_re_eval[n+N] *= i_h;
			beta_im_eval[n+N] *= i_h;
		}
#endif



		int Mk = 2*N+1;
		alpha_eval.resize(2*Mk);
		beta_eval.resize(2*Mk);
		TComplex I(0, 1);
		for (int i = 0; i < Mk; i++)
		{
			beta_eval[i] = (T)0.5*(beta_re_eval[i]+ I*beta_im_eval[i]);
			alpha_eval[i] = alpha_reim_eval[i];

			beta_eval[Mk+i] = (T)0.5*conj(-beta_re_eval[Mk-1-i] + I*beta_im_eval[Mk-1-i]);
			alpha_eval[Mk+i] = -alpha_reim_eval[i];
		}


		if (i_reduce_to_half)
		{
			// size of alpha
			int N = alpha_eval.size();

			// half the size
			int hN = N/2;

			// assure that it's odd
			assert(hN & 1);

			for (int i = 0; i < hN/2+1; i++)
			{
				alpha_eval[hN/2+1+i] = alpha_eval[hN+i];
				beta_eval[hN/2+1+i] = beta_eval[hN+i];
			}

			for (int i = 0; i < hN/2; i++)
			{
				beta_eval[i] *= (T)2.0;
				beta_eval[hN/2+1+i] *= (T)2.0;
			}

			alpha_eval.resize(hN+2);
			beta_eval.resize(hN+2);
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
						TStorageComplex b(beta_eval[n].real(), beta_eval[n].imag());
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
						sum += beta_eval[n]/alpha_eval[n];

					T normalization = sum.real();
					std::cout << "REXI Error: " << (double)(TStorage)((T)1.0-normalization) << std::endl;
					std::cout << "Using REXI normalization: " << (double)((TStorage)1.0/normalization) << std::endl;

					for (std::size_t n = 0; n < beta_eval.size(); n++)
						beta_eval[n] /= normalization;
				}
			}
		}


		alpha.resize(alpha_eval.size());
		beta.resize(beta_eval.size());
		for (std::size_t n = 0; n < alpha.size(); n++)
		{
			alpha[n] = -alpha_eval[n];
			beta[n] = beta_eval[n];
		}


#if 0
		alpha_reim.resize(alpha_reim_eval.size());
		beta_re.resize(beta_re_eval.size());
		beta_im.resize(beta_im_eval.size());
		for (std::size_t n = 0; n < alpha_reim.size(); n++)
		{
			alpha_reim[n] = alpha_reim_eval[n];
			beta_re[n] = beta_re_eval[n];
			beta_im[n] = beta_im_eval[n];
		}
#endif
	}


	void output()
	{
		int N = alpha.size();
		std::cout << "N: " << N << std::endl;

		for (int i = 0; i < N; i++)
			std::cout << "alpha[" << i << "] = " << alpha[i] << std::endl;

		for (int i = 0; i < N; i++)
			std::cout << "beta[" << i << "] = " << beta[i] << std::endl;
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
			TStorage i_x,
			bool i_linear_operator = false		///< set to true if applied to linear operator. Then a system of equations must be solved
	)
	{
		std::complex<TStorage> sum = 0;

		std::size_t S = alpha.size();

		if (!i_linear_operator)
		{
			for (std::size_t n = 0; n < S; n++)
				sum += (DQStuff::convertComplex<TStorage>(beta[n]) / (TStorageComplex(0, i_x) + DQStuff::convertComplex<TStorage>(alpha[n])));
		}
		else
		{
			TStorageComplex U0[2] = {1, 0};

			for (std::size_t n = 0; n < S; n++)
			{
				TStorageComplex L[2][2] = {{0, -i_x}, {i_x, 0}};
				L[0][0] += DQStuff::convertComplex<TStorage>(alpha[n]);
				L[1][1] += DQStuff::convertComplex<TStorage>(alpha[n]);

				TStorageComplex d = 1.0/(L[0][0]*L[1][1] - L[0][1]*L[1][0]);
				TStorageComplex invL[2][2] = {{L[1][1]*d, -L[0][1]*d}, {-L[1][0]*d, L[0][0]*d}};

				TStorageComplex ret[2] = {	invL[0][0]*U0[0] + invL[0][1]*U0[1], invL[1][0]*U0[0] + invL[1][1]*U0[1]};
				ret[0] *= beta[n];
				ret[1] *= beta[n];

				sum += DQStuff::Re(ret[0]) + DQStuff::Re(ret[1])*DQStuff::I((TStorage)1.0);
			}
		}

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
			sum += (DQStuff::convertComplex<T>(beta[n]) / (TComplex(0, i_x) + DQStuff::convertComplex<T>(alpha[n]))).real();

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
			sum += (DQStuff::convertComplex<T>(beta[n]) / (TComplex(0, i_x) + DQStuff::convertComplex<T>(alpha[n])));
		return sum;
	}
};


#endif
