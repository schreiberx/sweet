/*
 * REXI.hpp
 *
 *  Created on: 3 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_REXI_REXI_HPP_
#define SRC_INCLUDE_REXI_REXI_HPP_


#include <complex>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "ExponentialApproximation.hpp"
#include "GaussianApproximation.hpp"



class REXI
{
	typedef std::complex<double> complex;

	int L;
	int M;
	int N;

public:
	std::vector<complex> alpha;
	std::vector<complex> beta_re;
	std::vector<complex> beta_im;


public:
	REXI()	:
		L(-1),
		M(-1),
		N(-1)
	{
	}

public:
	REXI(
			double i_h,	///< sampling width
			int i_M		///< approximation area
	)
	{
		setup(i_h, i_M);
	}

public:
	void setup(
		double i_h,				///< sampling width
		int i_M,				///< approximation area
		bool i_reduce_to_half = true	///< reduce the number of poles to half
	)
	{
		GaussianApproximation ga;
		ExponentialApproximation ea(i_h, i_M);

		L = ga.L;
		N = i_M+L;
		M = i_M;

		alpha.resize(2*N+1);
		beta_re.resize(2*N+1);
		beta_im.resize(2*N+1);

#if 1
		for (int n = 0; n < 2*N+1; n++)
		{
			alpha[n] = {0,0};
			beta_re[n] = {0,0};
			beta_im[n] = {0,0};
		}

		for (int l = -L; l < L+1; l++)
		{
			for (int m = -M; m < M+1; m++)
			{
				int n = l+m;
				alpha[n+N] = i_h*(ga.mu + complex(0, n));

				beta_re[n+N] += ea.b[m+M].real()*i_h*ga.a[l+L];
				beta_im[n+N] += ea.b[m+M].imag()*i_h*ga.a[l+L];
			}
		}

#else
		for (int n = -N; n < N+1; n++)
		{
			alpha[n+N] = i_h*(ga.mu + complex(0, n));

			int L1 = std::max(-L, n-M);
			int L2 = std::min(L, n+M);

			beta_re[n+N] = 0;
			for (int k = L1; k < L2; k++)
			{
				assert(k+L >= 0);
				assert(k+L < 2*L+1);
				assert(n-k+M >= 0);
				assert(n-k+M < 2*M+1);

				beta_re[n+N] += ga.a[k+L]*ea.b[n-k+M].real();
			}

			beta_re[n+N] *= i_h;
		}
#endif


		if (i_reduce_to_half)
		{
			/**
			 * reduce the computational amount to its half,
			 * see understanding REXI in the documentation folder
			 */
			alpha.resize(N+1);
			beta_re.resize(N+1);
			beta_im.resize(N+1);

			// N+1 contains the pole and we don't rescale this one by 2 but all the other ones
			for (int i = 0; i < N; i++)
			{
				beta_re[i] *= 2.0;
				beta_im[i] *= 2.0;
			}
		}
	}


	/**
	 * \return \f$ cos(x) + i*sin(x) \f$
	 */
	complex eval_e_ix(
			double i_x	///< sampling position
	)
	{
		return std::exp(complex(0,1)*i_x);
	}


	/**
	 * compute the approximated value of e^{ix}
	 */
	complex approx_e_ix(
			double i_x	///< sampling position
	)
	{
		double sum_re = 0;
		double sum_im = 0;

		std::size_t S = alpha.size();

		// Split computation into real part of \f$ cos(x) \f$ and imaginary part \f$ sin(x) \f$
		for (std::size_t n = 0; n < S; n++)
		{
			complex denom = (complex(0, i_x) + alpha[n]);
			sum_re += (beta_re[n] / denom).real();
			sum_im += (beta_im[n] / denom).real();
		}

		return complex(sum_re, sum_im);
	}


	/**
	 * \return \f$ Re(cos(x) + i*sin(x)) = cos(x) \f$
	 */
	double approx_e_ix_returnReal(
			double i_x
	)
	{
		double sum = 0;

		std::size_t S = alpha.size();

		for (std::size_t n = 0; n < S; n++)
			sum += (beta_re[n] / (complex(0, i_x) + alpha[n])).real();

		return sum;
	}

	/**
	 * \return \f$ Im(cos(x) + i*sin(x)) = sin(x) \f$
	 *
	 * we simply use a phase shift of M_PI and use the returnReal variant
	 */
	double approx_e_ix_returnImag(
			double i_x		///< sampling position
	)
	{
#if 1
		return approx_e_ix_returnReal(i_x-M_PIl*0.5);
#else
		double sum = 0;
		for (int n = 0; n < 2*N+1; n++)
			sum += (beta_im[n] / (complex(0, i_x) + alpha[n])).real();
		return sum;
#endif
	}
};


#endif
