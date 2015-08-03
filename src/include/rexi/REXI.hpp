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

	std::vector<complex> alpha;
	std::vector<complex> beta;

public:
	REXI(
		double i_h,	///< sampling width
		int i_M		///< approximation area
	)
	{
		GaussianApproximation ga;
		ExponentialApproximation ea(i_h, i_M);


		L = ga.L;
		N = i_M+L;
		M = i_M;

		alpha.resize(2*N+1);
		beta.resize(2*N+1);

#if 1
		for (int n = 0; n < 2*N+1; n++)
		{
			alpha[n] = {0,0};
			beta[n] = {0,0};
		}

		for (int l = -L; l < L+1; l++)
		{
			for (int m = -M; m < M+1; m++)
			{
				int n = l+m;
				alpha[n+N] = i_h*(ga.mu + complex(0, n));

				beta[n+N] += ea.b[m+M].real()*i_h*ga.a[l+L];
			}
		}

#else
		for (int n = -N; n < N+1; n++)
		{
			alpha[n+N] = i_h*(ga.mu + complex(0, n));

			int L1 = std::max(-L, n-M);
			int L2 = std::min(L, n+M);

			beta[n+N] = 0;
			for (int k = L1; k < L2; k++)
			{
				assert(k+L >= 0);
				assert(k+L < 2*L+1);
				assert(n-k+M >= 0);
				assert(n-k+M < 2*M+1);

				beta[n+N] += ga.a[k+L]*ea.b[n-k+M].real();
			}

			beta[n+N] *= i_h;
		}
#endif
	}


	/**
	 * return
	 * cos(i_x) + i*sin(i_x)
	 */
	complex eval_e_ix(
			double i_x
	)
	{
		return std::exp(complex(0,1)*i_x);
	}


	complex approx_e_ix(
			double i_x
	)
	{
		complex sum = 0;

		// Split computation into real part of cos(i_x) and imaginary part sin(i_x)
		for (int n = 0; n < 2*N+1; n++)
		{
			sum.real(sum.real() + (beta[n] / (complex(0, i_x) + alpha[n])).real());
			sum.imag(sum.imag() + (beta[n] / (complex(0, i_x-M_PIl*0.5) + alpha[n])).real());
		}

		return sum;
	}


	/**
	 * return
	 * Re(cos(i_x) + i*sin(i_x)) = cos(i_x)
	 */
	double approx_e_ix_returnReal(
			double i_x
	)
	{
		double sum = 0;

		for (int n = 0; n < 2*N+1; n++)
			sum += (beta[n] / (complex(0, i_x) + alpha[n])).real();

		return sum;
	}

	/**
	 * return
	 * Im(cos(i_x) + i*sin(i_x)) = sin(i_x)
	 *
	 * we simply use a phase shift of M_PI and use the returnReal variant
	 */
	double approx_e_ix_returnImag(
			double i_x
	)
	{
		return approx_e_ix_returnReal(i_x-M_PIl*0.5);
	}
};


#endif /* SRC_INCLUDE_REXI_REXI_HPP_ */
