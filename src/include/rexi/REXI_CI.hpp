/*
 * REXI_CI.hpp
 *
 *  Created on: 18 Aug 2017
 *      Author: martin
 */

#ifndef SRC_INCLUDE_REXI_CI_HPP_
#define SRC_INCLUDE_REXI_CI_HPP_

#include <sweet/FatalError.hpp>
#include <libmath/DQStuff.hpp>
#include <rexi/REXIFunctions.hpp>
#include <vector>


template <typename T =  __float128>
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


	void setup(
			const std::string &i_function_name,
			int N,
			const std::string &i_primitive_name,
			T size_real,
			T size_imag,
			T mu
	)
	{
		alpha_eval.resize(N);
		beta_eval.resize(N);

		REXIFunctions<T> fun(i_function_name);

		if (i_primitive_name == "circle")
		{
			if (size_real != size_imag)
				FatalError("size along real and imaginary axis must be the same");

			T r = size_real*0.5;

			for (int j = 0; j < N; j++)
			{
				T theta_j = (T)pi2*((T)j+(T)0.5)/(T)N;

				// sampling position of support point
				TComplex pos = r*DQStuff::exp(I*theta_j);

				// shifted position
				TComplex gamma_j = pos + mu;

				beta_eval[j] = fun.eval(gamma_j)*pos;
				alpha_eval[j] = gamma_j;

				beta_eval[j] /= (T)N;
			}
		}
		else if (i_primitive_name == "rectangle")
		{
//			FatalError("Not yet working :-(");

			T SRe = size_real;
			T SIm = size_imag;

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

			std::cout << "Size: " << (double)size_real << " x " << (double)size_imag << std::endl;
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
				TComplex z = -I*(T)0.5*SIm + PRe(n) + mu;

				beta_eval[j] = fun.eval(z);
				beta_eval[j] *= SRe/NRe;
				beta_eval[j] /= I*pi2;

				alpha_eval[j] = z;
				j++;
			}

			// RIGHT
			for (int n = 0; n < NIm; n++)
			{
				TComplex z = (T)0.5*SRe + I*PIm(n) + mu;

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
				TComplex z = I*(T)0.5*SIm - PRe(n) + mu;

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
				TComplex z = -(T)0.5*SRe - I*PIm(n) + mu;

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
