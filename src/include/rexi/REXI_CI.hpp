/*
 * REXI_CI.hpp
 *
 *  Created on: 18 Aug 2017
 *      Author: martin
 */

#ifndef SRC_INCLUDE_REXI_REXI_CI_HPP_
#define SRC_INCLUDE_REXI_REXI_CI_HPP_

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
			T r,
			T mu
	)
	{
		alpha_eval.resize(N);
		beta_eval.resize(N);

		REXIFunctions<T> fun;
		fun.setup(i_function_name);

		for (int j = 0; j < N; j++)
		{
			T theta_j = (T)pi2*((T)j+(T)0.5)/(T)N;
			TComplex gamma_j = r*DQStuff::exp(I*theta_j) + mu;

			TComplex k = r*DQStuff::exp(I*theta_j);

			beta_eval[j] = -fun.eval(gamma_j)*k;
			beta_eval[j] /= (T)N;
			alpha_eval[j] = -(k + mu);
		}


		alpha.resize(N);
		beta.resize(N);
		for (int j = 0; j < N; j++)
		{
			alpha[j] = {(double)alpha_eval[j].real(), (double)alpha_eval[j].imag()};
			beta[j] = {(double)beta_eval[j].real(), (double)beta_eval[j].imag()};
		}
	}
};


#endif /* SRC_INCLUDE_REXI_REXI_CI_HPP_ */
