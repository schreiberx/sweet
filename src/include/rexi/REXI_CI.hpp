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

				beta_eval[j] = -fun.eval(gamma_j)*pos;
				alpha_eval[j] = -gamma_j;

				beta_eval[j] /= (T)N;
			}
		}
		else if (i_primitive_name == "rectangle")
		{
			FatalError("Not yet working :-(");

			int N_real = 0.5 * size_real * N / (size_real + size_imag);
			int N_imag = 0.5 * size_imag * N / (size_real + size_imag);

			/*
			 * Make N even numbered
			 */
			if (N_real & 1)
				N_real--;

			if (N_imag & 1)
				N_imag--;

			std::cout << "Points: " << N_real << " x " << N_imag << std::endl;
			std::cout << "Size: " << (double)size_real << " x " << (double)size_imag << std::endl;

			T d_real = size_real / (N_real-1);
			T d_imag = size_imag / (N_imag-1);

			int j = 0;
			for (int i = 0; i < N_real; i++)
			{
				TComplex pos;
				pos.real(-size_real*0.5 + d_real*i);

				{
					pos.imag(size_imag*0.5);

					TComplex gamma_j = pos + mu;

					beta_eval[j] = -fun.eval(gamma_j)*pos;
					beta_eval[j] /= size_real;
					alpha_eval[j] = -gamma_j;
					j++;
				}
				{
					pos.imag(-size_imag*0.5);

					TComplex gamma_j = pos + mu;

					beta_eval[j] = -fun.eval(gamma_j)*pos;
					beta_eval[j] /= size_real;
					alpha_eval[j] = -gamma_j;
					j++;
				}
			}

			for (int i = 1; i < N_imag-1; i++)
			{
				TComplex pos;
				pos.imag(-size_imag*0.5 + d_imag*i);

				{
					pos.real(-size_real*0.5);

					TComplex gamma_j = pos + mu;

					beta_eval[j] = -fun.eval(gamma_j)*pos;
					beta_eval[j] /= size_imag;
					alpha_eval[j] = -gamma_j;
					j++;
				}

				{
					pos.real(size_real*0.5);

					TComplex gamma_j = pos + mu;

					beta_eval[j] = -fun.eval(gamma_j)*pos;
					beta_eval[j] /= size_imag;
					alpha_eval[j] = -gamma_j;
					j++;
				}
			}

			N = j;
			alpha_eval.resize(N);
			beta_eval.resize(N);

			for (int j = 0; j < N; j++)
			{
				beta_eval[j] /= (T)N;
				beta_eval[j] *= std::complex<T>(0, 1);
			}
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
			beta[j] = {(double)beta_eval[j].real(), (double)beta_eval[j].imag()};
//			std::cout << j << ":\t" << alpha[j] << std::endl;
		}
	}
};


#endif /* SRC_INCLUDE_REXI_CI_HPP_ */
