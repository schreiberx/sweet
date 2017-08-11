/*
 * PhiApproximation.hpp
 *
 *  Created on: 10 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_INCLUDE_REXI_FUNAPPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_FUNAPPROXIMATION_HPP_

#include <sweet/sweetmath.hpp>
#include <libmath/DQStuff.hpp>
#include <libmath/GaussQuadrature.hpp>
#include "GaussianApproximation.hpp"



/**
 * This class computes an approximation of functions suitable for the
 * computation of exponential integrator formulations including non-linearities
 *
 * \f$
 * 	phi1(x) := (e^{ix}-1)/ix
 * \f$
 *
 * with a combination of Gaussians.
 *
 * For more information, see Terry Haut et al. paper
 * "A high-order time-parallel scheme for solving wave propagation problems via
 * the direct construction of an approximate time-evolution operator"
 *
 * We use a numerical quadrature to compute the following coefficients (see Terry's paper):
 *
 * \f$
 * 		c_m := \int_{-1/(2h)}^{1/2h} h * \exp(-2 \pi i m h \xi) \tilde\phi_N(\xi) / \tilde\psi_h(\xi) d \xi
 * \f$
 * with the Gaussian
 * \f$
 * 		\psi_h(\xi) = 1.0/\sqrt(4.0*\pi) \exp(-x*x/(4*h*h))
 * \f$
 * which is in Fourier space
 * [see https://www.wolframalpha.com/input/?i=Fourier+transform+sqrt(2*pi)*(1%2Fsqrt(4*pi)+*+exp(-x*x%2F(4*h*h))))) ]
 * \f$
 * 		\tilde\psi_h(\xi) = h*\exp(-h*h*\xi*\xi)
 * \f$
 *
 * Combining both equations leads to
 *
 * \f$
 * 		c_m := \int_{-1/(2h)}^{1/2h}  exp(-2 \pi i m h \xi) / (h* \exp(-h*h*\xi*\xi)) * phiN(\xi) d \xi
 * \f$
 */
template <
	typename TREXIComputation = __float128,	///< Accuracy for computation of REXI coefficients
	typename TStorage = double	///< Storage precision
>
class FunApproximation
{
	typedef std::complex<TREXIComputation> complexREXIComputation;
	typedef std::complex<TStorage> complexStorage;

	GaussianApproximation<TREXIComputation, TREXIComputation> ga;


	TREXIComputation h;
	int M;

	const TREXIComputation pi = DQStuff::fromString<TREXIComputation>("3.14159265358979323846264338327950288");
	const TREXIComputation pi2 = pi*DQStuff::fromString<TREXIComputation>("2.0");
	const TREXIComputation pi4 = pi*DQStuff::fromString<TREXIComputation>("4.0");
	const TREXIComputation sqrtpi4 = DQStuff::sqrt(pi4);
	TREXIComputation int_threshold;

public:
	std::vector<complexStorage> b;

	std::complex<TREXIComputation> rexi_int_fun(
			int m,
			TREXIComputation i_xi,
			std::function<std::complex<TREXIComputation>(TREXIComputation)> fun_phi_tilde	///< spectral representation of function to be approximated,
	)
	{
		/*
		 * Note that we had to change the equation from
		 *
		 * DQStuff::expIm(-pi2*(TREXIComputation)m*h*i_xi)/(h*DQStuff::exp(-h*h*i_xi*i_xi)) * fun_phi_tilde(i_xi);
		 *
		 * to
		 *
		 * DQStuff::expIm(pi2*(TREXIComputation)m*h*i_xi)/(h*DQStuff::exp(-h*h*i_xi*i_xi*pi2*pi2)) * fun_phi_tilde(i_xi*pi2) * h;
		 *
		 * to get correct coefficients. This is since FT's are differently defined.
		 */
		return DQStuff::expIm(pi2*(TREXIComputation)m*h*i_xi)/(h*DQStuff::exp(-h*h*i_xi*i_xi*pi2*pi2)) * fun_phi_tilde(i_xi*pi2) * h;
	};

	std::complex<TREXIComputation> computeREXICoefficient(
			int m,
			std::function<std::complex<TREXIComputation>(TREXIComputation)> fun_phi_tilde,	///< spectral representation of function to be approximated,
			TREXIComputation i_int_start,
			TREXIComputation i_int_end
	)
	{
		TREXIComputation real = GaussQuadrature::integrate5_intervals_adaptive_recursive<TREXIComputation>(
						i_int_start,	// start of quadrature
						i_int_end,	// end of quadrature
						[&](TREXIComputation xi) -> TREXIComputation
						{
							return rexi_int_fun(m, xi, fun_phi_tilde).real();
						},
						int_threshold
					);

		TREXIComputation imag = GaussQuadrature::integrate5_intervals_adaptive_recursive<TREXIComputation>(
						(TREXIComputation)i_int_start,	// start of quadrature
						(TREXIComputation)i_int_end,	// end of quadrature
						[&](TREXIComputation xi) -> TREXIComputation
						{
							return rexi_int_fun(m, xi, fun_phi_tilde).imag();
						},
						int_threshold
					);

		return std::complex<TStorage>(real, imag);
	}



	FunApproximation()
	{

	}

	FunApproximation(
			const std::string &i_function_name,
			TStorage i_h,
			int i_M
	)
	{
		setup(i_function_name, i_h, i_M);
	}

	void setup(
			const std::string &i_function_name,
			TStorage i_h,
			int i_M
	)
	{
		if (typeid(TREXIComputation) == typeid(double))
		{
			int_threshold = 1e-12;
		}
		else if (typeid(TREXIComputation) == typeid(__float128))
		{
			int_threshold = 1e-16;
		}
		else
		{
			FatalError("Type not supported!");
		}

		h = i_h;
		M = (i_M == -1 ? (TStorage)i_M/h : i_M);
		b.resize(i_M*2+1);

		__float128 asdf = DQStuff::fromString<TStorage>("3.14159265358979323846264338327950288");
		if (asdf - pi != 0)
			FatalError("Compiled constant not equal to string-induced constant!");

		if (i_function_name == "phi0")
		{
			/*
			 * This is only for testing purpose, since phi0 can be directly computed
			 *
			 * Since it's a dirac function, we have to compute the quadrature over a tiny interval at [1-eps;1+eps]
			 */

			auto fun = [&](TREXIComputation i_xi)	-> std::complex<TREXIComputation>
			{
				return pi2;
			};
			for (int m = -i_M; m <= i_M; m++)
			{
				b[m+M] = rexi_int_fun(m, 1.0/pi2, fun);
			}
		}
		else if (i_function_name == "phi1")
		{
			/*
			 * \f$
			 * 		\phi_1(\xi) = 2\pi for -1/(2*\pi) <= \xi <= 0
			 * \f$
			 *
			 * https://www.wolframalpha.com/input/?i=Fourier+transform+(exp(i*x%2F(2*pi))-1)%2F(i*x)*sqrt(2*pi)
			 */
			auto fun = [&](TREXIComputation i_xi)	-> std::complex<TREXIComputation>
			{
				return pi2;
			};

			for (int m = -i_M; m <= i_M; m++)
			{
				std::complex<TREXIComputation> t = computeREXICoefficient(m, fun, -1.0/pi2, 0);
				 b[m+M] = std::complex<TStorage>(t.real(), t.imag());
			}
		}
		else if (i_function_name == "phi2")
		{
			/*
			 * \f$
			 * 		\phi_2(\xi) = 2\pi for -1/(2*\pi) <= \xi <= 0
			 * \f$
			 *
			 * https://www.wolframalpha.com/input/?i=Fourier+transform+(exp(i*x%2F(2*pi))-1)%2F(i*x)*sqrt(2*pi)
			 */
			auto fun = [&](TREXIComputation i_xi)	-> std::complex<TREXIComputation>
			{
				return pi2;
			};

			for (int m = -i_M; m <= i_M; m++)
				b[m+M] = computeREXICoefficient(m, fun, -1.0/pi2, 0);
		}
		else
		{
			FatalError(std::string("Function ")+i_function_name+std::string(" is not supported"));
		}
	}



	void print()
	{
		for (std::size_t i = 0; i < b.size(); i++)
		{
			std::cout << "b[" << i << "]: " << b[i] << std::endl;
		}
	}


	complexREXIComputation approx(
			TREXIComputation i_x
	)
	{
		/// \f$ \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)} \f$

		complexREXIComputation sum = 0;
		for (int m = -M; m < M+1; m++)
		{
			sum += b[m+M] * ga.approxGaussian(i_x+((TREXIComputation)m)*h, h);
		}
		return sum;
	}
};



#endif
