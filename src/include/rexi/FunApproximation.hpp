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
 * This class computes an approximation of an exponential \f$ phi1(x) := (e^{ix}-1)/ix \f$ by
 * a combination of Gaussians.
 *
 * We use a numerical quadrature for this:
 *
 * \f$
 * 		c_m := \int_{-1/(2h)}^{1/2h}  exp(-2 \pi i m h \xi) phi1(\xi)/\psi_h(\xi) d \xi
 * \f$
 *
 * With
 *
 * \f$
 * 		\phi_1(\xi) = 2\pi for -1/(2*\pi) <= \xi <= 0
 * \f$
 *
 * https://www.wolframalpha.com/input/?i=Fourier+transform+(exp(i*x%2F(2*pi))-1)%2F(i*x)*sqrt(2*pi)
 *
 * For more information, see Terry Haut et al. paper
 * "A high-order time-parallel scheme for solving wave propagation problems via
 * the direct construction of an approximate time-evolution operator"
 */
template <
	typename TEvaluation = __float128,	///< evaluation accuracy of coefficients
	typename TStorageAndProcessing = double	///< storage precision of coefficients - use quad precision per default
>
class FunApproximation
{
	typedef std::complex<TEvaluation> complexEvaluation;
	typedef std::complex<TStorageAndProcessing> complexStorage;

	GaussianApproximation<TEvaluation, TStorageAndProcessing> ga;


	TStorageAndProcessing h;
	int M;

	const TStorageAndProcessing pi = DQStuff::fromString<TStorageAndProcessing>("3.14159265358979323846264338327950288");
	const TStorageAndProcessing pi2 = pi*DQStuff::fromString<TStorageAndProcessing>("2.0");
	const TStorageAndProcessing pi4 = pi*DQStuff::fromString<TStorageAndProcessing>("4.0");
	const TStorageAndProcessing sqrtpi4 = DQStuff::sqrt(pi4);

public:
	std::vector<complexStorage> b;

	FunApproximation(
			const std::string &i_function_name,
			TStorageAndProcessing i_h,
			int i_M
	)
	{
		h = i_h;
		M = (i_M == -1 ? (TStorageAndProcessing)i_M/h : i_M);

		b.resize(i_M*2+1);


		__float128 asdf = DQStuff::fromString<TStorageAndProcessing>("3.14159265358979323846264338327950288");
		if (asdf - pi != 0)
			FatalError("Compiled constant not equal to string-induced constant!");

		if (i_function_name == "phi1")
		{
			for (int m = -i_M; m < i_M+1; m++)
			{
				TStorageAndProcessing real = GaussQuadrature::integrate5_intervals_adaptive_recursive<TStorageAndProcessing>(
								(TStorageAndProcessing)-1.0/pi2,	// start of quadrature
								(TStorageAndProcessing)0.0,		// end of quadrature
								[&](TStorageAndProcessing xi) -> TStorageAndProcessing
								{
									return (h*DQStuff::expIm(-pi2*m*h*xi)*
											(TStorageAndProcessing)(
													pi2/(h*DQStuff::exp(-(pi2*h*xi)*(pi2*h*xi)))
											)).real();
								}
							);

				TStorageAndProcessing imag = GaussQuadrature::integrate5_intervals_adaptive_recursive<TStorageAndProcessing>(
								//std::max(-1.0/(2.0*h), -1.0/(2.0*M_PI)),	// start of quadrature
								-1.0/(2.0*M_PI),	// start of quadrature
								0.0,		// end of quadrature
								[&](TStorageAndProcessing xi) -> TStorageAndProcessing
								{
									return (h*DQStuff::expIm(-pi2*m*h*xi)*
											(TStorageAndProcessing)(
													pi2/(h*DQStuff::exp(-DQStuff::pow(pi2*h*xi, (TStorageAndProcessing)2.0)))
											)).imag();
								}
							);

				b[m+M] = std::complex<TStorageAndProcessing>(imag, real);	/// TODO: REAL AND IMAG PARTS ARE SWAPPED - WHY?!!!
			}
		}
		else if (i_function_name == "phi1")
		{

		}
	}


	void print()
	{
		for (std::size_t i = 0; i < b.size(); i++)
		{
			std::cout << "b[" << i << "]: " << b[i] << std::endl;
		}
	}

	static
	complexEvaluation eval(
			TEvaluation i_x
	)
	{
		return (DQStuff::expIm(i_x)-(TEvaluation)1.0)/i_x;
	}


	complexEvaluation approx(
			TEvaluation i_x
	)
	{
		complexEvaluation sum = 0;

		/// \f$ \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)} \f$

		for (int m = -M; m < M+1; m++)
		{
			sum += b[m+M] * ga.approxGaussian(i_x+((complexEvaluation)m)*h, h);
		}
		return sum;
	}
};



#endif
