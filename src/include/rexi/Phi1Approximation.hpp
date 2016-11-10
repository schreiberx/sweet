/*
 * PhiApproximation.hpp
 *
 *  Created on: 10 Nov 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_REXI_PHI1APPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_PHI1APPROXIMATION_HPP_

#include <sweet/sweetmath.hpp>
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
 */
class Phi1Approximation
{
	typedef std::complex<double> complex;

	GaussianApproximation ga;

	GaussQuadrature gQ;

	double h;
	int M;

public:
	std::vector<complex> b;

	Phi1Approximation(
			double i_h,
			int i_M
	)
	{
		h = i_h;
		M = (i_M == -1 ? (double)i_M/h : i_M);

		b.resize(i_M*2+1);

		for (int m = -i_M; m < i_M+1; m++)
		{
			double real = gQ.integrate5_intervals_adaptive(
							//std::max(-1.0/(2.0*h), -1.0/(2.0*M_PI)),	// start of quadrature
							-1.0/(2.0*M_PI),	// start of quadrature
							0,		// end of quadrature
							[&](double xi) -> double
							{
								return (h*std::exp(-2.0*M_PI*complex(0.0,1.0)*(double)m*h*xi)*
										(
												(2.0*M_PI)
												/(h*std::exp(-std::pow(2.0*h*M_PI*xi, 2.0)))
										)).real();
							}
						);

			double imag = gQ.integrate5_intervals_adaptive(
							//std::max(-1.0/(2.0*h), -1.0/(2.0*M_PI)),	// start of quadrature
							-1.0/(2.0*M_PI),	// start of quadrature
							0,		// end of quadrature
							[&](double xi) -> double
							{
								return (h*std::exp(-2.0*M_PI*complex(0.0,1.0)*(double)m*h*xi)*
										(
												(2.0*M_PI)
												/(h*std::exp(-std::pow(2.0*h*M_PI*xi, 2.0)))
										)).imag();
							}
						);

			b[m+M] = std::complex<double>(imag, real);	/// TODO: REAL AND IMAG PARTS ARE SWAPPED!!!
//			b[m+M] = std::complex<double>(real, imag);

		}
	}


	void print()
	{
		for (std::size_t i = 0; i < b.size(); i++)
		{
			std::cout << "b[" << i << "]: " << b[i] << std::endl;
		}
	}

	complex eval(
			double i_x
	)
	{
		return (std::exp(complex(0,1)*i_x)-1.0)/i_x;
	}


	complex approx(
			double i_x
	)
	{
		complex sum = 0;

		/// \f$ \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)} \f$

		for (int m = -M; m < M+1; m++)
		{
			sum += b[m+M] * ga.approxGaussian(i_x+((double)m)*h, h);
		}
		return sum;
	}
};



#endif
