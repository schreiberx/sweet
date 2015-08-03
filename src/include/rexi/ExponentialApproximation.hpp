/*
 * ExponentialApproximation.hpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_REXI_EXPONENTIALAPPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_EXPONENTIALAPPROXIMATION_HPP_

#include <cmath>
#include "GaussianApproximation.hpp"

/**
 * This class computes an approximation of an exponential e^{ix} by
 * a combination of Gaussians.
 *
 * The approximation is given by
 *
 * $$
 * 		e^{ix} - \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)}
 * $$
 *
 * see eq. (3.1) in
 * "A high-order time-parallel scheme for solving wave propagation problems via the direct construction of an approximate time-evolution operator", Haut et.al.
 *
 * Here, the coefficents b_m are given by the coefficents c_m (Yes, it's the c_m here!)
 * in the equation below equation (3.4):
 *
 * 		c_m := int_{-1/(2h)}^{1/2h}  exp(-2 pi i m h xi) F(xi)/Psi_h(xi) d xi
 *
 * with F and Psi the functions f and psi in Fourier space
 */
class ExponentialApproximation
{
	typedef std::complex<double> complex;

	GaussianApproximation ga;

	double h;
	int M;

public:
	std::vector<complex> b;

	ExponentialApproximation(
			double i_h,
			int i_M = -1
	)
	{
		h = i_h;
		M = (i_M == -1 ? (double)i_M/h : i_M);

		b.resize(i_M*2+1);

		for (int m = -i_M; m < i_M; m++)
		{
			/*
			 * See Section 3.1 for approx. of general function
			 *
			 * Note, that in the following we use the Fourier transformation
			 *  F(f(x)) := (xi) -> int_-inf^inf {  exp(-i*2*Pi*x*xi)  } dx
			 * hence, with 2*Pi in the exponent
			 *
			 * STEP 1: specialize on f(x) = e^{ix}
			 *
			 * In Fourier space, the function f(x) is given by
			 * 	F(xi) = \delta( xi-1.0/(2 Pi) )
			 * with \delta the Kronecker delta
			 *
			 * This simplifies the equation
			 *    c_m := h * int_{-1/(2h)}^{1/(2h)}  exp(-2 pi i m h xi) F(xi)/Psi_h(xi)   d xi
			 * to
			 *    c_m := h * exp(-i m h) / Psi_h( 1/(2*Pi) )
			 *
			 *
			 * STEP 2: Compute Psi_h( 1/(2*Pi) )
			 *
			 * Furthermore, Psi_h(xi) is given by
			 *
			 *    Psi_h(xi) := 1/sqrt(2) * e^{-(2*pi*h*xi)^2}
			 *
			 * and restricting it to xi=1/(2*pi) (see Kronecker delta above), yields
			 *
			 *    Psi_h(1/2Pi) := 1/sqrt(2) * e^{-h^2}
			 *
			 * URL:
			 * http://www.wolframalpha.com/input/?i=FourierTransform[1%2Fsqrt%284*pi%29*exp%28-x^2%2F%284*h^2%29%29%2C+x%2C+\[Omega]%2C+FourierParameters+-%3E+{0%2C+2+Pi}]
			 * RESULT: exp(-4 h^2 π^2 ω^2)/sqrt(1/h^2)
			 *         = h * exp(-(2 h π ω)^2)
			 * and with xi = 1/(2 Pi):
			 * 			h * exp(-h^2)
			 *
			 * STEP 3: combine (1) and (2):
			 *
			 * This simplifies c_m to
			 *
			 *    c_m := e^{h^2} * e^{-i*m*h}
			 *
			 * Let's hope, that these equations are right.
			 */

			b[m+M] = std::exp(i_h*i_h)*std::exp(-complex(0, 1)*((double)m*(double)h));
		}
	}


	complex evalExponential(
			double i_x
	)
	{
		return std::exp(complex(0,1)*i_x);
	}


	complex approxExponentialEvalGaussian(
			double i_x
	)
	{
		complex sum = 0;

		// \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)}

		for (int m = -M; m < M; m++)
		{
			sum += b[m+M] * ga.evalGaussian(i_x+((double)m)*h, h);
		}
		return sum;
	}

	complex approxExponentialApproxGaussian(
			double i_x
	)
	{
		complex sum = 0;

		// \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)}

		for (int m = -M; m < M; m++)
		{
			sum += b[m+M] * ga.approxGaussian(i_x+((double)m)*h, h);
		}
		return sum;
	}
};



#endif /* SRC_INCLUDE_REXI_EXPONENTIALAPPROXIMATION_HPP_ */
