/*
 * ExponentialApproximation.hpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_INCLUDE_REXI_REXI_TERRY_EXPONENTIALAPPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_REXI_TERRY_EXPONENTIALAPPROXIMATION_HPP_

#include <libmath/DQStuff.hpp>
#include <cmath>
#if SWEET_QUADMATH
	#include <quadmath.h>
#endif
#include <rexi/REXI_Terry_GaussianApproximation.hpp>



/**
 * This class computes an approximation of an exponential \f$ e^{ix} \f$ by
 * a combination of Gaussians.
 *
 * The approximation is given by
 *
 * \f$
 * 		e^{ix} - \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)}
 * \f$
 *
 * see eq. (3.1) in
 * "A high-order time-parallel scheme for solving wave propagation problems via the direct construction of an approximate time-evolution operator", Haut et.al.
 *
 * Here, the coefficents b_m are given by the coefficents \f$ c_m \f$ (Yes, it's the \f$ c_m \f$ here!)
 * in the equation below equation (3.4):
 *
 * \f$
 * 		c_m := \int_{-1/(2h)}^{1/2h}  exp(-2 \pi i m h \xi) F(\xi)/\psi_h(\xi) d \xi
 * \f$
 *
 * with F and \f$ \psi \f$ the functions f and \f$ \psi \f$ in Fourier space
 */
template <typename T>
class REXI_Terry_ExponentialApproximation
{
	typedef std::complex<T> TComplex;

	REXI_Terry_GaussianApproximation<T> ga;

	T h;
	int M;

public:
	std::vector<TComplex> b;



	REXI_Terry_ExponentialApproximation(
			T i_h,
			int i_M
	)
	{
		h = i_h;
		M = (i_M == -1 ? (T)(i_M/h+1) : i_M);

		b.resize(i_M*2+1);

		for (int m = -i_M; m < i_M+1; m++)
		{
			/*
			 * See Section 3.1 for approx. of general function
			 *
			 * Note, that in the following we use the Fourier transformation
			 *  \f$ F(f(x)) := (xi) -> \int_{-\inf}^\inf {  exp(-i*2*\pi*x*\xi)  } dx \f$
			 * hence, with 2*Pi in the exponent
			 *
			 * STEP 1: specialize on \f$ f(x) = e^{i x} \f$
			 *
			 * In Fourier space, the function f(x) is given by
			 *
			 * \f$
			 * 		F(\xi) = \delta( \xi-1.0/(2 \pi) )
			 * \f$
			 *
			 * with \delta the Kronecker delta
			 *
			 * This simplifies the equation
			 *
			 * \f$
			 *    c_m := h * int_{-1/(2 h)}^{1/(2 h)}  exp(-2 \pi i m h \xi) F(\xi)/Psi_h(\xi)   d \xi
			 * \f$
			 *
			 * to
			 *
			 * \f$
			 *    c_m := h * exp(-i m h) / Psi_h( 1/(2*\pi) )
			 * \f$
			 *
			 *
			 * STEP 2: Compute \f$ Psi_h( 1/(2*\pi) ) \f$
			 *
			 * Furthermore, \psi_h(\xi) is given by
			 *
			 * \f$
			 *    \psi_h(xi) := 1/sqrt(2) * e^{-(2*\pi*h*\xi)^2}
			 * \f$
			 *
			 * and restricting it to xi=1/(2*pi) (see Kronecker delta above), yields
			 *
			 * \f$
			 *    \psi_h(1/2 \pi) := 1/\sqrt{2} * e^{-h^2}
			 * \f$
			 *
			 * URL:
			 * http://www.wolframalpha.com/input/?i=FourierTransform[1%2Fsqrt%284*pi%29*exp%28-x^2%2F%284*h^2%29%29%2C+x%2C+\[Omega]%2C+FourierParameters+-%3E+{0%2C+2+Pi}]
			 * RESULT: \f$ exp(-4 h^2 \pi^2 \omega^2)/\sqrt(1/h^2)
			 *         = h * exp(-(2 h \pi \omega)^2) \f$
			 * and with \f$ \xi = 1/(2 \pi) \f$:
			 * 			\f$ h * exp(-h^2) \f$
			 *
			 * STEP 3: combine (1) and (2):
			 *
			 * This simplifies c_m to
			 *
			 * \f$
			 *    c_m := e^{h^2} * e^{-i*m*h}
			 * \f$
			 */

			b[m+M] = DQStuff::exp(
					TComplex(h*h, -(T)m*h)
				);
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
	TComplex eval(
			T i_x
	)
	{
		return DQStuff::expIm(i_x);
	}


	TComplex approx(
			T i_x
	)
	{
		TComplex sum = 0;

		/// \f$ \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)} \f$

		for (int m = -M; m < M+1; m++)
		{
			sum += DQStuff::convertComplex<T>(b[m+M]) * ga.approxGaussian(i_x+(T)m*h, h);
		}
		return sum;
	}
};



#endif
