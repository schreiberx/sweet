/*
 * GaussianApproximationNew.hpp
 *
 *  Created on: 29 Nov 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_NEW_HPP_
#define SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_NEW_HPP_

#include <sweet/sweetmath.hpp>
#include <iostream>
#include <complex>
#include <vector>



/**
 * This class provides the weights and coefficients for the
 * approximation of a Gaussian
 *
 * \f$
 * 	  exp(-(x*x)/(4*h*h))/sqrt(4 \pi)
 * \f$
 *
 * with a sum over complex rational functions.
 *
 * See e.g. Near optimal rational approximations of large data sets, Damle et. al.
 */
template <typename T = double>
class GaussianApproximationNew
{
	typedef std::complex<T> complex;
	T pi = (T)3.14159265358979323846264338327950288;
	T pi4 = (T)3.14159265358979323846264338327950288*(T)4;
	T sqrtpi4 = std::sqrt(pi4);


public:
	T mu;				///< average
	std::vector<complex> a;	///< weights for approximation
	int L;					///< 2*L+1 = number of weights

	GaussianApproximationNew(
			int i_L = 11	///< L
	)
	{
		L = i_L;

		a.resize(L*2+1);
		for (int i = 0; i < L*2+1; i++)
			a[i] = {0,0};

		mu = 0;

		do
		{
			// iterate over all L's
			for (int l = 0; l < L*2+1; l++)
			{
				/**
				 * Compute residual for current optimization pole
				 */
				T res = evalGaussian(x, 1)
			}
		} while (error > 1e-10);
	}

	T psi(T x)
	{
		return ((T)1.0/std::sqrt(pi4))*std::exp((-x*x)/((T)4.0*h*h));
	}



	/**
	 * directly evaluate basis function which is to be approximated
	 */
	T evalGaussian(
			T x,	///< x-coefficient for Gaussian basis function
			T h	///< h-coefficient for Gaussian basis function
	)
	{
		return std::exp(-(x*x)/(4*h*h))/std::sqrt(4.0*M_PIl);
	}


	/**
	 * evaluate approximation of Gaussian basis function
	 *
	 * with sum of complex rational functions
	 */
	T approxGaussian(
			T x,	///< x-coefficient for Gaussian basis function
			T h	///< h-coefficient for Gaussian basis function
	)
	{
		// scale x, since it depends linearly on h:
		// x^2 ~ h^2
		x /= h;

		T sum = 0;

		for (int l = 0; l < 2*L+1; l++)
		{
			int j = l-L;

			// WORKS with max error 7.15344e-13
			sum += (a[l]/(complex(0, x) + mu + complex(0, j))).real();
		}

		return sum;
	}

	void print()
	{
		std::cout << "mu: " << mu << std::endl;
		for (std::size_t i = 0; i < a.size(); i++)
		{
			std::cout << "a[" << i << "]: " << a[i] << std::endl;
		}
	}
};


#endif /* SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_ */
