/*
 * SPHIdentities.hpp
 *
 *  Created on: 26 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SPHERE_SPHERESPHIDENTITIES_HPP_
#define SRC_INCLUDE_SPHERE_SPHERESPHIDENTITIES_HPP_

#include <cassert>

class SphereSPHIdentities
{
public:

	inline
	static double D(double k, double m)
	{
		double n=k+1;
		assert(n >= 0);
//		if (n < 0)
//			n = -n-1;
		return ((2.0*n+1.0)*std::sqrt((n*n-m*m)/(4.0*n*n-1.0)));
	}

	inline
	static double E(double n, double m)
	{
		assert(n >= 0);
		return -n;
	}


	inline
	static double G(double n, double m)
	{
		n = n+1;
		return D(n-1,m) + E(n,m)*R(n-1,m);
	}


	inline
	static double H(double n, double m)
	{
		n = n-1;
		return E(n,m)*S(n+1,m);
	}

	inline
	static double R(double k, double m)
	{
		double n=k+1;

		if (n < 0)
			n = -n-1;

		if (std::abs(m) > n)
			return 0;

		assert (n*n-m*m >= 0);

		assert(n >= 0);
		return std::sqrt((n*n-m*m)/(4.0*n*n-1.0));
	}


	inline
	static double S(double k, double m)
	{
		double n=k-1;

		if (n < 0)
			n = -n-1;

		if (std::abs(m) > n)
			return 0;

		assert (n*n-m*m >= 0);

		assert(n >= 0);
		return std::sqrt(((n+1.0)*(n+1.0)-m*m)/((2.0*n+1.0)*(2.0*n+3.0)));
	}
	inline
	static double A(double k, double m)
	{
		double n = k+2.0;
		return R(n-1,m)*R(n-2,m);
	}

	inline
	static double B(double n, double m)
	{
		return R(n-1,m)*S(n,m) + S(n+1,m)*R(n,m);
	}

	inline
	static double C(double k, double m)
	{
		double n = k-2.0;
		return S(n+1,m)*S(n+2,m);
	}


#if 0
	inline
	std::complex<double>
	static csqrt(const double &i_value)
	{
		if (i_value >= 0)
			return std::complex<double>(std::sqrt(i_value), 0);
		else
			return std::complex<double>(0, std::sqrt(-i_value));
	}


	inline
	static std::complex<double> Rc(double k, double m)
	{
		if (k < 0)
			k = -k-1;

		if (std::abs(m) > k)
			return 0;

		double n = k+1;

		return csqrt((n*n-m*m)/(4.0*n*n-1.0));
	}


	inline
	static std::complex<double> Sc(double k, double m)
	{
		if (k < 0)
			k = -k-1;

		if (std::abs(m) > k)
			return 0;


		double n = k-1;

//		if (n < 0)
//			n = -n-1;

		return csqrt(((n+1.0)*(n+1.0)-m*m)/((2.0*n+1.0)*(2.0*n+3.0)));
	}


	inline
	static std::complex<double> Ac(double k, double m)
	{
		double n = k+2.0;
		return Rc(n-1,m)*Rc(n-2,m);
	}

	inline
	static std::complex<double> Bc(double n, double m)
	{
		return Rc(n-1,m)*Sc(n,m) + Sc(n+1,m)*Rc(n,m);
	}

	inline
	static std::complex<double> Cc(double k, double m)
	{
		double n = k-2.0;
		return Sc(n+1,m)*Sc(n+2,m);
	}
#endif

};


#endif /* SRC_INCLUDE_SPHERE_SPHERESPHIDENTITIES_HPP_ */
