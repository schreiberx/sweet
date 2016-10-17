/*
 * SPHIdentities.hpp
 *
 *  Created on: 26 Aug 2016
 *      Author: martin
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
//		if (n < 0)
//			n = -n-1;
		if (n < 0)
			return 0;
		if (n < std::abs(m))
			return 0;
		if (n*n-m*m < 0)
			return 0;

		assert(n >= 0);
		return std::sqrt((n*n-m*m)/(4.0*n*n-1.0));
	}

	inline
	static double S(double k, double m)
	{
		double n=k-1;
//		if (n < 0)
//			n = -n-1;
		if (n < 0)
			return 0;
		if (n < std::abs(m))
			return 0;
		if (n*n-m*m < 0)
			return 0;

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





	inline
	static std::complex<double> cR(
			const std::complex<double> &k,
			const std::complex<double> &m
	)
	{
		std::complex<double> n=k+1.0;

		return std::sqrt((n*n-m*m)/(4.0*n*n-1.0));
	}

	inline
	static std::complex<double> cS(
			const std::complex<double> &k,
			const std::complex<double> &m
	)
	{
		std::complex<double> n=k-1.0;

		return std::sqrt(((n+1.0)*(n+1.0)-m*m)/((2.0*n+1.0)*(2.0*n+3.0)));
	}


	inline
	static std::complex<double> cA(
			const std::complex<double> &k,
			const std::complex<double> &m
	)
	{
		std::complex<double> n = k+2.0;

		return cR(n-1.0, m)*cR(n-2.0,m);
	}

	inline
	static std::complex<double> cB(
			const std::complex<double> &n,
			const std::complex<double> &m
	)
	{
		return cR(n-1.0,m)*cS(n,m) + cS(n+1.0, m)*cR(n,m);
	}

	inline
	static std::complex<double> cC(
			const std::complex<double> &k,
			const std::complex<double> &m
	)
	{
		std::complex<double> n = k-2.0;
		return cS(n+1.0,m)*cS(n+2.0,m);
	}
};


#endif /* SRC_INCLUDE_SPHERE_SPHERESPHIDENTITIES_HPP_ */
