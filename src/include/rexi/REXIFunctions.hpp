/*
 * REXIFunctions.hpp
 *
 *  Created on: 10 Aug 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_REXI_REXIFUNCTIONS_HPP_
#define SRC_INCLUDE_REXI_REXIFUNCTIONS_HPP_

#include <iostream>
#include <typeinfo>
#include <quadmath.h>


/**
 * This class implements various REXI functions.
 *
 * See Cox-Matthews paper for the provided REXI functions
 */
template <typename T>
class REXIFunctions
{
public:
	int phi_id;

	T eps_phi;
	T eps_ups;
	T pi2;

public:
	REXIFunctions()	:
		phi_id(0)
	{
		if (typeid(T) == typeid(__float128))
		{
			eps_phi = 1e-10;
			eps_ups= 1e-10;

			static char *sp;
			pi2 = (__float128)2.0*strtoflt128("3.1415926535897932384626433832795029", &sp);
		}
		else if (typeid(T) == typeid(double))
		{
			eps_phi = 1e-10;
			eps_ups= 1e-10;
			pi2 = (double)2.0*(double)M_PI;
		}
		else
		{
			FatalError("Type not supported");
		}
	}



public:
	void setup(
			const std::string &i_function_name
	)
	{
		if (i_function_name  == "phi0")
			phi_id = 0;
		else if (i_function_name  == "phi1")
			phi_id = 1;
		else if (i_function_name  == "phi2")
			phi_id = 2;
		else if (i_function_name  == "phi3")
			phi_id = 3;

		else if (i_function_name  == "ups1")
			phi_id = 101;
		else if (i_function_name  == "ups2")
			phi_id = 102;
		else if (i_function_name  == "ups3")
			phi_id = 103;
	}


	/**************************************************************
	 * __float128 TYPES
	 **************************************************************/
	static std::complex<__float128> l_expcplx(std::complex<__float128> &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = cexpq(z);

		return std::complex<__float128>(crealq(val), cimagq(val));
	}

	inline
	static T l_sqrt(const __float128 &i_value)
	{
		return sqrtq(i_value);
	}

	static const std::complex<__float128> l_sqrtcplx(const std::complex<__float128> &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = csqrtq(z);

		return std::complex<double>(crealq(val), cimagq(val));
	}



	/**************************************************************
	 * double TYPES
	 * Might suffer of numerical double precision limited effects
	 **************************************************************/

	std::complex<double> l_expcplx(std::complex<double> &i_value)
	{
		return std::exp(i_value);
	};

	T l_sqrt(double &i_value)
	{
		return l_sqrt(i_value);
	};

	std::complex<double> l_sqrtcplx(std::complex<double> &i_value)
	{
		return std::exp(i_value);
	};





	/**
	 * Evaluate the function (phi/ups)
	 */
	std::complex<T> eval(
			const std::complex<T> &i_K
	)
	{
		std::complex<T> K = i_K;
		T lamdt = K.imag();
		if (lamdt < 0)
			lamdt = -lamdt;

		switch(phi_id)
		{
		case 0:	// phi0
			K = l_expcplx(K);
			break;


		case 1:	// phi1
			// http://www.wolframalpha.com/input/?i=(exp(i*x)-1)%2F(i*x)
			if (lamdt < eps_phi)
			{
				K = 1.0;
			}
			else
			{
				K = (l_expcplx(K) - std::complex<T>(1.0))/K;
			}
			break;


		case 2:	// phi2
			// http://www.wolframalpha.com/input/?i=(exp(i*x)-1-i*x)%2F(i*x*i*x)
			if (lamdt*lamdt < eps_phi)
			{
				K = 1.0/2.0;
			}
			else
			{
				K = (l_expcplx(K) - std::complex<T>(1.0) - K)/(K*K);
			}
			break;


		case 3:	// phi3
			if (lamdt*lamdt*lamdt < eps_phi)
			{
				K = 1.0/(2.0*3.0);
			}
			else
			{
				K = (l_expcplx(K) - std::complex<T>(1.0) - K - K*K)/(K*K*K);
			}
			break;


		case 101:	// ups1
			// http://www.wolframalpha.com/input/?i=(-4-K%2Bexp(K)*(4-3K%2BK*K))%2F(K*K*K)
			if (lamdt*lamdt*lamdt < eps_ups)
			{
				K = 1.0/(2.0*3.0);
			}
			else
			{
				K = (std::complex<T>(-4.0)-K+l_expcplx(K)*(std::complex<T>(4.0)-std::complex<T>(3.0)*K+K*K))/(K*K*K);
			}
			break;


		case 102:	// ups2
			// http://www.wolframalpha.com/input/?i=(2%2BK%2Bexp(K)*(-2%2BK))%2F(K*K*K)
			if (lamdt*lamdt*lamdt < eps_ups)
			{
				K = 1.0/(2.0*3.0);
			}
			else
			{
				K = (std::complex<T>(2.0)+K+l_expcplx(K)*(std::complex<T>(-2.0)+K))/(K*K*K);
			}
			break;


		case 103:	// ups3
			// http://www.wolframalpha.com/input/?i=(-4-3*K-K*K%2Bexp(K)*(4-K))%2F(K*K*K)
			if (lamdt*lamdt*lamdt < eps_ups)
			{
				K = 1.0/(2.0*3.0);
			}
			else
			{
				K = (std::complex<T>(-4.0)-std::complex<T>(3.0)*K-K*K+l_expcplx(K)*(std::complex<T>(4.0)-K))/(K*K*K);
			}
			break;


		default:
			FatalError("This phi is not yet supported");
		}

		return K;
	}

};



#endif /* SRC_INCLUDE_REXI_REXIFunctions_HPP_ */
