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
#include <libmath/DQStuff.hpp>



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

	void setup_constvars()
	{
#if SWEET_QUADMATH
		if (typeid(T) == typeid(__float128))
		{
			eps_phi = 1e-10;
			eps_ups= 1e-10;

			static char *sp;

			pi2 = (__float128)2.0*strtoflt128("3.1415926535897932384626433832795029", &sp);
		}
		else
#endif
		if (typeid(T) == typeid(double))
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
	REXIFunctions()	:
		phi_id(0)
	{
		setup_constvars();
	}


public:
	REXIFunctions(const std::string &i_function_name)	:
		phi_id(0)
	{
		setup_constvars();
		setup(i_function_name);
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
		else if (i_function_name  == "phi4")
			phi_id = 4;
		else if (i_function_name  == "phi5")
			phi_id = 5;

		else if (i_function_name  == "ups1")
			phi_id = 101;
		else if (i_function_name  == "ups2")
			phi_id = 102;
		else if (i_function_name  == "ups3")
			phi_id = 103;

		//Semi-Lag phi functions (phi0 factored out) - see sl-rexi paper
		else if (i_function_name  == "psi1")
			phi_id = 1001;
		else if (i_function_name  == "psi2")
			phi_id = 1002;
		else if (i_function_name  == "psi3")
			phi_id = 1003;

		else
			FatalError("This phi function is not supported!");

		if (phi_id == 1 || phi_id == 2 || phi_id == 3 || phi_id == 101 || phi_id == 102 || phi_id == 103 || phi_id == 1002 || phi_id == 1003)
		{
			if (typeid(T) == typeid(double))
			{
				std::cout << "**************************************************************" << std::endl;
				std::cout << "* WARNING: " << i_function_name << " typically requires __float128 precision!" << std::endl;
				std::cout << "**************************************************************" << std::endl;
//				FatalError("Seriously, you shouldn't use me with only double precision!");
			}
		}
	}


#if SWEET_QUADMATH
	/**************************************************************
	 * __float128 TYPES
	 **************************************************************/
	static std::complex<__float128> l_expcplx(const std::complex<__float128> &i_value)
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
#endif


	/**************************************************************
	 * double TYPES
	 * Might suffer of numerical double precision limited effects
	 **************************************************************/

	std::complex<double> l_expcplx(const std::complex<double> &i_value)
	{
		return std::exp(i_value);
	};

	T l_sqrt(const double &i_value)
	{
		return std::sqrt(i_value);
	};

	std::complex<double> l_sqrtcplx(const std::complex<double> &i_value)
	{
		return std::sqrt(i_value);
	};



	template <typename D>
	static
	std::complex<D> convert(
			const std::complex<T> &i_value
	)
	{
		std::complex<D> data;
		data.real(i_value.real());
		data.imag(i_value.imag());

		return data;
	}

	int factorial(int N)
	{
		int retval = 1;

		for (; N > 0; N--)
			retval = retval * N;

		return retval;
	}

	std::complex<T> phi(
			int n,
			const std::complex<T> &z
	)
	{
		if (n == 0)
			return l_expcplx(z);

        if (std::abs((double)z.real()) < eps_phi && std::abs((double)z.imag()) < eps_phi)
				return (T)1.0/(T)factorial(n);

        return (phi(n-1, z) - (T)1.0/(T)factorial(n-1))/z;
	}

	/**
	 * Evaluate the function (phi/ups)
	 */
	std::complex<T> eval(
		const std::complex<T> &i_K
	)
	{
		std::complex<T> K = i_K;

		// for tests of numerical instability
		T lamdt = K.imag();
		if (lamdt < 0)
			lamdt = -lamdt;

		switch(phi_id)
		{
#if 0
		case 0:	// \phi_0
			K = l_expcplx(K);
			break;


		case 1:	// \phi_1
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


		case 2:	// \phi_2
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


		case 3:	// \phi_3
			if (lamdt*lamdt*lamdt < eps_phi)
			{
				K = 1.0/(2.0*3.0);
			}
			else
			{
				K = (l_expcplx(K) - std::complex<T>(1.0) - K - K*K)/(K*K*K);
			}
			break;
#endif

		case 101:	// \ups_1
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


		case 102:	// \ups_2
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


		case 103:	// \ups_3
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
			return phi(phi_id, i_K);
			//FatalError("This phi is not yet supported");
		}

		return K;
	}

};



#endif /* SRC_INCLUDE_REXI_REXIFunctions_HPP_ */
