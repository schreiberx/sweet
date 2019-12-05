/*
 * REXIFunctions.hpp
 *
 *  Created on: 10 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
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
template <typename T = double>
class REXIFunctions
{
	enum fun_id_enum
	{
		INVALID,

		PHI0,
		PHI1,
		PHI2,
		PHI3,
		PHI4,
		PHI5,

		UPS1,
		UPS2,
		UPS3,

		PSI1,
		PSI2,
		PSI3
	};

	fun_id_enum function_id;


public:
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
		function_id(INVALID)
	{
		setup_constvars();
	}


public:
	REXIFunctions(const std::string &i_function_name)	:
		function_id(INVALID)
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
			function_id = PHI0;
		else if (i_function_name  == "phi1")
			function_id = PHI1;
		else if (i_function_name  == "phi2")
			function_id = PHI2;
		else if (i_function_name  == "phi3")
			function_id = PHI3;
		else if (i_function_name  == "phi4")
			function_id = PHI4;
		else if (i_function_name  == "phi5")
			function_id = PHI5;

		else if (i_function_name  == "ups1")
			function_id = UPS1;
		else if (i_function_name  == "ups2")
			function_id = UPS2;
		else if (i_function_name  == "ups3")
			function_id = UPS3;

		// Semi-Lag phi functions (phi0 factored out) - see sl-rexi paper
		else if (i_function_name  == "psi1")
			function_id = PSI1;
		else if (i_function_name  == "psi2")
			function_id = PSI2;
		else if (i_function_name  == "psi3")
			function_id = PSI3;

		else
		{
			std::cout<<"REXIFunctions error: i_function_name ="<< i_function_name <<std::endl;
			FatalError("This phi function is not supported! : ");
		}

		
/*
		if (function_id == 1 || function_id == 2 || function_id == 3 || function_id == 101 || function_id == 102 || function_id == 103 || function_id == 1002 || function_id == 1003)
		{
			if (typeid(T) == typeid(double))
			{
				std::cout << "**************************************************************" << std::endl;
				std::cout << "* WARNING: " << i_function_name << " typically requires __float128 precision!" << std::endl;
				std::cout << "**************************************************************" << std::endl;
//				FatalError("Seriously, you shouldn't use me with only double precision!");
			}
		}
*/
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


	int factorial(
			int N
	)
	{
		int retval = 1;

		for (; N > 0; N--)
			retval = retval * N;

		return retval;
	}



	/*
	 * \phi_N: Recursive computations
	 *
	 * ATTENTION: There's a singularity close to 0!
	 */
	std::complex<T> phiNRec(
		int n,
		const std::complex<T> &z
	)
	{
		if (n == 0)
			return l_expcplx(z);

        return (phiNRec(n-1, z) - (T)1.0/(T)factorial(n-1))/z;
	}




	/*
	 * \phi_N: Series-based computation
	 *
	 * This avoids the singularity close to 0
	 */
	std::complex<T> phiNSeries(
		int n,
		const std::complex<T> &z
	)
	{
		std::complex<T> powz = 1.0;
		T facn = factorial(n);

		std::complex<T> retval = powz/facn;

		for (int i = 1; i < 20; i++)
		{
			powz *= z;
			facn *= (n+i);

			retval += powz/facn;
		}

		return retval;
	}



	/*
	 * \phi_N: default caller switching between recursive formulation and series
	 */
	std::complex<T> phiN(
			int N,
			const std::complex<T> &z
	)
	{
		T linf = z.real()*z.real() + z.imag()*z.imag();
        if (linf < 0.2)
			return phiNSeries(N, z);

        return phiNRec(N, z);
	}



#if 0
	std::complex<T> phi(
			int n,
			const std::complex<T> &z
	)
	{
		if (n == 0)
			return l_expcplx(z);

		T linf = z.real()*z.real() + z.imag()*z.imag();
        if (linf < eps_phi)
			return (T)1.0/(T)factorial(n);

        return (phi(n-1, z) - (T)1.0/(T)factorial(n-1))/z;
	}
#endif


	std::complex<T> upsNDirect(
			int n,
			const std::complex<T> &z
    )
	{
		switch(n)
		{
		case 1:
			// http://www.wolframalpha.com/input/?i=(-4-K%2Bexp(K)*(4-3K%2BK*K))%2F(K*K*K)
			return (-4.0-z+l_expcplx(z)*(4.0-3.0*z+z*z)) / (z*z*z);

		case 2:
			// http://www.wolframalpha.com/input/?i=(2%2BK%2Bexp(K)*(-2%2BK))%2F(K*K*K)
			return (2.0+z+l_expcplx(z)*(-2.0+z)) / (z*z*z);

		case 3:
			// http://www.wolframalpha.com/input/?i=(-4-3*K-K*K%2Bexp(K)*(4-K))%2F(K*K*K)
			return (-4.0-3.0*z-z*z+l_expcplx(z)*(4.0-z)) / (z*z*z);
		}

		FatalError("ups number not supported!");
	}


	std::complex<T> upsNSeries(
			int n,
			const std::complex<T> &z
    )
	{
		std::complex<T> powz;
		T facn;

		std::complex<T> retval;

		int niters = 20;

		switch(n)
		{
		case 1:
			powz = 1.0;
			facn = factorial(3);

			retval = powz/facn;
			for (int l = 1; l < niters; l++)
			{
		        powz *= z;
		        facn *= (l+3);
		        retval += powz*(l+1.0)*(l+1.0)/facn;
			}
			return retval;

		case 2:
            retval = 1.0/2.0;

            powz = 1.0;
            facn = factorial(3);

            retval += (z-2.0)*powz/facn;
			for (int l = 1; l < niters; l++)
			{
				powz *= z;
				facn *= (l+3);
				retval += (z-2.0)*powz/facn;
			}

			return retval;

		case 3:
            retval = -1.0/2.0;

            powz = 1.0;
            facn = factorial(3);

            retval += (4.0-z)*powz/facn;
			for (int l = 1; l < niters; l++)
			{
				powz *= z;
				facn *= (l+3);
				retval += (4.0-z)*powz/facn;
			}

			return retval;
		}

		FatalError("This Upsilon function is not supported!");
	}


	/*
	 * \phi_N: default caller switching between recursive formulation and series
	 */
	std::complex<T> upsN(
			int N,
			const std::complex<T> &z
	)
	{
		T linf = z.real()*z.real() + z.imag()*z.imag();
        if (linf < 0.2)
			return upsNSeries(N, z);

        return upsNDirect(N, z);
	}




	/*
	 * Semi-Lagrangian psi functions
	 */
	std::complex<T> psiN(
			int n,
			const std::complex<T> &i_K
	)
	{
		// for tests of numerical instability
//		T lamdt = i_K.real()*i_K.real() + i_K.imag()*i_K.imag();

		switch(n)
		{
		case 1:	// psi1
#if 1
			// psi1(z) = phi1(-z)
			return phiN(1, -i_K);
#else
			//
			if (lamdt < eps_phi)
			{
				return 1.0;
			}
			else
			{
				//psi1(z)=phi(-z)
				return (l_expcplx(i_K) - std::complex<T>(1.0))/i_K;
			}
#endif
			break;


		case 2:	// psi2
#if 1
			// psi2(z)=-phi2(-z)+phi1(-z)
			return -phiN(2, -i_K) + phiN(1, -i_K);
#else
			if (lamdt < eps_phi)
				//					if (lamdt*lamdt < rexiFunctions.eps_phi)
			{
				return 1.0/2.0;
			}
			else
			{
				//psi2(z)=-phi2(-z)+phi1(-z)
				return -(l_expcplx(i_K) - std::complex<T>(1.0) - i_K)/(i_K*i_K)
								+(l_expcplx(i_K) - std::complex<T>(1.0))/i_K;
			}
#endif
			break;


		case 3:	// psi3
#if 1
			FatalError("TODO: Redo this with e.g. series treatment");
#else
			if (lamdt < eps_phi)
				//					if (lamdt*lamdt*lamdt < rexiFunctions.eps_phi)
			{
				return 1.0/(2.0*3.0);
			}
			else
			{
				return (l_expcplx(i_K) - std::complex<T>(1.0) - i_K - i_K*i_K)/(i_K*i_K*i_K)
						- (l_expcplx(i_K) - std::complex<T>(1.0) - i_K)/(i_K*i_K)
						+ std::complex<T>(0.5)*(l_expcplx(i_K) - std::complex<T>(1.0))/i_K;
			}
#endif
			break;


		default:
			FatalError("This psi number is not yet supported");
		}

		return -1;
	}


	/**
	 * Evaluate the function (phi/ups)
	 */
	std::complex<T> eval(
		const std::complex<T> &i_K
	)
	{
		std::complex<T> K = i_K;

		switch(function_id)
		{
		case PHI0:
			return phiN(0, i_K);

		case PHI1:
			return phiN(1, i_K);

		case PHI2:
			return phiN(2, i_K);

		case PHI3:
			return phiN(3, i_K);

		case PHI4:
			return phiN(4, i_K);

		case PHI5:
			return phiN(5, i_K);



		case UPS1:
			return upsN(1, i_K);

		case UPS2:
			return upsN(2, i_K);

		case UPS3:
			return upsN(3, i_K);



		case PSI1:
			return psiN(1, i_K);

		case PSI2:
			return psiN(2, i_K);

		case PSI3:
			return psiN(3, i_K);

		default:
			FatalError("This phi is not yet supported");
		}

		return K;
	}

};



#endif /* SRC_INCLUDE_REXI_REXIFunctions_HPP_ */
