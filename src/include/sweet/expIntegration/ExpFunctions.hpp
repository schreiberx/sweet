/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_EXPFUNCTIONS_HPP_
#define SRC_INCLUDE_EXPFUNCTIONS_HPP_

#include <iostream>
#include <complex>
#include <typeinfo>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/libmath/DQStuff.hpp>

#define EXP_FUNCTIONS_MAX_ITERS_DEFAULT 20

/*
 * 0: Series
 * 1: Cauchy
 */
#define EXP_FUNCTIONS_PHI_SPECIAL_EVAL_TYPE	0

#define EXP_FUNCTIONS_PHI_SPECIAL_THRESHOLD	1

/**
 * This class implements various EXP functions.
 *
 * See Cox-Matthews paper for the provided EXP functions
 */

namespace sweet
{

template <typename T = double>
class ExpFunctions
{
public:
	ErrorBase error;

private:
	typedef std::complex<T> CT;

	enum fun_id_enum
	{
		INVALID = 0,

		PHI0 = 1,
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

public:
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
			SWEETError("Type not supported");
		}
	}



public:
	ExpFunctions()	:
		function_id(INVALID)
	{
		setup_constvars();
	}


public:
	bool setup(
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

		/*
		 * Semi-Lag phi functions (phi0 factored out) - see sl-rexi paper
		 * This is not the psi from some of Martin's presentation slides!
		 */
		else if (i_function_name  == "psi1")
			function_id = PSI1;
		else if (i_function_name  == "psi2")
			function_id = PSI2;
		else if (i_function_name  == "psi3")
			function_id = PSI3;

		else
		{
			std::ostringstream ss;
			ss << "The function '" << i_function_name << "' is not supported!";
			error.set(ss.str());
			return false;
		}

		return true;
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
		const std::complex<T> &K,
		int max_iters = EXP_FUNCTIONS_MAX_ITERS_DEFAULT
	)
	{
		std::complex<T> pow_K = 1;
		T facn = factorial(n);
		std::complex<T> retval = pow_K/facn;

		for (int i = 1; i < max_iters; i++)
		{
			pow_K *= K;
			facn *= (n+i);

			retval += pow_K/facn;
		}

		return retval;
	}


	/*
	 * \phi_N: default caller switching between recursive formulation and series
	 */
	std::complex<T> phiN(
			int N,
			const std::complex<T> &z,
			int max_iters = EXP_FUNCTIONS_MAX_ITERS_DEFAULT
	)
	{
        if (std::abs(z) < EXP_FUNCTIONS_PHI_SPECIAL_THRESHOLD)
			return phiNSeries(N, z, max_iters);

        return phiNRec(N, z);
	}



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

		SWEETError("ups number not supported!");
		return 0;
	}

	std::complex<T> ups1Series(
			const std::complex<T> &K,
			int max_iters = EXP_FUNCTIONS_MAX_ITERS_DEFAULT
    )
	{
		std::complex<T> pow_K = 1;
		T fac_denom = 6;	//factorial(3);
		std::complex<T> retval = pow_K/fac_denom;

		for (int l = 1; l < max_iters; l++)
		{
			pow_K *= K;
			fac_denom *= (l+3);
			retval += pow_K*(l+1.0)*(l+1.0)/fac_denom;
		}
		return retval;
	}

	std::complex<T> ups2Series(
			const std::complex<T> &K,
			int max_iters = EXP_FUNCTIONS_MAX_ITERS_DEFAULT
    )
	{
		std::complex<T> pow_K = 1;
		T facn = 6;//factorial(3);
		std::complex<T> retval = 1.0/2.0;

		retval += (K-2.0)*pow_K/facn;

		for (int l = 1; l < max_iters; l++)
		{
			pow_K *= K;
			facn *= (l+3);
			retval += (K-2.0)*pow_K/facn;
		}
		return retval;
	}

	std::complex<T> ups3Series(
			const std::complex<T> &K,
			int max_iters = EXP_FUNCTIONS_MAX_ITERS_DEFAULT
    )
	{
		std::complex<T> pow_K = 1;
		T facn = 6;//factorial(3);
		std::complex<T> retval = -1.0/2.0;

		retval += (4.0-K)*pow_K/facn;

		for (int l = 1; l < max_iters; l++)
		{
		    pow_K *= K;
		    facn *= (l+3);
		    retval += (4.0-K)*pow_K/facn;
		}
		return retval;
	}

	std::complex<T> upsNSeries(
			int n,
			const std::complex<T> &z,
			int max_iters = EXP_FUNCTIONS_MAX_ITERS_DEFAULT
    )
	{
		switch(n)
		{
		case 1:
			return ups1Series(z, max_iters);

		case 2:
			return ups2Series(z, max_iters);

		case 3:
			return ups3Series(z, max_iters);
		}

		SWEETError("This Upsilon function is not supported!");
		return 0;
	}


	/*
	 * \phi_N: default caller switching between recursive formulation and series
	 */
	std::complex<T> upsN(
			int N,
			const std::complex<T> &z,
			int max_iters = EXP_FUNCTIONS_MAX_ITERS_DEFAULT
	)
	{
		T linf = z.real()*z.real() + z.imag()*z.imag();
		if (linf < EXP_FUNCTIONS_PHI_SPECIAL_THRESHOLD)
			return upsNSeries(N, z, max_iters);

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
		switch(n)
		{
		case 1:	// SL psi1
			// psi1(z) = phi1(-z)
			return phiN(1, -i_K);
			break;


		case 2:	// SL psi2
			// psi2(z)=-phi2(-z)+phi1(-z)
			return -phiN(2, -i_K) + phiN(1, -i_K);
			break;


		case 3:	// SL psi3
#if 1
			SWEETError("TODO: Redo this with e.g. series treatment");
#else
			if (lamdt < eps_phi)
				//					if (lamdt*lamdt*lamdt < expFunctions.eps_phi)
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
			SWEETError("This psi number is not yet supported");
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
			{
				std::ostringstream ss;
				ss << "The phi function with id '" << function_id << "' is not supported!";
				error.set(ss.str());
				SWEETError(ss.str().c_str());
			}
		}

		return K;
	}

};

}

#endif
