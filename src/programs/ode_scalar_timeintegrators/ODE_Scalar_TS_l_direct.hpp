/*
 * ODE_Scalar_TS_l_direct.hpp
 *
 *  Created on: 04 Oct 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_DIRECT_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_DIRECT_HPP_

#include <rexi/EXPFunctions.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "ODE_Scalar_TS_interface.hpp"

/*
#if SWEET_QUADMATH
	#include <quadmath.h>
#endif
*/

template <typename T>
class ODE_Scalar_TS_l_direct	: public ODE_Scalar_TS_interface<T>
{
	SimulationVariables &simVars;


	int phi_id;

#if 0
#if SWEET_QUADMATH && 0
	typedef __float128 T;
#else
	typedef double T;
#endif
#endif

	EXPFunctions<double> rexiFunctions;

#if 0

#if SWEET_QUADMATH

	typedef __float128 T;

	static std::complex<T> l_expcplx(std::complex<T> &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = cexpq(z);

		return std::complex<double>(crealq(val), cimagq(val));
	}

	inline
	static T l_sqrt(T &i_value)
	{
		return sqrtq(i_value);
	}

	static const std::complex<T> l_sqrtcplx(const std::complex<T> &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = csqrtq(z);

		return std::complex<double>(crealq(val), cimagq(val));
	}

	inline
	static T eps_phi()
	{
		return 1e-10;
	}

	inline
	static T eps_ups()
	{
		return 1e-10;
	}

	inline
	static T pi2()
	{
		static char *sp;
		static T retval = (T)2.0*strtoflt128("3.1415926535897932384626433832795029", &sp);
		return retval;
	}

#else

	/*
	 * Double precision
	 * Might suffer of numerical double precision limited effects
	 */
	typedef double T;

	std::complex<T> l_expcplx(std::complex<T> &i_value)
	{
		return std::exp(i_value);
	};

	T l_sqrt(T &i_value)
	{
		return l_sqrt(i_value);
	};

	std::complex<T> l_sqrtcplx(std::complex<T> &i_value)
	{
		return std::exp(i_value);
	};


	static T eps_phi()
	{
		return 1e-10;
	}

	static T eps_ups()
	{
		return 1e-10;
	}

	static T pi2()
	{
		return (T)2.0*(T)M_PI;
	}

#endif

#endif

public:
	ODE_Scalar_TS_l_direct(
			SimulationVariables &i_simVars
		)
		:
		simVars(i_simVars)
	{
	}

	void run_timestep(
			///T &io_u,	///< prognostic variables
			ScalarDataArray &io_u,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	)
	{

		std::size_t N = io_u.number_of_elements;
		ScalarDataArray lL = this->lambda_L(i_dt, i_simulation_timestamp);

		for (std::size_t i = 0; i < N; i++)
		{

			///std::complex<double> lambda = this->param_function_L;
			std::complex<double> lambda = lL.get(i);

			std::complex<double> K = rexiFunctions.eval(lambda * i_dt);


#if !SWEET_SCALAR_COMPLEX
			double u = io_u.get(i);
			u *= K.real();
#else
			std::complex<double> u = io_u.get(i);
			u *= K;
#endif
			io_u.set(i, u);
		}
	}

	void setup(
			const std::string &i_function_name = "phi0"
	)
	{
		rexiFunctions.setup(i_function_name);
	}

	void setup(
			std::vector<double> i_L,
			std::vector<double> i_N,
			std::vector<double> i_extra,
			std::string i_model
		)
	{
		this->param_function_L = i_L;
		this->param_function_N = i_N;
		this->param_function_extra = i_extra;
		this->model = i_model;
	}

	virtual ~ODE_Scalar_TS_l_direct()
	{
	}
};

#endif /* SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_DIRECT_HPP_ */
