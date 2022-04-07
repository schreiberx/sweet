/*
 * SWE_Plane_TS_l_direct.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_DIRECT_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_DIRECT_HPP_

#include <rexi/EXPFunctions.hpp>
#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataGridMapping.hpp>
#include <sweet/plane/PlaneStaggering.hpp>

#include "SWE_Plane_TS_interface.hpp"

/*
#if SWEET_QUADMATH
	#include <quadmath.h>
#endif
*/

class SWE_Plane_TS_l_direct	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int phi_id;

#if SWEET_QUADMATH && 0
	typedef __float128 T;
#else
	typedef double T;
#endif

	EXPFunctions<T> rexiFunctions;

	PlaneDataGridMapping planeDataGridMapping;

	// Precompute Z = Z(k1, k2, Dt) = Q * e^{\Lambda*Dt} * Q^{-1}
	std::vector<std::vector<std::array<std::array<std::complex<T>, 3>, 3>>> Z;  // Z[k1][k2][0,1,2][0,1,2];
	double dt_precompute_phin= 0.; // store dt for which Z has been precomputed in order to check if it is necessary to recompute it.

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
	SWE_Plane_TS_l_direct(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);


	void run_timestep(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_cgrid(
			PlaneData_Spectral &io_h_pert,	///< prognostic variables
			PlaneData_Spectral &io_u,		///< prognostic variables
			PlaneData_Spectral &io_v,		///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid_planedata(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);



	void run_timestep_agrid_planedatacomplex(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

	void setup(
			const std::string &i_function_name = "phi0"
	);

	virtual ~SWE_Plane_TS_l_direct();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_DIRECT_HPP_ */
