/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_Sphere_TS_LG_EXP_DIRECT_HPP_
#define SRC_PROGRAMS_SWE_Sphere_TS_LG_EXP_DIRECT_HPP_



#include <sweet/expIntegration/ExpFunctions.hpp>
#include <complex>
#include <sweet/core/SimulationVariables.hpp>
#include <string.h>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>
#include <sweet/core/sphere/SphereOperators_SphereDataComplex.hpp>
#include "helpers/SWERexiTerm_SPH.hpp"

#include "SWE_Sphere_TS_interface.hpp"


#ifndef SWEET_BENCHMARK_TIMINGS
	#define SWEET_BENCHMARK_TIMINGS 1
#endif

#ifndef SWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS
	#define SWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS 1
#endif



#if SWEET_BENCHMARK_TIMINGS
#include <sweet/core/Stopwatch.hpp>
#endif

#if SWEET_MPI
	#include <mpi.h>
#endif



class SWE_Sphere_TS_lg_exp_direct	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(
			const std::string &i_timestepping_method
	);

	std::string string_id();
	void setup_auto();


private:
	typedef std::complex<double> complex;


	/// Simulation variables
	SimulationVariables &simVars;

public:
	// WARNING: Do NOT use a reference to this to get more flexibility by overriding certain things in here
	SimulationCoefficients simCoeffs;

private:
	/// Sphere operators
	SphereOperators_SphereData &op;


	const sweet::SphereData_Config *sphereDataConfig;

	ExpFunctions<double> rexiFunctions;


private:
	/*
	 * Function name to be used by REXI
	 */
	std::string function_name;


private:
	void reset();


public:
	SWE_Sphere_TS_lg_exp_direct(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	/**
	 * setup the REXI
	 */
public:
	void setup(
			const std::string &i_function_name
	);

	void run_timestep_lg_exp(
		SphereData_Spectral &io_prog_phi,
		SphereData_Spectral &io_prog_vrt,
		SphereData_Spectral &io_prog_div,

		double i_dt,
		double i_simulation_timestamp
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			SphereData_Spectral &o_h,	///< prognostic variables
			SphereData_Spectral &o_u,	///< prognostic variables
			SphereData_Spectral &o_v,	///< prognostic variables

			double i_fixed_dt,
			double i_simulation_timestamp
	);


	/**
	 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
	 *
	 * See
	 * 		doc/rexi/understanding_rexi.pdf
	 * for further information
	 */
public:
	bool run_timestep_rexi_vectorinvariant_progphivortdiv(
		SphereData_Spectral &io_phi0,
		SphereData_Spectral &io_u0,
		SphereData_Spectral &io_v0,

		double i_timestep_size,	///< timestep size

		const SimulationVariables &i_parameters
	);


	virtual ~SWE_Sphere_TS_lg_exp_direct();
};


#endif
