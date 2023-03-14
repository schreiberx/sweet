/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDESWESphereTS_LG_EXP_DIRECT_HPP_
#define SRC_PROGRAMS_PDESWESphereTS_LG_EXP_DIRECT_HPP_



#include <sweet/expIntegration/ExpFunctions.hpp>
#include <complex>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <string.h>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperatorsComplex.hpp>
#include "helpers/SWERexiTerm_SPH.hpp"

#include "PDESWESphereTS_BaseInterface.hpp"


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



class PDESWESphereTS_lg_exp_direct	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

public:
	bool setup_main(
			sweet::SphereOperators *io_ops,
			const std::string &i_function_name
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;

	std::string getIDString() override;


private:
	typedef std::complex<double> complex;

private:

	const sweet::SphereData_Config *sphereDataConfig;

	sweet::ExpFunctions<double> expFunctions;


private:
	/*
	 * Function name to be used by REXI
	 */
	std::string function_name;


private:
	void reset();


public:

	void run_timestep_lg_exp(
		sweet::SphereData_Spectral &io_prog_phi,
		sweet::SphereData_Spectral &io_prog_vrt,
		sweet::SphereData_Spectral &io_prog_div,

		double i_dt,
		double i_simulation_timestamp
	);

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vort,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	void runTimestep(
			const sweet::SphereData_Spectral &i_h,
			const sweet::SphereData_Spectral &i_u,
			const sweet::SphereData_Spectral &i_v,

			sweet::SphereData_Spectral &o_h,
			sweet::SphereData_Spectral &o_u,
			sweet::SphereData_Spectral &o_v,

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
		sweet::SphereData_Spectral &io_phi0,
		sweet::SphereData_Spectral &io_u0,
		sweet::SphereData_Spectral &io_v0,

		double i_timestep_size,	///< timestep size

		const sweet::ShackDictionary &i_parameters
	);

	PDESWESphereTS_lg_exp_direct();

	virtual ~PDESWESphereTS_lg_exp_direct();
};


#endif
