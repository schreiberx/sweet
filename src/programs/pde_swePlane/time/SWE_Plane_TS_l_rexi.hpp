/*
 * SWE_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TS_L_REXI_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TS_L_REXI_HPP_

#include <limits>
#include <string>
#include <complex>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>

#include "SWE_Plane_TS_interface.hpp"
#include "SWE_Plane_TS_l_direct.hpp"


#ifndef SWEET_BENCHMARK_TIMINGS
	#define SWEET_BENCHMARK_TIMINGS	0
#endif

#if SWEET_BENCHMARK_TIMINGS
#	include <sweet/core/Stopwatch.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

class SWE_Plane_TS_l_rexi	: public SWE_Plane_TS_interface
{
	sweet::ShackDictionary *shackDict;
	EXP_SimulationVariables *rexiSimVars;

	sweet::PlaneOperators &op;

	std::vector<std::complex<double>> rexi_alphas;
	std::vector<std::complex<double>> rexi_betas;
	std::complex<double> rexi_gamma;

	/// simulation domain size
	double domain_size[2];

	std::size_t block_size;

	sweet::PlaneDataConfig *planeDataConfig;

#if SWEET_BENCHMARK_TIMINGS
	Stopwatch stopwatch_preprocessing;
	Stopwatch stopwatch_broadcast;
	Stopwatch stopwatch_reduce;
	Stopwatch stopwatch_solve_rexi_terms;
#endif

	class PerThreadVars
	{
	public:
		sweet::PlaneOperatorsComplex op;

		sweet::PlaneData_SpectralComplex eta;

		sweet::PlaneData_SpectralComplex eta0;
		sweet::PlaneData_SpectralComplex u0;
		sweet::PlaneData_SpectralComplex v0;

		sweet::PlaneData_SpectralComplex h_sum;
		sweet::PlaneData_SpectralComplex u_sum;
		sweet::PlaneData_SpectralComplex v_sum;
	};

	/// per-thread allocated variables to avoid NUMA domain effects
	std::vector<PerThreadVars*> perThreadVars;

	/// number of threads to be used
	int num_local_rexi_par_threads;

	/// number of mpi ranks to be used
	int mpi_rank;

	/// MPI ranks
	int num_mpi_ranks;

	/// number of threads to be used
	int num_global_threads;

public:
	/// final time step
	bool final_timestep;

	/// use direct solution instead of REXI
	bool rexi_use_direct_solution;

	/// Direct solution for linear parts
	SWE_Plane_TS_l_direct ts_l_direct;

public:
	SWE_Plane_TS_l_rexi(
			sweet::ShackDictionary *shackDict,
			sweet::PlaneOperators &i_op
		);

	void setup(
			EXP_SimulationVariables &i_rexi,
			const std::string &i_function_name,
			double i_timestep_size
	);
/*
	void setup_REXI(
			double i_h,						///< sampling size
			int i_M,						///< number of sampling points
			int i_L,						///< number of sampling points for Gaussian approximation
											///< set to 0 for auto detection

			bool i_rexi_half,				///< use half-pole reduction
			bool i_rexi_normalization,		///< REXI normalization
			bool i_rexi_next_generation
	);
*/

	void run_timestep(
			const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
			sweet::PlaneData_Spectral &o_u,	///< prognostic variables
			sweet::PlaneData_Spectral &o_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_real(
			const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
			sweet::PlaneData_Spectral &o_u,	///< prognostic variables
			sweet::PlaneData_Spectral &o_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

	void run_timestep(
			sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	void cleanup();



public:
	static
	void MPI_quitWorkers(
			sweet::PlaneDataConfig *i_planeDataConfig
	);



	virtual ~SWE_Plane_TS_l_rexi();
};

#endif
