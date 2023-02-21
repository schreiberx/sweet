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
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#include "../swe_plane_timeintegrators/SWE_Plane_TS_interface.hpp"
#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_direct.hpp"


#ifndef SWEET_BENCHMARK_TIMINGS
	#define SWEET_BENCHMARK_TIMINGS	0
#endif

#if SWEET_BENCHMARK_TIMINGS
#	include <sweet/Stopwatch.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

class SWE_Plane_TS_l_rexi	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	EXP_SimulationVariables *rexiSimVars;

	PlaneOperators &op;

	std::vector<std::complex<double>> rexi_alphas;
	std::vector<std::complex<double>> rexi_betas;
	std::complex<double> rexi_gamma;

	/// simulation domain size
	double domain_size[2];

	std::size_t block_size;

	PlaneDataConfig *planeDataConfig;

#if SWEET_BENCHMARK_TIMINGS
	Stopwatch stopwatch_preprocessing;
	Stopwatch stopwatch_broadcast;
	Stopwatch stopwatch_reduce;
	Stopwatch stopwatch_solve_rexi_terms;
#endif

	class PerThreadVars
	{
	public:
		PlaneOperatorsComplex op;

		PlaneData_SpectralComplex eta;

		PlaneData_SpectralComplex eta0;
		PlaneData_SpectralComplex u0;
		PlaneData_SpectralComplex v0;

		PlaneData_SpectralComplex h_sum;
		PlaneData_SpectralComplex u_sum;
		PlaneData_SpectralComplex v_sum;
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
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
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
			const PlaneData_Spectral &i_h_pert,	///< prognostic variables
			const PlaneData_Spectral &i_u,	///< prognostic variables
			const PlaneData_Spectral &i_v,	///< prognostic variables

			PlaneData_Spectral &o_h_pert,	///< prognostic variables
			PlaneData_Spectral &o_u,	///< prognostic variables
			PlaneData_Spectral &o_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_real(
			const PlaneData_Spectral &i_h_pert,	///< prognostic variables
			const PlaneData_Spectral &i_u,	///< prognostic variables
			const PlaneData_Spectral &i_v,	///< prognostic variables

			PlaneData_Spectral &o_h_pert,	///< prognostic variables
			PlaneData_Spectral &o_u,	///< prognostic variables
			PlaneData_Spectral &o_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

	void run_timestep(
			PlaneData_Spectral &io_h_pert,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	void cleanup();



public:
	static
	void MPI_quitWorkers(
			PlaneDataConfig *i_planeDataConfig
	);



	virtual ~SWE_Plane_TS_l_rexi();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
