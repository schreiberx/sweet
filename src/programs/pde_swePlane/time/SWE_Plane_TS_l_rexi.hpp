/*
 * SWE_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TS_L_REXI_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TS_L_REXI_HPP_

#include <sweet/core/defaultPrecompilerValues.hpp>
#include <limits>
#include <string>
#include <complex>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>

#include "PDESWEPlaneTS_BaseInterface.hpp"
#include "SWE_Plane_TS_l_direct.hpp"


#if SWEET_BENCHMARK_TIMINGS
#	include <sweet/core/Stopwatch.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

class SWE_Plane_TS_l_rexi	:
		public PDESWEPlaneTS_BaseInterface
{
	std::vector<std::complex<double>> rexi_alphas;
	std::vector<std::complex<double>> rexi_betas;
	std::complex<double> rexi_gamma;

	/// simulation domain size
	double domain_size[2];

	std::size_t block_size;

#if SWEET_BENCHMARK_TIMINGS
	sweet::Stopwatch stopwatch_preprocessing;
	sweet::Stopwatch stopwatch_broadcast;
	sweet::Stopwatch stopwatch_reduce;
	sweet::Stopwatch stopwatch_solve_rexi_terms;
#endif

	class PerThreadVars
	{
	public:
		sweet::PlaneOperatorsComplex ops;

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
	bool setup(
			sweet::PlaneOperators *io_ops,
			const std::string &i_function_name
	);

	bool setup(
			sweet::PlaneOperators *io_ops
	);

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
