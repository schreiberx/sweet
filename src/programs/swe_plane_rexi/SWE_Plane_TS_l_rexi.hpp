/*
 * SWE_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TS_L_REXI_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TS_L_REXI_HPP_

#include <limits>
#include <complex>
#include <rexi/REXI.hpp>
#include <rexi/RexiNG.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include "SWE_Plane_TS_interface.hpp"
#include "SWE_Plane_TS_l_direct.hpp"


#define SWEET_BENCHMARK_REXI	0

#if SWEET_BENCHMARK_REXI
#	include <sweet/Stopwatch.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

class SWE_Plane_TS_l_rexi	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	/// REXI parameter h
	double h;

	/// REXI parameter M
	int M;

	/// simulation domain size
	double domain_size[2];

	std::size_t block_size;

	PlaneDataConfig *planeDataConfig;

#if SWEET_BENCHMARK_REXI
	Stopwatch stopwatch_preprocessing;
	Stopwatch stopwatch_broadcast;
	Stopwatch stopwatch_reduce;
	Stopwatch stopwatch_solve_rexi_terms;
#endif

	class PerThreadVars
	{
	public:
		PlaneOperatorsComplex op;

		PlaneDataComplex eta;

		PlaneDataComplex eta0;
		PlaneDataComplex u0;
		PlaneDataComplex v0;

		PlaneDataComplex h_sum;
		PlaneDataComplex u_sum;
		PlaneDataComplex v_sum;
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


public:
	/// REXI stuff
	REXI<> rexi;

public:
	/// REXI next generation stuff
	RexiNG<> rexiNG;

	/// Direct solution for linear parts
	SWE_Plane_TS_l_direct ts_l_direct;

public:
	SWE_Plane_TS_l_rexi(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			REXI_SimulationVariables &i_rexi
	);

	void setup(
			double i_h,						///< sampling size
			int i_M,						///< number of sampling points
			int i_L,						///< number of sampling points for Gaussian approximation
											///< set to 0 for auto detection

			bool i_rexi_half,				///< use half-pole reduction
			bool i_rexi_normalization,		///< REXI normalization
			bool i_rexi_next_generation
	);

	void run_timestep(
			PlaneData &io_h,	///< prognostic variables
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
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
