/*
 * SWE_Sphere_REXI.hpp
 *
 *  Created on: 25 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_REXI_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_REXI_HPP_


#ifndef SWEET_MPI
#define SWEET_MPI 1
#endif


#include <complex>
#include <rexi/REXITerry.hpp>
#include <sweet/SimulationVariables.hpp>
#include <string.h>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/app_swe/SWERexiTerm_SPHRobert.hpp>

#include "SWE_Sphere_TS_interface.hpp"


#if SWEET_MPI
	#include <mpi.h>
#endif



class SWE_Sphere_TS_l_rexi	: public SWE_Sphere_TS_interface
{
	typedef std::complex<double> complex;


	/// Simulation variables
	SimulationVariables &simVars;

	/// Sphere operators
	SphereOperators &op;


	double h;
	int M;
	bool normalization;

	SphereDataConfig *sphereDataConfig;
	SphereDataConfig *sphereDataConfigSolver;

	/// This class is only setp and used in case of added modes
	SphereDataConfig sphereDataConfigInstance;

	/// Extend modes for REXI time stepping?
	int rexi_use_extended_modes;

#if SWEET_MPI
public:
	bool final_timestep;

#endif

private:


	/*
	 * Time step size of REXI
	 */
	double timestep_size;

	bool use_f_sphere;

	bool use_rexi_preallocation;


	std::size_t block_size;

#if SWEET_BENCHMARK_REXI
	Stopwatch stopwatch_preprocessing;
	Stopwatch stopwatch_broadcast;
	Stopwatch stopwatch_reduce;
	Stopwatch stopwatch_solve_rexi_terms;
#endif

	class PerThreadVars
	{
	public:
		std::vector<SWERexiTerm_SPHRobert> rexiSPHRobert_vector;

		std::vector< std::complex<double> > alpha;
		std::vector< std::complex<double> > beta_re;

		SphereData accum_phi;
		SphereData accum_vort;
		SphereData accum_div;
	};

	// per-thread allocated variables to avoid NUMA domain effects
	std::vector<PerThreadVars*> perThreadVars;

	// number of threads to be used
	int num_local_rexi_par_threads;

	// number of mpi ranks to be used
	int mpi_rank;

	// MPI ranks
	int num_mpi_ranks;

	// number of threads to be used
	int num_global_threads;


public:
	// REXI stuff
	REXITerry<> rexi;


private:
	void reset();


public:
	SWE_Sphere_TS_l_rexi(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
		);

private:
	void update_coefficients();

	/**
	 * setup the REXI
	 */
public:
	void setup(
			double i_h,		///< sampling size
			int i_M,		///< number of sampling points
			int i_L,		///< number of sampling points for Gaussian approx

			double i_timestep_size,
			bool i_rexi_half,				///< use half-pole reduction
			int i_rexi_use_extended_modes,
			int i_rexi_normalization,

			bool i_use_f_sphere,

			bool i_use_rexi_sphere_preallocation
	);

	void run_timestep(
			SphereData &io_h,	///< prognostic variables
			SphereData &io_u,	///< prognostic variables
			SphereData &io_v,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	);

	void get_workload_start_end(
			std::size_t &o_start,
			std::size_t &o_end
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
		SphereData &io_phi0,
		SphereData &io_u0,
		SphereData &io_v0,

		double i_timestep_size,	///< timestep size

		const SimulationVariables &i_parameters
	);



public:
	static
	void MPI_quitWorkers(
			SphereDataConfig *i_sphereDataConfig
	);

	virtual ~SWE_Sphere_TS_l_rexi();
};


#endif /* SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_REXI_HPP_ */
