/*
 * SWE_Sphere_REXI.hpp
 *
 *  Created on: 25 Oct 2016
 *      Author: martin
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_REXI_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_REXI_HPP_

#include <complex>
#include <rexi/REXI.hpp>
#include <sweet/SimulationVariables.hpp>
#include <string.h>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <rexi/SWERexi_SPH.hpp>
#include <rexi/SWERexi_SPHRobert.hpp>


class SWE_Sphere_REXI
{
	typedef std::complex<double> complex;

	double h;
	int M;

	SphereDataConfig *sphereDataConfig;
	SphereDataConfig *sphereDataConfigRexi;

	/*
	 * Simulation coefficients
	 */
	SimulationVariables::Coefficients *simCoeffs;

	bool use_robert_functions;


	/*
	 * Extend modes for REXI time stepping?
	 */
	int rexi_use_extended_modes;


	/*
	 * This class is only used in case of added modes
	 */
	SphereDataConfig sphereConfigRexiAddedModes;


	/*
	 * Time step size of REXI
	 */
	double timestep_size;

	bool use_coriolis_rexi_formulation;

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
		std::vector<SWERexi_SPHRobert> rexiSPHRobert_vector;
		std::vector<SWERexi_SPH> rexiSPH_vector;

		std::vector< std::complex<double> > alpha;
		std::vector< std::complex<double> > beta_re;

		SphereData accum_phi;
		SphereData accum_u;
		SphereData accum_v;

		SphereDataConfig sphereDataConfig;
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
	//REXI stuff
	REXI rexi;


private:
	void cleanup();

public:
	SWE_Sphere_REXI();


	/**
	 * setup the REXI
	 */
public:
	void setup(
			double i_h,		///< sampling size
			int i_M,		///< number of sampling points
			int i_L,		///< number of sampling points for Gaussian approx

			SphereDataConfig *i_sphereDataConfig,
			SimulationVariables::Coefficients *i_simCoeffs,
			double i_timestep_size,

			bool i_rexi_half,				///< use half-pole reduction
			bool i_use_robert_functions,	///< use Robert functions
			int i_rexi_use_extended_modes,
			bool i_use_coriolis_rexi_formulation
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
	bool run_timestep_rexi(
		SphereData &io_h,
		SphereData &io_u,
		SphereData &io_v,

		double i_timestep_size,	///< timestep size

		const SimulationVariables &i_parameters
	);



public:
	inline
	static
	void MPI_quitWorkers(
			SphereDataConfig *i_sphereDataConfig
	)
	{
#if SWEET_MPI
	SphereData dummyData(i_sphereDataConfig);
	dummyData.set_all(NAN);

	MPI_Bcast(dummyData.physical_space_data, dummyData.sphereDataConfig->physical_array_data_number_of_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif
	}

	~SWE_Sphere_REXI();
};


#endif /* SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_REXI_HPP_ */
