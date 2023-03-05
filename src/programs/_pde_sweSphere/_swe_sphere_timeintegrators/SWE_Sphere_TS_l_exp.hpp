/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_EXP_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_EXP_HPP_


#ifndef SWEET_MPI
#define SWEET_MPI 1
#endif


#include <complex>
#include <rexi/REXI_Terry.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <string.h>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperatorsComplex.hpp>
#include <sweet/expIntegration/ExpFunctions.hpp>
#include "helpers/SWERexiTerm_SPH.hpp"
#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_exp_direct_special.hpp"
#include "SWE_Sphere_TS_lg_exp_lc_taylor.hpp"


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



class SWE_Sphere_TS_l_exp	: public SWE_Sphere_TS_interface
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
	sweet::ShackDictionary &shackDict;

public:
	// WARNING: Do NOT use a reference to this to get more flexibility by overriding certain things in here
	SimulationCoefficients simCoeffs;

private:
	/// Sphere operators
	sweet::SphereOperators &ops;


public:
	std::vector<std::complex<double>> rexi_alphas;
	std::vector<std::complex<double>> rexi_betas;
	std::complex<double> rexi_gamma;


	const sweet::SphereDataConfig *sphereDataConfig;

	/// This class is only setup and used in case of added modes
	sweet::SphereDataConfig sphereDataConfigInstance;

	ExpFunctions<double> expFunctions;


private:

	EXP_sweet::ShackDictionary *rexiSimVars;

	/*
	 * Time step size of REXI
	 */
	double timestep_size;

	/*
	 * Function name to be used by REXI
	 */
	std::string function_name;

	/*
	 * Don't use any Coriolis effect (reduction to very simple Helmholtz problem)
	 */
	bool no_coriolis;

	/*
	 * Assume f-sphere (reduction to Helmholtz problem)
	 */
	bool use_f_sphere;

	/*
	 * Preallocate the REXI matrices
	 */
	bool use_rexi_sphere_solver_preallocation;

	/*
	 * True, if EXP method is "direct"
	 */
	bool use_exp_method_direct_solution;

	/*
	 * True, if EXP method is "ss_taylor"
	 */
	bool use_exp_method_strang_split_taylor;

	/*
	 * True, if a REXI method is used
	 */
	bool use_exp_method_rexi;


	std::size_t block_size;

	class PerThreadVars
	{
	public:
		std::vector<SWERexiTerm_SPH> rexiTermSolvers;

		std::vector< std::complex<double> > alpha;
		std::vector< std::complex<double> > beta;

		sweet::SphereData_Spectral accum_phi;
		sweet::SphereData_Spectral accum_vrt;
		sweet::SphereData_Spectral accum_div;
	};

	// per-thread allocated variables to avoid NUMA domain effects
	std::vector<PerThreadVars*> perThreadVars;

	// number of threads to be used
	int num_local_rexi_par_threads;

	// number of threads to be used
	int num_global_threads;


#if SWEET_MPI
	// MPI communicator
	MPI_Comm mpi_comm;

	// number of mpi ranks to be used
	int mpi_rank;

	// MPI ranks
	int num_mpi_ranks;
#endif

	SWE_Sphere_TS_l_exp_direct_special *timestepping_method_l_exp_direct_special;

	SWE_Sphere_TS_lg_exp_lc_taylor *timestepping_method_lg_exp_lc_exp;


private:
	void reset();


public:
	SWE_Sphere_TS_l_exp(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

private:
	void p_update_coefficients();

	void p_get_workload_start_end(
			std::size_t &o_start,
			std::size_t &o_end,
			int i_local_thread_id
	);


	/**
	 * setup the REXI
	 */
public:
	void setup(
			EXP_sweet::ShackDictionary &i_rexi,
			const std::string &i_function_name,
			double i_timestep_size,
			bool i_use_f_sphere,
			bool i_no_coriolis,
			int i_timestepping_order
	);


	void run_timestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h,	///< prognostic variables
			sweet::SphereData_Spectral &o_u,	///< prognostic variables
			sweet::SphereData_Spectral &o_v,	///< prognostic variables

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


	virtual ~SWE_Sphere_TS_l_exp();
};


#endif
