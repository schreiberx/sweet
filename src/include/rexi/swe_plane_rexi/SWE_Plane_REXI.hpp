/*
 * rexi_swe.hpp
 *
 *  Created on: 24 Jul 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_PROGRAMS_REXISWE_HPP_
#define SRC_PROGRAMS_REXISWE_HPP_


#include <complex>
#include <rexi/REXI.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>


#define SWEET_BENCHMARK_REXI	0

#if SWEET_BENCHMARK_REXI
#	include <sweet/Stopwatch.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

/**
 * This class implements the REXI (rational approximation of exponential integrator) solver for the SWE,
 * see High-order time-parallel approximation of evolution operators, T. Haut et al.
 *
 * We split this file into the header and cpp file.
 *
 * This allows using a OpenMP parallelization only in the RexiSWE class to check this degree of
 * parallelism.
 */
class SWE_Plane_REXI
{
	double h;
	int M;

	// ID of solver to use, specify -1 to print a list of solvers
//	int helmholtz_solver;

	// simulation domain size
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

//		PlaneDataComplex op_diff_c_x, op_diff_c_y;
//		PlaneDataComplex opc.diff2_c_x, opc.diff2_c_y;

		PlaneDataComplex eta;

		PlaneDataComplex eta0;
		PlaneDataComplex u0;
		PlaneDataComplex v0;

		PlaneDataComplex h_sum;
		PlaneDataComplex u_sum;
		PlaneDataComplex v_sum;
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
	REXI<> rexi;

	// Interpolation stuff
	//Sampler2D sampler2D;

	// Semi-Lag stuff
	//SemiLagrangian semiLagrangian;

private:
	void cleanup();

public:
	SWE_Plane_REXI();


	/**
	 * setup the REXI
	 */
public:
	void setup(
			double i_h,		///< sampling size
			int i_M,		///< number of sampling points
			int i_L,		///< number of sampling points for Gaussian approx

			PlaneDataConfig *i_planeDataConfig,
			const double *i_domain_size,		///< size of domain
			bool i_rexi_half = true,			///< use half-pole reduction
			bool i_rexi_normalization = true
	);



	/**
	 * Solve complex-valued Helmholtz problem with a spectral solver,
	 * values are given in Cartesian space
	 *
	 * (kappa - gh*D2) X = B
	 *
	 * This is here to compare speedups and such a solver cannot be applied in general,
	 * e.g. with special general boundary values
	 */
public:
	void helmholtz_spectral_solver_cart(
			std::complex<double> i_kappa,
			double i_gh0,
			const PlaneDataComplex &i_rhs,
			PlaneDataComplex &io_x,
			int i_thread_id = 0
	)
	{
		// compute
		// 		kappa - g * eta_bar * D2
		// NOTE!!! We add kappa in Cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
		// This is *NOT* straightforward and different to adding a constant for computations.
		// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.

		PlaneDataComplex lhs = ((-i_gh0)*(perThreadVars[i_thread_id]->op.diff2_c_x + perThreadVars[i_thread_id]->op.diff2_c_y)).spectral_addScalarAll(i_kappa);

		io_x =i_rhs.spectral_div_element_wise(lhs);
	}


	/**
	 * Solve complex-valued Helmholtz problem with a spectral solver,
	 * values are given in Spectral space
	 *
	 * (kappa - gh*D2) X = B
	 *
	 * This is here to compare speedups and such a solver cannot be applied in general,
	 * e.g. with special general boundary values
	 */
public:
	void helmholtz_spectral_solver_spec(
			std::complex<double> i_kappa,
			double i_gh0,
			const PlaneDataComplex &i_rhs,
			PlaneDataComplex &io_x,
			int i_thread_id = 0
	)
	{
		// compute
		// 		kappa - g * eta_bar * D2
		// NOTE!!! We add kappa in Cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
		// This is *NOT* straightforward and different to adding a constant for computations.
		// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.

		PlaneDataComplex lhs = (-i_gh0*(perThreadVars[i_thread_id]->op.diff2_c_x + perThreadVars[i_thread_id]->op.diff2_c_y)).spectral_addScalarAll(i_kappa);

		io_x = i_rhs.spectral_div_element_wise(lhs);
	}

	/**
	 * Solve real-valued Helmholtz problem with a spectral solver,
	 * values are given in Spectral space
	 *
	 * (kappa - gh*D2) X = B
	 *
	 * This is here to compare speedups and such a solver cannot be applied in general,
	 * e.g. with special general boundary values
	 */
public:
	void helmholtz_spectral_solver(
			double i_kappa,
			double i_gh0,
			const PlaneData &i_rhs,
			PlaneData &io_x,
			PlaneOperators &op     ///< Operator class
	)
	{

		// compute
		// 		kappa - g * eta_bar * D2
		// NOTE!!! We add kappa to all frequency components to account for scaling all frequencies!!!
		// This is *NOT* straightforward and different to adding a constant for computations.
		// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.


#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PlaneData laplacian = -i_gh0*op.diff2_c_x -i_gh0*op.diff2_c_y;
		PlaneData lhs = laplacian.spectral_addScalarAll(i_kappa);

		io_x = i_rhs.spectral_div_element_wise(lhs);
#else
		FatalError("Cannot use helmholtz_spectral_solver if spectral space not enable in compilation time");
#endif
	}


	/**
	 * Solve U_t = L U via implicit solver:
	 *
	 *
	 */
public:
	bool run_timestep_implicit_ts(
		PlaneData &io_h,
		PlaneData &io_u,
		PlaneData &io_v,

		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_parameters
	);


	/**
	 * Solve U_t = L U via Crank-Nicolson:
	 * with (semi)-implicit semi-lagrangian solver
	 */
public:
	bool run_timestep_cn_sl_ts(
			PlaneData &io_h,  ///< Current and past fields
			PlaneData &io_u,
			PlaneData &io_v,
			PlaneData &io_h_prev,
			PlaneData &io_u_prev,
			PlaneData &io_v_prev,

			ScalarDataArray &i_posx_a, //Arrival point positions in x and y (this is basically the grid)
			ScalarDataArray &i_posy_a,

			double i_timestep_size,	///< timestep size
			int i_param_nonlinear, ///< degree of nonlinearity (0-linear, 1-full nonlinear, 2-only nonlinear adv)

			const SimulationVariables &i_simVars, ///< Parameters for simulation

			PlaneOperators &op,     ///< Operator class
			PlaneDataSampler &sampler2D, ///< Interpolation class
			SemiLagrangian &semiLagrangian  ///< Semi-Lag class
	);

	/**
	 * Solve  SWE with the novel Semi-Lagrangian Exponential Integrator
	 *  SL-REXI
	 *
	 *  See documentation for details
	 *
	 */
public:
	bool run_timestep_slrexi(
		PlaneData &io_h,  ///< Current and past fields
		PlaneData &io_u,
		PlaneData &io_v,
		PlaneData &io_h_prev,
		PlaneData &io_u_prev,
		PlaneData &io_v_prev,

		ScalarDataArray &i_posx_a, //Arrival point positions in x and y (this is basically the grid)
		ScalarDataArray &i_posy_a,

		double i_timestep_size,	///< timestep size
		int i_param_nonlinear, ///< degree of nonlinearity (0-linear, 1-full nonlinear, 2-only nonlinear adv)

		bool i_linear_exp_analytical, //

		const SimulationVariables &i_simVars, ///< Parameters for simulation

		PlaneOperators &op,     ///< Operator class
		PlaneDataSampler &sampler2D, ///< Interpolation class
		SemiLagrangian &semiLagrangian  ///< Semi-Lag class
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
		PlaneData &io_h,
		PlaneData &io_u,
		PlaneData &io_v,

		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_parameters
	);


	/**
	 * This method computes the analytical solution based on the given initial values.
	 *
	 * See Embid/Madja/1996, Terry/Beth/2014, page 16
	 * and
	 * 		doc/swe_solution_for_L/sympy_L_spec_decomposition.py
	 * for the dimensionful formulation.
	 *
	 * Don't use this function to frequently, since it always computes
	 * the required coefficients on-the-fly which is expensive.
	 */
public:
	void run_timestep_direct_solution(
			PlaneData &io_h,
			PlaneData &io_u,
			PlaneData &io_v,

			double i_timestep_size,	///< timestep size

			PlaneOperators &op,
			const SimulationVariables &i_parameters
	);



public:
	void run_timestep_direct_solution_geopotential_formulation(
			PlaneData &io_phi,
			PlaneData &io_u,
			PlaneData &io_v,

			double i_timestep_size,	///< timestep size

			PlaneOperators &op,
			const SimulationVariables &i_parameters
	);



public:
	inline
	static
	void MPI_quitWorkers(
			PlaneDataConfig *i_planeDataConfig
	)
	{
#if SWEET_MPI
	PlaneData dummyData(i_planeDataConfig);
	dummyData.physical_set_all(NAN);

	MPI_Bcast(dummyData.physical_space_data, dummyData.planeDataConfig->physical_array_data_number_of_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif
	}

	~SWE_Plane_REXI();
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
