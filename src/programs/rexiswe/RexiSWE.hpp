/*
 * rexi_swe.hpp
 *
 *  Created on: 24 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_PROGRAMS_REXISWE_HPP_
#define SRC_PROGRAMS_REXISWE_HPP_


#include <complex>
#include <rexi/REXI.hpp>
#include <sweet/DataArray.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/Complex2DArrayFFT.hpp>
#include <sweet/Sampler2D.hpp>
#include <sweet/SemiLagrangian.hpp>
#include "../rexiswe/RexiSWE_HelmholtzSolver.hpp"
//#include <sweet/SWEValidationBenchmarks.hpp>

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
class RexiSWE
{
	double h;
	int M;

	// ID of solver to use, specify -1 to print a list of solvers
	int helmholtz_solver;

	// simulation domain size
	double domain_size[2];

	// error threshold, e.g. for iterative solvers
	double eps;

	// use finite differences
	bool use_spec_diffs;

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
		Complex2DArrayFFT op_diff_c_x, op_diff_c_y;
		Complex2DArrayFFT op_diff2_c_x, op_diff2_c_y;

		Complex2DArrayFFT eta;

		Complex2DArrayFFT eta0;
		Complex2DArrayFFT u0;
		Complex2DArrayFFT v0;

		Complex2DArrayFFT h_sum;
		Complex2DArrayFFT u_sum;
		Complex2DArrayFFT v_sum;
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

	// Interpolation stuff
	//Sampler2D sampler2D;

	// Semi-Lag stuff
	//SemiLagrangian semiLagrangian;

private:
	void cleanup();

public:
	RexiSWE();


	/**
	 * setup the REXI
	 */
public:
	void setup(
			double i_h,		///< sampling size
			int i_M,		///< number of sampling points
			int i_L,		///< number of sampling points for Gaussian approx

			std::size_t *i_resolution,			///< resolution of domain
			const double *i_domain_size,		///< size of domain
			bool i_rexi_half = true,			///< use half-pole reduction
			bool i_use_finite_differences = false,		///< use finite-differences for derivatives,	///< use finite differences for REXI approximation
			int i_helmholtz_solver = 0,		///< Helmholtz solver ID to use (0: spectral space)
			double i_eps = 1e-7		///< error threshold
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
			const Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			int i_thread_id = 0
	)
	{
		// compute
		// 		kappa - g * eta_bar * D2
		// NOTE!!! We add kappa in Cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
		// This is *NOT* straightforward and different to adding a constant for computations.
		// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.

		Complex2DArrayFFT lhs = ((-i_gh0)*(perThreadVars[i_thread_id]->op_diff2_c_x + perThreadVars[i_thread_id]->op_diff2_c_y)).addScalar_Cart(i_kappa);

		io_x = ((i_rhs.toSpec()).spec_div_element_wise(lhs)).toCart();
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
			const Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			int i_thread_id = 0
	)
	{
		// compute
		// 		kappa - g * eta_bar * D2
		// NOTE!!! We add kappa in Cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
		// This is *NOT* straightforward and different to adding a constant for computations.
		// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.

		Complex2DArrayFFT lhs = (-i_gh0*(perThreadVars[i_thread_id]->op_diff2_c_x + perThreadVars[i_thread_id]->op_diff2_c_y)).addScalar_Cart(i_kappa);
		//std::cout << lhs << std::endl;
		io_x = i_rhs.spec_div_element_wise(lhs);
	}

	/**
	 * Solve real-valued Helmholtz problem with a spectral solver,
	 * values are given in Spectral space
	 * ***under development! Bugged!! ************
	 * (kappa - gh*D2) X = B
	 *
	 * This is here to compare speedups and such a solver cannot be applied in general,
	 * e.g. with special general boundary values
	 */
public:
	void helmholtz_spectral_solver(
			double i_kappa,
			double i_gh0,
			const DataArray<2> &i_rhs,
			DataArray<2> &io_x,
			Operators2D &op     ///< Operator class
	)
	{

		// compute
		// 		kappa - g * eta_bar * D2
		// NOTE!!! We add kappa in Cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
		// This is *NOT* straightforward and different to adding a constant for computations.
		// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.

		//std::cout<<"gh0: " << i_gh0 <<std::endl;
		//std::cout<<"kappa: " << i_kappa<<std::endl;
		//This works
#if SWEET_USE_SPECTRAL_SPACE
		DataArray<2> laplacian = -i_gh0*op.diff2_c_x -i_gh0*op.diff2_c_y;
		DataArray<2> lhs = laplacian.spec_addScalarAll(i_kappa);
		std::cout << "Correct lhs" << std::endl;
		std::cout << lhs << std::endl;
		std::cout << "Correct lhs spec" << std::endl;
		lhs.printSpectrum();

		// This does not work
		DataArray<2> lhs2 = -i_gh0*(op.diff2_c_x +op.diff2_c_y).spec_addScalarAll(i_kappa);
		std::cout << "Nonsense lhs" << std::endl;
		std::cout << lhs2 << std::endl;
		std::cout << std::endl;
		std::cout << "Nonsense lhs spec" << std::endl;
		lhs2.printSpectrum();
		std::cout << "Diff" << std::endl;
		std::cout << (lhs-lhs2).reduce_maxAbs() << std::endl;
		io_x = i_rhs.spec_div_element_wise(lhs);
#else
		std::cerr << "Cannot use helmholtz_spectral_solver if spectral space not enable in compilation time" << std::endl;
		exit(1);
#endif
		exit(1);
	}


	/**
	 * Solve U_t = L U via implicit solver:
	 *
	 *
	 */
public:
	bool run_timestep_implicit_ts(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		double i_timestep_size,	///< timestep size

		Operators2D &op,
		const SimulationVariables &i_parameters
	);


	/**
	 * Solve U_t = L U via Crank-Nicolson:
	 * with (semi)-implicit semi-lagrangian solver
	 */
public:
	bool run_timestep_cn_sl_ts(
			DataArray<2> &io_h,  ///< Current and past fields
			DataArray<2> &io_u,
			DataArray<2> &io_v,
			DataArray<2> &io_h_prev,
			DataArray<2> &io_u_prev,
			DataArray<2> &io_v_prev,

			DataArray<2> &i_posx_a, //Arrival point positions in x and y (this is basically the grid)
			DataArray<2> &i_posy_a,

			double i_timestep_size,	///< timestep size
			int i_param_nonlinear, ///< degree of nonlinearity (0-linear, 1-full nonlinear, 2-only nonlinear adv)

			const SimulationVariables &i_simVars, ///< Parameters for simulation

			Operators2D &op,     ///< Operator class
			Sampler2D &sampler2D, ///< Interpolation class
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
	bool run_timestep(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		double i_timestep_size,	///< timestep size

		Operators2D &op,
		const SimulationVariables &i_parameters,
		bool i_iterative_solver_zero_solution = false
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
			DataArray<2> &io_h,
			DataArray<2> &io_u,
			DataArray<2> &io_v,

			double i_timestep_size,	///< timestep size

			Operators2D &op,
			const SimulationVariables &i_parameters
	);


public:
	inline
	static
	void MPI_quitWorkers(
			std::size_t i_resolution[2]
	)
	{
#if SWEET_MPI
	DataArray<2> dummyData(i_resolution);
	dummyData.set_all(NAN);

	MPI_Bcast(dummyData.array_data_cartesian_space, dummyData.resolution[0]*dummyData.resolution[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif
	}

	~RexiSWE();
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
