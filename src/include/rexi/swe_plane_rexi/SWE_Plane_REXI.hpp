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



#if 0
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
		/*
		 * compute
		 * 		kappa - g * eta_bar * D2
		 *
		 * NOTE!!! We add kappa in Cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
		 * This is *NOT* straightforward and different to adding a constant for computations.
		 * We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.
		 */
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
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PlaneData laplacian = -i_gh0*op.diff2_c_x -i_gh0*op.diff2_c_y;
		PlaneData lhs = laplacian.spectral_addScalarAll(i_kappa);

		io_x = i_rhs.spectral_div_element_wise(lhs);
#else
		FatalError("Cannot use helmholtz_spectral_solver if spectral space not enable in compilation time");
#endif
	}
#endif


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
			PlaneDataSemiLagrangian &semiLagrangian  ///< Semi-Lag class
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
		PlaneDataSemiLagrangian &semiLagrangian  ///< Semi-Lag class
	);


};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
