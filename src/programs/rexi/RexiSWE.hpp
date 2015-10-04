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
#include "RexiSWE_HelmholtzSolver.hpp"



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
	double tau;
	double h;
	int M;
	double f;

	// ID of solver to use, specify -1 to print a list of solvers
	int helmholtz_solver;

	// simulation domain size
	double domain_size[2];

	// error threshold, e.g. for iterative solvers
	double eps;

	// use finite differences
	bool use_finite_differences;

	std::size_t block_size;

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
	int num_threads;


public:
	REXI rexi;

private:
	void cleanup();

public:
	RexiSWE();


	/**
	 * setup the REXI
	 */
public:
	void setup(
			double i_tau,	///< time step size
			double i_h,		///< sampling size
			int i_M,		///< number of sampling points
			int i_L,		///< number of sampling points for Gaussian approx
			double i_f,		///< Coriolis force
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

		Complex2DArrayFFT lhs = (-i_gh0*(perThreadVars[i_thread_id]->op_diff2_c_x + perThreadVars[i_thread_id]->op_diff2_c_y)).addScalar_Cart(i_kappa);

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

		io_x = i_rhs.spec_div_element_wise(lhs);
	}



	/**
	 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
	 *
	 * See
	 * 		doc/rexi/understanding_rexi.pdf
	 * for further information
	 */
	void run_timestep(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		Operators2D &op,
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
	void run_timestep_direct_solution(
			DataArray<2> &io_h,
			DataArray<2> &io_u,
			DataArray<2> &io_v,

			double i_timestep_size,

			Operators2D &op,
			const SimulationVariables &i_parameters
	);

	~RexiSWE();
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
