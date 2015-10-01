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

	bool use_iterative_solver;

	double domain_size[2];

	std::size_t block_size;

	class PerThreadVars
	{
	public:
		Complex2DArrayFFT op_diff_c_x, op_diff_c_y;
		Complex2DArrayFFT op_diff2_c_x, op_diff2_c_y;

		Complex2DArrayFFT eta0;
		Complex2DArrayFFT u0;
		Complex2DArrayFFT v0;

		Complex2DArrayFFT h_sum;
		Complex2DArrayFFT u_sum;
		Complex2DArrayFFT v_sum;
	};

	std::vector<PerThreadVars*> perThreadVars;

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
			bool i_use_iterative_solver = false		///< Use iterative solver instead of direct solving it in spectral space
	);



	/**
	 * Solve complex-valued Helmholtz problem
	 *
	 * (kappa - gh*D2) X = B
	 */
	bool helmholtz_iterative_solve_X(
			std::complex<double> i_kappa,
			double i_gh0,
			Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double i_error_threshold = 0.000001,
			int i_max_iters = 999999999,
			int i_thread_id = 0,
			double i_omega = 0.1
	)
	{
		double helm_h[2];
		helm_h[0] = (double)domain_size[0] / (double)i_rhs.resolution[0];
		helm_h[1] = (double)domain_size[1] / (double)i_rhs.resolution[1];

		if (perThreadVars.size() == 0)
		{
			std::cerr << "RexiSWE: Setup not executed!" << std::endl;
			exit(-1);
		}
		PerThreadVars *t = perThreadVars[i_thread_id];

		std::complex<double> D2_kernel[3*3];
		{
			D2_kernel[0] = { 0, 0 };
			D2_kernel[1] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[2] = { 0, 0 };

			D2_kernel[3] = { 1.0/(helm_h[0]*helm_h[0]), 0 };
			D2_kernel[4] = { -(2.0/(helm_h[0]*helm_h[0]) + 2.0/(helm_h[1]*helm_h[1])), 0};
			D2_kernel[5] = { 1.0/(helm_h[0]*helm_h[0]), 0 };

			D2_kernel[6] = { 0, 0 };
			D2_kernel[7] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[8] = { 0, 0 };
		}

		double check = (t->op_diff2_c_x(i_rhs.toSpec()).toCart() + t->op_diff2_c_y(i_rhs.toSpec()).toCart() - i_rhs.op_stencil_3x3(D2_kernel)).reduce_rms();
		std::cout << check << std::endl;

		std::complex<double> iD2_kernel[3*3];
		{
			iD2_kernel[0] = D2_kernel[0];
			iD2_kernel[1] = D2_kernel[1];
			iD2_kernel[2] = D2_kernel[2];

			iD2_kernel[3] = D2_kernel[3];
			iD2_kernel[4] = {0,0};
			iD2_kernel[5] = D2_kernel[5];

			iD2_kernel[6] = D2_kernel[6];
			iD2_kernel[7] = D2_kernel[7];
			iD2_kernel[8] = D2_kernel[8];
		}


#if 0
		{
			double residual_pre = (-i_gh0*(t->op_diff2_c_x(io_x.toSpec()).toCart()+t->op_diff2_c_y(io_x.toSpec()).toCart()).addScalar_Cart(i_kappa) - i_rhs).reduce_rms();
			Complex2DArrayFFT lhs = (-i_gh0*(t->op_diff2_c_x + t->op_diff2_c_y)).addScalar_Cart(i_kappa);
			Complex2DArrayFFT eta = i_rhs.spec_div_element_wise(lhs);
			double error_pre = (io_x-eta).reduce_rms();
//			std::cout << "Iter START, residual_pre: " << residual_pre << ", error_pre: " << error_pre << std::endl;
		}
#endif

		std::complex<double> inv_diag = 1.0/(i_kappa - i_gh0*D2_kernel[4]);

		Complex2DArrayFFT new_x(i_rhs.resolution);
		int i = 0;
		double residual = 0;
		for (i = 0; i < i_max_iters; i++)
		{
			new_x = (i_rhs + i_gh0*io_x.op_stencil_3x3(iD2_kernel))*inv_diag;
			io_x = (1.0-i_omega)*io_x + i_omega*new_x;

			residual = ((io_x*i_kappa - i_gh0*io_x.op_stencil_3x3(D2_kernel)) - i_rhs).reduce_rms();

//			std::cout << residual << std::endl;

//			double residual = (-i_gh0*(t->op_diff2_c_x(io_x.toSpec()).toCart()+t->op_diff2_c_y(io_x.toSpec()).toCart()).addScalar_Cart(i_kappa) - i_rhs).reduce_rms();
#if 0
			Complex2DArrayFFT lhs = (-i_gh0*(t->op_diff2_c_x + t->op_diff2_c_y)).addScalar_Cart(i_kappa);
			Complex2DArrayFFT eta = i_rhs.toSpec().spec_div_element_wise(lhs).toCart();
			double error = (io_x-eta).reduce_rms();
			std::cout << "Iter " << i << ", residual: " << residual << ", error: " << error << std::endl;
#endif

			if (residual < i_error_threshold)
			{
				std::cout << "Required iterations: " << i << " with residual " << residual << std::endl;
				return true;
			}
		}

		std::cout << "ERROR, iterations exceeded: " << i << " with residual " << residual << std::endl;
		return false;
	}



	/**
	 * Solve complex-valued Helmholtz problem
	 *
	 * (kappa - gh*D2) X = B
	 */
	static
	bool helmholtz_iterative_smoother(
			std::complex<double> i_kappa,
			double i_gh0,
			Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double *i_domain_size,
			double i_error_threshold = 0.000001,
			int i_max_iters = 999999999,
			double i_omega = 0.1
	)
	{
		double helm_h[2];
		helm_h[0] = (double)i_domain_size[0] / (double)i_rhs.resolution[0];
		helm_h[1] = (double)i_domain_size[1] / (double)i_rhs.resolution[1];

		double scalar_Dx = 1.0/(helm_h[0]*helm_h[0]);
		double scalar_Dy = 1.0/(helm_h[1]*helm_h[1]);
		double center_D = -(2.0/(helm_h[0]*helm_h[0]) + 2.0/(helm_h[1]*helm_h[1]));

		std::complex<double> D2_kernel[3*3];
		{
			D2_kernel[0] = { 0, 0 };
			D2_kernel[1] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[2] = { 0, 0 };

			D2_kernel[3] = { 1.0/(helm_h[0]*helm_h[0]), 0 };
			D2_kernel[4] = { -(2.0/(helm_h[0]*helm_h[0]) + 2.0/(helm_h[1]*helm_h[1])), 0};
			D2_kernel[5] = { 1.0/(helm_h[0]*helm_h[0]), 0 };

			D2_kernel[6] = { 0, 0 };
			D2_kernel[7] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[8] = { 0, 0 };
		}

		std::complex<double> inv_diag = 1.0/(i_kappa - i_gh0*center_D);

		Complex2DArrayFFT new_x(i_rhs.resolution);
		int i = 0;
		double residual = 0;
		for (i = 0; i < i_max_iters; i++)
		{
			// do more iterations to compensate overheads of residual computation
			int j = i+9;
			for (; i < j; i++)
				io_x = (1.0-i_omega)*io_x + i_omega*((i_rhs + i_gh0*io_x.op_stencil_Re_X(scalar_Dx, scalar_Dy))*inv_diag);
			io_x = (1.0-i_omega)*io_x + i_omega*((i_rhs + i_gh0*io_x.op_stencil_Re_X(scalar_Dx, scalar_Dy))*inv_diag);

			residual = ((io_x*i_kappa - i_gh0*io_x.op_stencil_3x3(D2_kernel)) - i_rhs).reduce_rms();

			if (residual < i_error_threshold)
				return true;
		}

		return false;
	}


	static
	bool helmholtz_iterative_solve_mg(
			std::complex<double> i_kappa,
			double i_gh0,
			Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double *i_domain_size,
			double i_error_threshold = 0.000001,
			int i_max_iters = 999999999,
			double i_omega = 0.1,
			int level_restriction = 0
	)
	{
		int levels = 1;
		while ((std::size_t)(1 << levels) < std::min(i_rhs.resolution[0], i_rhs.resolution[1]))
			levels++;

		std::vector<Complex2DArrayFFT> rhs;
		std::vector<Complex2DArrayFFT> x;
		rhs.resize(levels);
		x.resize(levels);

		if (level_restriction < 0)
			level_restriction = levels + level_restriction-1;

//		std::cout << "USING " << levels << " levels" << std::endl;

		rhs[0] = i_rhs;
		x[0] = io_x;

		for (int l = 0; l < levels-1-level_restriction; l++)
		{
//			std::cout << "DOWN -> LEVEL: " << l << std::endl;

			// presmoother
#if 0
			helmholtz_iterative_smoother(
					i_kappa,
					i_gh0,
					rhs[l],		// rhs
					x[l],		// x
					i_error_threshold,
					3,
					i_omega
			);
#endif
//			std::cout << "      downscale to level " << l+1 << std::endl;

			rhs[l+1] = rhs[l].scale_down();
			x[l+1] = x[l].scale_down();
		}

		for (int l = levels-2-level_restriction; l >= 1; l--)
		{
//			std::cout << " upscale to level " << l-1 << std::endl;
			x[l-1] = x[l].scale_up();
			// reuse old RHS, so use this only for debugging purpose
//			rhs[l-1] = rhs[l].scale_up();

//			std::cout << "UP -> LEVEL: " << l << std::endl;

			bool retval = helmholtz_iterative_smoother(
					i_kappa,
					i_gh0,
					rhs[l-1],	// rhs
					x[l-1],		// x
					i_domain_size,
					i_error_threshold,
					i_max_iters,
					i_omega
			);

			if (!retval)
			{
				std::cout << "No convergence for iterative solver" << std::endl;
				return false;
//				exit(1);
			}
		}

		io_x = x[0];
		return true;
	}


	double helmholtz_iterative_get_residual(
			std::complex<double> i_kappa,
			double i_gh0,
			Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x
	)
	{
		double helm_h[2];
		helm_h[0] = (double)domain_size[0] / (double)i_rhs.resolution[0];
		helm_h[1] = (double)domain_size[1] / (double)i_rhs.resolution[1];

		std::complex<double> D2_kernel[3*3];
		{
			D2_kernel[0] = { 0, 0 };
			D2_kernel[1] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[2] = { 0, 0 };

			D2_kernel[3] = { 1.0/(helm_h[0]*helm_h[0]), 0 };
			D2_kernel[4] = { -(2.0/(helm_h[0]*helm_h[0]) + 2.0/(helm_h[1]*helm_h[1])), 0};
			D2_kernel[5] = { 1.0/(helm_h[0]*helm_h[0]), 0 };

			D2_kernel[6] = { 0, 0 };
			D2_kernel[7] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[8] = { 0, 0 };
		}

		return ((io_x*i_kappa - i_gh0*io_x.op_stencil_3x3(D2_kernel)) - i_rhs).reduce_rms();
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
