/*
 * rexi_swe.hpp
 *
 *  Created on: 24 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#include "RexiSWE.hpp"
#include <cmath>



#ifndef SWEET_REXI_THREAD_PARALLEL_SUM
#	define SWEET_REXI_THREAD_PARALLEL_SUM 1
#endif

/**
 * Compute the REXI sum massively parallel *without* a parallelization with parfor in space
 */
#if SWEET_REXI_THREAD_PARALLEL_SUM
#	include <omp.h>
#endif


#ifndef SWEET_MPI
#	define SWEET_MPI 1
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif


RexiSWE::RexiSWE()	:
	helmholtz_solver(0)
{
#if !SWEET_USE_LIBFFT
	std::cerr << "Spectral space required for solvers, use compile option --libfft=enable" << std::endl;
	exit(-1);
#endif


#if SWEET_REXI_THREAD_PARALLEL_SUM

//#if SWEET_THREADING || 1
	num_local_rexi_par_threads = omp_get_max_threads();
//#else
//	num_local_rexi_par_threads = 1;
//#endif

	if (num_local_rexi_par_threads == 0)
	{
		std::cerr << "FATAL ERROR: omp_get_max_threads == 0" << std::endl;
		exit(-1);
	}
#else
	num_local_rexi_par_threads = 1;
#endif

#if SWEET_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
#else
	mpi_rank = 0;
	num_mpi_ranks = 1;
#endif

	num_global_threads = num_local_rexi_par_threads * num_mpi_ranks;
}



void RexiSWE::cleanup()
{
	for (std::vector<PerThreadVars*>::iterator iter = perThreadVars.begin(); iter != perThreadVars.end(); iter++)
	{
		PerThreadVars* p = *iter;
		delete p;
	}

	perThreadVars.resize(0);
}



RexiSWE::~RexiSWE()
{
	cleanup();

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
	{
		std::cout << "STOPWATCH broadcast: " << stopwatch_broadcast() << std::endl;
		std::cout << "STOPWATCH preprocessing: " << stopwatch_preprocessing() << std::endl;
		std::cout << "STOPWATCH reduce: " << stopwatch_reduce() << std::endl;
		std::cout << "STOPWATCH solve_rexi_terms: " << stopwatch_solve_rexi_terms() << std::endl;
	}
#endif
}



/**
 * setup the REXI
 */
void RexiSWE::setup(
		double i_h,						///< sampling size
		int i_M,						///< number of sampling points
		int i_L,						///< number of sampling points for Gaussian approx.
										///< set to 0 for auto detection

		std::size_t *i_resolution,		///< resolution of domain
		const double *i_domain_size,	///< size of domain

		bool i_rexi_half,				///< use half-pole reduction
		bool i_use_spec_diffs_for_complex_array,	///< use spectral differences for complex arrays in REXI approximation
		int i_helmholtz_solver,			///< Use iterative solver instead of direct solving it in spectral space
		double i_eps					///< Error threshold
)
{
	M = i_M;
	h = i_h;

	helmholtz_solver = i_helmholtz_solver;
	eps = i_eps;
	use_spec_diffs = i_use_spec_diffs_for_complex_array;

	domain_size[0] = i_domain_size[0];
	domain_size[1] = i_domain_size[1];

	rexi.setup(h, M, i_L, i_rexi_half);

	std::size_t N = rexi.alpha.size();
	block_size = N/num_global_threads;
	if (block_size*num_global_threads != N)
		block_size++;

	cleanup();

	perThreadVars.resize(num_local_rexi_par_threads);

	/**
	 * We split the setup from the utilization here.
	 *
	 * This is necessary, since it has to be assured that
	 * the FFTW plans are initialized before using them.
	 */
	if (num_local_rexi_par_threads == 0)
	{
		std::cerr << "FATAL ERROR B: omp_get_max_threads == 0" << std::endl;
		exit(-1);
	}

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
	if (omp_in_parallel())
	{
		std::cerr << "FATAL ERROR X: in parallel region" << std::endl;
		exit(-1);
	}
#endif

	// use a kind of serialization of the input to avoid threading conflicts in the ComplexFFT generation
	for (int j = 0; j < num_local_rexi_par_threads; j++)
	{
#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none) shared(i_resolution,std::cout,j)
#endif
		for (int i = 0; i < num_local_rexi_par_threads; i++)
		{
			if (i != j)
				continue;

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
			if (omp_get_thread_num() != i)
			{
				// leave this dummy std::cout in it to avoid the intel compiler removing this part
				std::cout << "ERROR: thread " << omp_get_thread_num() << " number mismatch " << i << std::endl;
				exit(-1);
			}
#endif

			perThreadVars[i] = new PerThreadVars;

			perThreadVars[i]->op_diff_c_x.setup(i_resolution);
			perThreadVars[i]->op_diff_c_y.setup(i_resolution);
			perThreadVars[i]->op_diff2_c_x.setup(i_resolution);
			perThreadVars[i]->op_diff2_c_y.setup(i_resolution);
			perThreadVars[i]->eta.setup(i_resolution);
			perThreadVars[i]->eta0.setup(i_resolution);
			perThreadVars[i]->u0.setup(i_resolution);
			perThreadVars[i]->v0.setup(i_resolution);
			perThreadVars[i]->h_sum.setup(i_resolution);
			perThreadVars[i]->u_sum.setup(i_resolution);
			perThreadVars[i]->v_sum.setup(i_resolution);
		}
	}

	if (num_local_rexi_par_threads == 0)
	{
		std::cerr << "FATAL ERROR C: omp_get_max_threads == 0" << std::endl;
		exit(-1);
	}

	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
		if (perThreadVars[i]->op_diff_c_x.data == nullptr)
		{
			std::cerr << "ARRAY NOT INITIALIZED!!!!" << std::endl;
			exit(-1);
		}
	}

#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none)  shared(i_domain_size,i_use_spec_diffs_for_complex_array, std::cout)
#endif
	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
//		int global_thread_id = omp_get_thread_num() + mpi_rank*num_local_rexi_par_threads;
		if (omp_get_thread_num() != i)
		{
			// leave this dummy std::cout in it to avoid the intel compiler removing this part
			std::cout << "ERROR: thread " << omp_get_thread_num() << " number mismatch " << i << std::endl;
			exit(-1);
		}
#else

//#if SWEET_REXI_THREAD_PARALLEL_SUM
//		int global_thread_id = 0;
//#endif

#endif

		if (perThreadVars[i]->op_diff_c_x.data == nullptr)
		{
			std::cout << "ERROR, data == nullptr" << std::endl;
			exit(-1);
		}

		perThreadVars[i]->op_diff_c_x.op_setup_diff_x(i_domain_size, i_use_spec_diffs_for_complex_array);
		perThreadVars[i]->op_diff_c_y.op_setup_diff_y(i_domain_size, i_use_spec_diffs_for_complex_array);
		perThreadVars[i]->op_diff2_c_x.op_setup_diff2_x(i_domain_size, i_use_spec_diffs_for_complex_array);
		perThreadVars[i]->op_diff2_c_y.op_setup_diff2_y(i_domain_size, i_use_spec_diffs_for_complex_array);

		// initialize all values to account for first touch policy reason
		perThreadVars[i]->eta.setAll(0, 0);
		perThreadVars[i]->eta0.setAll(0, 0);
		perThreadVars[i]->u0.setAll(0, 0);
		perThreadVars[i]->v0.setAll(0, 0);

		perThreadVars[i]->h_sum.setAll(0, 0);
		perThreadVars[i]->u_sum.setAll(0, 0);
		perThreadVars[i]->v_sum.setAll(0, 0);
	}


#if SWEET_BENCHMARK_REXI
	stopwatch_preprocessing.reset();
	stopwatch_broadcast.reset();
	stopwatch_reduce.reset();
	stopwatch_solve_rexi_terms.reset();
#endif
}


/**
 * Solve REXI with implicit time stepping
 *
 * U_t = L U(0)
 *
 * (U(tau) - U(0)) / tau = L U(tau)
 *
 * <=> U(tau) - U(0) = L U(tau) tau
 *
 * <=> U(tau) - L tau U(tau) = U(0)
 *
 * <=> (1 - L tau) U(tau) = U(0)
 *
 * <=> (1/tau - L) U(tau) = U(0)/tau
 */
bool RexiSWE::run_timestep_implicit_ts(
	DataArray<2> &io_h,
	DataArray<2> &io_u,
	DataArray<2> &io_v,

	double i_timestep_size,	///< timestep size

	Operators2D &op,
	const SimulationVariables &i_simVars
)
{
	Complex2DArrayFFT eta(io_h.resolution);

	Complex2DArrayFFT eta0(io_h.resolution);
	Complex2DArrayFFT u0(io_u.resolution);
	Complex2DArrayFFT v0(io_v.resolution);

	eta0.loadRealFromDataArray(io_h);
	u0.loadRealFromDataArray(io_u);
	v0.loadRealFromDataArray(io_v);

	double alpha = 1.0/i_timestep_size;

	eta0 = eta0.toSpec() * alpha;
	u0 = u0.toSpec() * alpha;
	v0 = v0.toSpec() * alpha;

	// load kappa (k)
	double kappa = alpha*alpha + i_simVars.sim.f0*i_simVars.sim.f0;

	double eta_bar = i_simVars.setup.h0;
	double g = i_simVars.sim.g;
	std::cout << "hi" <<std::endl;
	//@MArtin, this seems to be buggy
	Complex2DArrayFFT &op_diff_c_x = perThreadVars[0]->op_diff_c_x;
	Complex2DArrayFFT &op_diff_c_y = perThreadVars[0]->op_diff_c_y;
	std::cout << kappa << std::endl;
	std::cout << eta_bar << std::endl;
	Complex2DArrayFFT rhs =
			(kappa/alpha) * eta0
			- eta_bar*(op_diff_c_x(u0) + op_diff_c_y(v0))
			- (i_simVars.sim.f0*eta_bar/alpha) * (op_diff_c_x(v0) - op_diff_c_y(u0))
		;
	std::cout << kappa << std::endl;
	std::cout << "hi1" <<std::endl;
	helmholtz_spectral_solver_spec(kappa, g*eta_bar, rhs, eta, 0);
	std::cout << "hi2" <<std::endl;
	Complex2DArrayFFT uh = u0 - g*op_diff_c_x(eta);
	Complex2DArrayFFT vh = v0 - g*op_diff_c_y(eta);

	Complex2DArrayFFT u1 = alpha/kappa * uh     + i_simVars.sim.f0/kappa * vh;
	Complex2DArrayFFT v1 =    -i_simVars.sim.f0/kappa * uh + alpha/kappa * vh;

	eta.toCart().toDataArrays_Real(io_h);
	u1.toCart().toDataArrays_Real(io_u);
	v1.toCart().toDataArrays_Real(io_v);

	return true;
}


/**
 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
 *
 * See
 * 		doc/rexi/understanding_rexi.pdf
 * for further information
 */
bool RexiSWE::run_timestep(
	DataArray<2> &io_h,
	DataArray<2> &io_u,
	DataArray<2> &io_v,

	double i_timestep_size,	///< timestep size

	Operators2D &op,
	const SimulationVariables &i_parameters,
	bool i_iterative_solver_always_init_zero_solution
)
{
	typedef std::complex<double> complex;

	std::size_t N = rexi.alpha.size();

	io_h.requestDataInCartesianSpace();
	io_u.requestDataInCartesianSpace();
	io_v.requestDataInCartesianSpace();


#if SWEET_MPI

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.start();
#endif

	std::size_t data_size = io_h.resolution[0]*io_h.resolution[1];
	MPI_Bcast(io_h.array_data_cartesian_space, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (std::isnan(io_h.get(0,0)))
		return false;


	MPI_Bcast(io_u.array_data_cartesian_space, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(io_v.array_data_cartesian_space, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.stop();
#endif

#endif




#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none) shared(i_parameters, i_timestep_size, io_h, io_u, io_v, N, std::cout, std::cerr, i_iterative_solver_always_init_zero_solution)
#endif
	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
#if SWEET_BENCHMARK_REXI
		bool stopwatch_measure = false;
	#if SWEET_REXI_THREAD_PARALLEL_SUM
		if (omp_get_thread_num() == 0)
	#endif
			if (mpi_rank == 0)
				stopwatch_measure = true;
#endif

#if SWEET_BENCHMARK_REXI
		if (stopwatch_measure)
			stopwatch_preprocessing.start();
#endif

		double eta_bar = i_parameters.setup.h0;
		double g = i_parameters.sim.g;

		Complex2DArrayFFT &op_diff_c_x = perThreadVars[i]->op_diff_c_x;
		Complex2DArrayFFT &op_diff_c_y = perThreadVars[i]->op_diff_c_y;
//		Complex2DArrayFFT &op_diff2_c_x = perThreadVars[i]->op_diff2_c_x;
//		Complex2DArrayFFT &op_diff2_c_y = perThreadVars[i]->op_diff2_c_y;

		Complex2DArrayFFT &eta0 = perThreadVars[i]->eta0;
		Complex2DArrayFFT &u0 = perThreadVars[i]->u0;
		Complex2DArrayFFT &v0 = perThreadVars[i]->v0;

		Complex2DArrayFFT &h_sum = perThreadVars[i]->h_sum;
		Complex2DArrayFFT &u_sum = perThreadVars[i]->u_sum;
		Complex2DArrayFFT &v_sum = perThreadVars[i]->v_sum;

		Complex2DArrayFFT &eta = perThreadVars[i]->eta;


		/*
		 * INITIALIZATION - THIS IS THE NON-PARALLELIZABLE PART!
		 */
		h_sum.setAll(0, 0);
		u_sum.setAll(0, 0);
		v_sum.setAll(0, 0);

		eta0.loadRealFromDataArray(io_h);
		u0.loadRealFromDataArray(io_u);
		v0.loadRealFromDataArray(io_v);

		if (helmholtz_solver == 0)
		{
			/**
			 * SPECTRAL SOLVER - DO EVERYTHING IN SPECTRAL SPACE
			 */
			// convert to spectral space
			// scale with inverse of tau
			eta0 = eta0.toSpec()*(1.0/i_timestep_size);
			u0 = u0.toSpec()*(1.0/i_timestep_size);
			v0 = v0.toSpec()*(1.0/i_timestep_size);

#if SWEET_REXI_THREAD_PARALLEL_SUM || SWEET_MPI

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
			int local_thread_id = omp_get_thread_num();
#else
			int local_thread_id = 0;
#endif
			int global_thread_id = local_thread_id + num_local_rexi_par_threads*mpi_rank;

			std::size_t start = block_size*global_thread_id;
			std::size_t end = std::min(N, start+block_size);
#else
			std::size_t start = 0;
			std::size_t end = N;
#endif


			// reuse result from previous computations
			// this significantly speeds up the process
			// initial guess
			eta.setAll(0,0);

			/*
			 * DO SUM IN PARALLEL
			 */

			// precompute a bunch of values
			// this would belong to a serial part according to Amdahls law
			//
			// (kappa + lhs_a)\eta = kappa/alpha*\eta_0 - (i_parameters.sim.f0*eta_bar/alpha) * rhs_b + rhs_a
			//
			Complex2DArrayFFT rhs_a = eta_bar*(op_diff_c_x(u0) + op_diff_c_y(v0));
			Complex2DArrayFFT rhs_b = (op_diff_c_x(v0) - op_diff_c_y(u0));

			Complex2DArrayFFT lhs_a = (-g*eta_bar)*(perThreadVars[i]->op_diff2_c_x + perThreadVars[i]->op_diff2_c_y);

#if SWEET_BENCHMARK_REXI
			if (stopwatch_measure)
				stopwatch_preprocessing.stop();
#endif

#if SWEET_BENCHMARK_REXI
			if (stopwatch_measure)
				stopwatch_solve_rexi_terms.start();
#endif

			for (std::size_t n = start; n < end; n++)
			{
				// load alpha (a) and scale by inverse of tau
				// we flip the sign to account for the -L used in exp(\tau (-L))
				complex alpha = rexi.alpha[n]/i_timestep_size;
				complex beta = rexi.beta_re[n];

				// load kappa (k)
				complex kappa = alpha*alpha + i_parameters.sim.f0*i_parameters.sim.f0;

				/*
				 * TODO: we can even get more performance out of this operations
				 * by partly using the real Fourier transformation
				 */
				Complex2DArrayFFT rhs =
						(kappa/alpha) * eta0
						+ (-i_parameters.sim.f0*eta_bar/alpha) * rhs_b
						+ rhs_a
					;

#if 0
				helmholtz_spectral_solver_spec(kappa, g*eta_bar, rhs, eta, i);
#else
				Complex2DArrayFFT lhs = lhs_a.addScalar_Cart(kappa);
				rhs.spec_div_element_wise(lhs, eta);
#endif

				Complex2DArrayFFT uh = u0 + g*op_diff_c_x(eta);
				Complex2DArrayFFT vh = v0 + g*op_diff_c_y(eta);

				Complex2DArrayFFT u1 = (alpha/kappa) * uh     - (i_parameters.sim.f0/kappa) * vh;
				Complex2DArrayFFT v1 = (i_parameters.sim.f0/kappa) * uh + (alpha/kappa) * vh;

				DataArray<2> tmp(h_sum.resolution);

				h_sum += eta.toCart()*beta;
				u_sum += u1.toCart()*beta;
				v_sum += v1.toCart()*beta;
			}

#if SWEET_BENCHMARK_REXI
			if (stopwatch_measure)
				stopwatch_solve_rexi_terms.stop();
#endif

		}
		else
		{
			/*
			 * Use FD solver in CARTESIAN SPACE ONLY
			 */
			if (use_spec_diffs)
			{
				std::cerr << "Using FD solvers only makes sense if spectral diffs is NOT activated for REXI!" << std::endl;
				exit(-1);
			}

			// convert to spectral space
			// scale with inverse of tau
			eta0 = eta0.toSpec()*(1.0/i_timestep_size);
			u0 = u0.toSpec()*(1.0/i_timestep_size);
			v0 = v0.toSpec()*(1.0/i_timestep_size);

#if SWEET_REXI_THREAD_PARALLEL_SUM || SWEET_MPI

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
			int local_thread_id = omp_get_thread_num();
#else
			int local_thread_id = 0;
#endif
			int global_thread_id = local_thread_id + num_local_rexi_par_threads*mpi_rank;

			std::size_t start = block_size*global_thread_id;
			std::size_t end = std::min(N, start+block_size);

#else
			std::size_t start = 0;
			std::size_t end = N;
#endif

			// only at the beginning
			if (!i_iterative_solver_always_init_zero_solution)
				eta.setAll(0,0);

			/*
			 * DO SUM IN PARALLEL
			 */
			Complex2DArrayFFT rhs_a = eta_bar*(op_diff_c_x(u0) + op_diff_c_y(v0));
			Complex2DArrayFFT rhs_b = (op_diff_c_x(v0) - op_diff_c_y(u0));

#if SWEET_BENCHMARK_REXI
		if (stopwatch_measure)
			stopwatch_preprocessing.stop();
#endif

#if SWEET_BENCHMARK_REXI
			if (stopwatch_measure)
				stopwatch_solve_rexi_terms.start();
#endif

			for (std::size_t n = start; n < end; n++)
			{
				// load alpha (a) and scale by inverse of tau
				// we flip the sign to account for the -L used in exp(\tau (-L))
				complex alpha = -rexi.alpha[n]/i_timestep_size;
				complex beta = -rexi.beta_re[n];

				// load kappa (k)
				complex kappa = alpha*alpha + i_parameters.sim.f0*i_parameters.sim.f0;

				/*
				 * TODO: we can even get more performance out of this operations
				 * by partly using the real Fourier transformation
				 */
				Complex2DArrayFFT rhs =
						(kappa/alpha) * eta0
						- rhs_a
						- (i_parameters.sim.f0*eta_bar/alpha) * rhs_b
					;

				rhs = rhs.toCart();

				Complex2DArrayFFT eta(rhs.resolution);

				// don't reuse old solution?
				if (i_iterative_solver_always_init_zero_solution)
					eta.setAll(0,0);

				int max_iters = 2000;

				bool retval = true;
				switch (helmholtz_solver)
				{
				default:
					{
						{
						std::cout << "Helmholtz solver IDs:" << std::endl;
						std::cout << "	0: Spectral solver" << std::endl;
						std::cout << "	1: Jacobi" << std::endl;
						std::cout << "	2: CG" << std::endl;
						std::cout << "	3: MG: Jacobi" << std::endl;
						std::cout << "	4: MG: CG" << std::endl;
						}
					}
					exit(-1);
					break;
/*

				case 0:
					helmholtz_spectral_solver_spec(kappa, g*eta_bar, rhs, eta, i);
					break;
*/

				case 1:
					max_iters *= 10;
//					std::cout << "kappa: " << kappa << std::endl;
					retval = RexiSWE_HelmholtzSolver::smoother_jacobi(
							kappa,
							g*eta_bar,
							rhs,
							eta,
							domain_size,
							eps,
							max_iters,
							0.5,
							0
						);

					break;

				case 2:
					retval = RexiSWE_HelmholtzSolver::smoother_conjugate_gradient(
							kappa,
							g*eta_bar,
							rhs,
							eta,
							domain_size,
							eps,
							max_iters,
							1.0,
							0
						);

					break;

				case 3:
					retval = RexiSWE_HelmholtzSolver::multigrid(	// MG with jacobi
						kappa,
						g*eta_bar,
						rhs,
						eta,
						RexiSWE_HelmholtzSolver::smoother_jacobi,
						domain_size,
						eps,
						max_iters,
						1.0,
						2,
						0//simVars.misc.verbosity
					);

					break;

				case 4:
					retval = RexiSWE_HelmholtzSolver::multigrid(	// MG with jacobi
						kappa,
						g*eta_bar,
						rhs,
						eta,
						RexiSWE_HelmholtzSolver::smoother_conjugate_gradient,
						domain_size,
						eps,
						max_iters,
						1.0,
						2,
						0//simVars.misc.verbosity
					);

					break;
				}

				if (!retval)
				{
					double residual = RexiSWE_HelmholtzSolver::helmholtz_iterative_get_residual_rms(
							kappa,
							g*eta_bar,
							rhs,
							eta,
							domain_size
						);

					std::cerr << "NO CONVERGENCE AFTER " << max_iters << " iterations for K=" << kappa << ", gh=" << g*eta_bar << ", eps=" << eps << ", residual=" << residual << std::endl;
					std::cerr << "alpha: " << alpha << std::endl;
					std::cerr << "kappa: " << kappa << std::endl;
					exit(-1);
				}

				eta = eta.toSpec();

				Complex2DArrayFFT uh = u0 - g*op_diff_c_x(eta);
				Complex2DArrayFFT vh = v0 - g*op_diff_c_y(eta);

				Complex2DArrayFFT u1 = (alpha/kappa) * uh	+ (i_parameters.sim.f0/kappa) * vh;
				Complex2DArrayFFT v1 =    (-i_parameters.sim.f0/kappa) * uh	+ (alpha/kappa) * vh;

				h_sum += eta.toCart()*beta;
				u_sum += u1.toCart()*beta;
				v_sum += v1.toCart()*beta;
			}
#if SWEET_BENCHMARK_REXI
			if (stopwatch_measure)
				stopwatch_solve_rexi_terms.stop();
#endif
		}
	}

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_reduce.start();
#endif


#if SWEET_REXI_THREAD_PARALLEL_SUM
	io_h.set_all(0);
	io_u.set_all(0);
	io_v.set_all(0);

	for (int n = 0; n < num_local_rexi_par_threads; n++)
	{
		// sum real-valued elements
		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < io_h.array_data_cartesian_length; i++)
			io_h.array_data_cartesian_space[i] += perThreadVars[n]->h_sum.data[i<<1];

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < io_h.array_data_cartesian_length; i++)
			io_u.array_data_cartesian_space[i] += perThreadVars[n]->u_sum.data[i<<1];

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < io_h.array_data_cartesian_length; i++)
			io_v.array_data_cartesian_space[i] += perThreadVars[n]->v_sum.data[i<<1];
	}

#else

	io_h = perThreadVars[0]->h_sum.getRealWithDataArray();
	io_u = perThreadVars[0]->u_sum.getRealWithDataArray();
	io_v = perThreadVars[0]->v_sum.getRealWithDataArray();

#endif


#if SWEET_MPI
	DataArray<2> tmp(io_h.resolution);

	int retval = MPI_Reduce(io_h.array_data_cartesian_space, tmp.array_data_cartesian_space, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (retval != MPI_SUCCESS)
	{
		std::cerr << "MPI FAILED!" << std::endl;
		exit(1);
	}

	std::swap(io_h.array_data_cartesian_space, tmp.array_data_cartesian_space);

	MPI_Reduce(io_u.array_data_cartesian_space, tmp.array_data_cartesian_space, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(io_u.array_data_cartesian_space, tmp.array_data_cartesian_space);

	MPI_Reduce(io_v.array_data_cartesian_space, tmp.array_data_cartesian_space, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(io_v.array_data_cartesian_space, tmp.array_data_cartesian_space);
#endif


#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_reduce.stop();
#endif

	return true;
}



inline std::complex<double> conj(const std::complex<double> &v)
{
	return std::complex<double>(v.real(), -v.imag());
}



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
void RexiSWE::run_timestep_direct_solution(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		double i_timestep_size,	///< timestep size

		Operators2D &op,
		const SimulationVariables &i_simVars
)
{
	typedef std::complex<double> complex;

	double eta_bar = i_simVars.setup.h0;
	double g = i_simVars.sim.g;
	double f = i_simVars.sim.f0;
	complex I(0.0,1.0);

	Complex2DArrayFFT i_h(io_h.resolution);
	Complex2DArrayFFT i_u(io_h.resolution);
	Complex2DArrayFFT i_v(io_h.resolution);

	Complex2DArrayFFT o_h(io_h.resolution);
	Complex2DArrayFFT o_u(io_h.resolution);
	Complex2DArrayFFT o_v(io_h.resolution);

	i_h.loadRealFromDataArray(io_h);
	i_h = i_h.toSpec();

	i_u.loadRealFromDataArray(io_u);
	i_u = i_u.toSpec();

	i_v.loadRealFromDataArray(io_v);
	i_v = i_v.toSpec();

	double s0 = i_simVars.sim.domain_size[0];
	double s1 = i_simVars.sim.domain_size[1];

	for (std::size_t ik1 = 0; ik1 < i_h.resolution[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < i_h.resolution[0]; ik0++)
		{
			if (ik0 == i_h.resolution[0]/2 || ik1 == i_h.resolution[1]/2)
			{
				o_h.set(ik1, ik0, 0, 0);
				o_u.set(ik1, ik0, 0, 0);
				o_v.set(ik1, ik0, 0, 0);
			}

			complex U_hat[3];
			U_hat[0] = i_h.get(ik1, ik0);
			U_hat[1] = i_u.get(ik1, ik0);
			U_hat[2] = i_v.get(ik1, ik0);

			double k0, k1;
			if (ik0 < i_h.resolution[0]/2)
				k0 = (double)ik0;
			else
				k0 = (double)((int)ik0-(int)i_h.resolution[0]);

			if (ik1 < i_h.resolution[1]/2)
				k1 = (double)ik1;
			else
				k1 = (double)((int)ik1-(int)i_h.resolution[1]);

			/*
			 * dimensionful formulation
			 * see doc/swe_solution_for_L
			 */

			double H0 = eta_bar;

			//////////////////////////////////////
			// GENERATED CODE START
			//////////////////////////////////////
			complex eigenvalues[3];
			complex eigenvectors[3][3];

			if (k0 == 0 && k1 == 0)
			{
//					complex wg = std::sqrt((complex)f*f*s0*s0*s1*s1);

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*f;
				eigenvalues[2] = f;

				eigenvectors[0][0] = 1.00000000000000;
				eigenvectors[0][1] = 0.0;
				eigenvectors[0][2] = 0.0;
				eigenvectors[1][0] = 0.0;
				eigenvectors[1][1] = -1.0*I;
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 0.0;
				eigenvectors[2][1] = I;
				eigenvectors[2][2] = 1.00000000000000;
			}
			else if (k0 == 0)
			{
//					complex wg = std::sqrt((complex)s0*s0*(f*f*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1));

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvalues[2] = -1.0*I*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k1*k1 - 1.0*f*f*s1*s1);

				eigenvectors[0][0] = (1.0/2.0)*I*1.0/M_PI*f*1.0/g*1.0/k1*s1;
				eigenvectors[0][1] = 1.00000000000000;
				eigenvectors[0][2] = 0.0;
				eigenvectors[1][0] = -2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][1] = -1.0*I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][1] = I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][2] = 1.00000000000000;
			}
			else if (k1 == 0)
			{
//					complex wg = std::sqrt((complex)s1*s1*(f*f*s0*s0 + 4.0*M_PI*M_PI*g*g*k0*k0));

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
				eigenvalues[2] = -1.0*I*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);

				eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
				eigenvectors[0][1] = 0.0;
				eigenvectors[0][2] = 1.00000000000000;
				eigenvectors[1][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
				eigenvectors[1][1] = -1.0*I*1.0/f*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
				eigenvectors[2][1] = 1.0/f*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);
				eigenvectors[2][2] = 1.00000000000000;
			}
			else
			{
//					complex K2 = M_PI*M_PI*k0*k0 + M_PI*M_PI*k1*k1;
				complex w = std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);

//					complex wg = std::sqrt((complex)f*f*s0*s0*s1*s1 + 4.0*M_PI*M_PI*g*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1*s0*s0);
				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);
				eigenvalues[2] = -1.0*I*1.0/s0*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 - 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 - 1.0*f*f*s0*s0*s1*s1);

				eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
				eigenvectors[0][1] = -1.0*1.0/k0*k1*s0*1.0/s1;
				eigenvectors[0][2] = 1.00000000000000;
				eigenvectors[1][0] = 2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 + f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 - 1.0*I*f*w);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = -2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 - 1.0*f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 + I*f*w);
				eigenvectors[2][2] = 1.00000000000000;
			}




			//////////////////////////////////////
			// GENERATED CODE END
			//////////////////////////////////////


			if (f == 0)
			{
				/*
				 * override if f == 0, see ./sympy_L_spec_decomposition.py executed with LNr=4
				 */
				if (k0 != 0 || k1 != 0)
				{
					double K2 = K2;

					eigenvalues[0] = 0.0;
					eigenvalues[1] = -2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);
					eigenvalues[2] = 2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);

					eigenvectors[0][0] = 0.0;
					eigenvectors[0][1] = -1.0*k1/sqrt(k0*k0 + k1*k1);
					eigenvectors[0][2] = k0/sqrt(k0*k0 + k1*k1);
					eigenvectors[1][0] = -1.0*sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[1][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[1][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][0] = sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
				}
				else
				{

					eigenvalues[0] = 0.0;
					eigenvalues[1] = 0.0;
					eigenvalues[2] = 0.0;

					eigenvectors[0][0] = 1.00000000000000;
					eigenvectors[0][1] = 0.0;
					eigenvectors[0][2] = 0.0;
					eigenvectors[1][0] = 0.0;
					eigenvectors[1][1] = 1.00000000000000;
					eigenvectors[1][2] = 0.0;
					eigenvectors[2][0] = 0.0;
					eigenvectors[2][1] = 0.0;
					eigenvectors[2][2] = 1.00000000000000;
				}
			}


			/*
			 * Compute inverse of Eigenvectors.
			 * This generalizes to the case that the Eigenvectors are not orthonormal.
			 */
			complex eigenvectors_inv[3][3];

			eigenvectors_inv[0][0] =  (eigenvectors[1][1]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][1]);
			eigenvectors_inv[0][1] = -(eigenvectors[0][1]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][1]);
			eigenvectors_inv[0][2] =  (eigenvectors[0][1]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][1]);

			eigenvectors_inv[1][0] = -(eigenvectors[1][0]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][0]);
			eigenvectors_inv[1][1] =  (eigenvectors[0][0]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][0]);
			eigenvectors_inv[1][2] = -(eigenvectors[0][0]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][0]);

			eigenvectors_inv[2][0] =  (eigenvectors[1][0]*eigenvectors[2][1] - eigenvectors[1][1]*eigenvectors[2][0]);
			eigenvectors_inv[2][1] = -(eigenvectors[0][0]*eigenvectors[2][1] - eigenvectors[0][1]*eigenvectors[2][0]);
			eigenvectors_inv[2][2] =  (eigenvectors[0][0]*eigenvectors[1][1] - eigenvectors[0][1]*eigenvectors[1][0]);

			complex s = eigenvectors[0][0]*eigenvectors_inv[0][0] + eigenvectors[0][1]*eigenvectors_inv[1][0] + eigenvectors[0][2]*eigenvectors_inv[2][0];

			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++)
					eigenvectors_inv[j][i] /= s;


			// check
			for (int j = 0; j < 3; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					if (
							std::isnan(eigenvectors[j][i].real()) || std::isinf(eigenvectors[j][i].real())	||
							std::isnan(eigenvectors[j][i].imag()) || std::isinf(eigenvectors[j][i].imag())
					)
					{
						std::cerr << "Invalid number in Eigenvector " << j << " detected: " << eigenvectors[j][0] << ", " << eigenvectors[j][1] << ", " << eigenvectors[j][2] << std::endl;
					}

					if (
							std::isnan(eigenvectors_inv[j][i].real()) || std::isinf(eigenvectors_inv[j][i].real())	||
							std::isnan(eigenvectors_inv[j][i].imag()) || std::isinf(eigenvectors_inv[j][i].imag())
					)
					{
						std::cerr << "Invalid number in inverse of Eigenvector " << j << " detected: " << eigenvectors_inv[j][0] << ", " << eigenvectors_inv[j][1] << ", " << eigenvectors_inv[j][2] << std::endl;
					}
				}
			}

			/*
			 * Solve based on previously computed data.
			 * Note, that this data can be also precomputed and reused every time.
			 */
			complex UEV0_sp[3];
			for (int k = 0; k < 3; k++)
			{
				UEV0_sp[k] = {0, 0};
				for (int j = 0; j < 3; j++)
					UEV0_sp[k] += eigenvectors_inv[j][k] * U_hat[j];
			}

			complex omega[3];
			omega[0] = std::exp(-I*eigenvalues[0]*i_timestep_size);
			omega[1] = std::exp(-I*eigenvalues[1]*i_timestep_size);
			omega[2] = std::exp(-I*eigenvalues[2]*i_timestep_size);

			complex U_hat_sp[3];
			for (int k = 0; k < 3; k++)
			{
				U_hat_sp[k] = {0, 0};
				for (int j = 0; j < 3; j++)
					U_hat_sp[k] += eigenvectors[j][k] * omega[j] * UEV0_sp[j];
			}

			o_h.set(ik1, ik0, U_hat_sp[0]);
			o_u.set(ik1, ik0, U_hat_sp[1]);
			o_v.set(ik1, ik0, U_hat_sp[2]);
		}
	}

	io_h = o_h.toCart().getRealWithDataArray();
	io_u = o_u.toCart().getRealWithDataArray();
	io_v = o_v.toCart().getRealWithDataArray();
}


