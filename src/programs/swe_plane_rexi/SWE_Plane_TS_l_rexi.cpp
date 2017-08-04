/*
 * SWE_Plane_TS_l_rexi.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Plane_TS_l_rexi.hpp"
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>

#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>

#define FOOBAR123 0

#if FOOBAR123
std::vector<std::complex<double>> x_rexi_alpha;
std::vector<std::complex<double>> x_rexi_beta_re;
#endif

void SWE_Plane_TS_l_rexi::setup(
		REXI_SimulationVariables &i_rexi,
		const std::string &i_function_name
)
{
	rexiSimVars = &i_rexi;

	domain_size[0] = simVars.sim.domain_size[0];
	domain_size[1] = simVars.sim.domain_size[1];

	if (rexiSimVars->use_direct_solution)
	{
		ts_l_direct.setup(i_function_name);
		return;
	}

	if (rexiSimVars->use_next_generation)
	{
		bool retval = rexiNG.auto_load(
				i_function_name,
				rexiSimVars->ng_N,	/// N
				rexiNG.None(),	/// max_error
				rexiSimVars->ng_max_error_double_precision,			/// max_error_double_precision
				rexiSimVars->ng_test_min,
				rexiSimVars->ng_test_max,
				rexiNG.None(),	/// basis_function_scaling
				rexiSimVars->ng_h, //rexiNG.None(),	/// basis_function_spacing
				rexiNG.None(),	/// basis_function rat shift

				rexiSimVars->use_half_poles,

				rexiSimVars->ng_faf_dir
			);

		if (!retval)
			FatalError(std::string("Not able to find coefficients for given constraints for function "+i_function_name));

		if (simVars.misc.verbosity > 0)
			std::cout << "Loaded REXI coefficients from file '" << rexiNG.fafcoeffs.filename << "'" << std::endl;

		if (simVars.misc.verbosity > 3)
		{
			rexiNG.fafcoeffs.output();
			rexiNG.fafcoeffs.outputWeights();
			rexiNG.output();
		}

		rexi_alpha = rexiNG.alpha;
		rexi_beta_re = rexiNG.beta_re;

#if 1
		for (std::size_t i = 0; i < rexi_alpha.size(); i++)
		{
			rexi_alpha[i] = -std::conj(rexi_alpha[i]);
			rexi_beta_re[i] = -std::conj(rexi_beta_re[i]);
		}
#endif
	}
	else
	{
		rexi.setup(0, rexiSimVars->h, rexiSimVars->M, rexiSimVars->L, rexiSimVars->use_half_poles, rexiSimVars->normalization);

		rexi_alpha = rexi.alpha;
		rexi_beta_re = rexi.beta_re;

		if (simVars.misc.verbosity > 3)
			rexi.output();

#if FOOBAR123
		// fake setup
		//		rexi.setup(0, rexiSimVars->h, rexiSimVars->M, rexiSimVars->L, !rexiSimVars->use_half_poles, rexiSimVars->normalization);
		rexi.setup(0, rexiSimVars->h, rexiSimVars->M, rexiSimVars->L, rexiSimVars->use_half_poles, rexiSimVars->normalization);
		x_rexi_alpha = rexi.alpha;
		x_rexi_beta_re = rexi.beta_re;
#endif

	}

	std::cout << "Halving rule = " << rexiSimVars->use_half_poles << std::endl;
	std::cout << "Number of total REXI coefficients N = " << rexi_alpha.size() << std::endl;

	std::size_t N = rexi_alpha.size();
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

	PlaneDataConfig *planeDataConfig_local = this->planeDataConfig;

	// use a kind of serialization of the input to avoid threading conflicts in the ComplexFFT generation
	for (int j = 0; j < num_local_rexi_par_threads; j++)
	{
#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none) shared(planeDataConfig_local, std::cout,j)
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

			perThreadVars[i]->op.setup(planeDataConfig_local, domain_size);

			perThreadVars[i]->eta.setup(planeDataConfig_local);
			perThreadVars[i]->eta0.setup(planeDataConfig_local);
			perThreadVars[i]->u0.setup(planeDataConfig_local);
			perThreadVars[i]->v0.setup(planeDataConfig_local);
			perThreadVars[i]->h_sum.setup(planeDataConfig_local);
			perThreadVars[i]->u_sum.setup(planeDataConfig_local);
			perThreadVars[i]->v_sum.setup(planeDataConfig_local);
		}
	}

	if (num_local_rexi_par_threads == 0)
	{
		std::cerr << "FATAL ERROR C: omp_get_max_threads == 0" << std::endl;
		exit(-1);
	}

	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
		if (perThreadVars[i]->op.diff_c_x.physical_space_data == nullptr)
		{
			std::cerr << "ARRAY NOT INITIALIZED!!!!" << std::endl;
			exit(-1);
		}
	}

#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none)  shared(planeDataConfig_local, std::cout)
#endif
	for (int i = 0; i < num_local_rexi_par_threads; i++)
	{
#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
		if (omp_get_thread_num() != i)
		{
			// leave this dummy std::cout in it to avoid the intel compiler removing this part
			std::cout << "ERROR: thread " << omp_get_thread_num() << " number mismatch " << i << std::endl;
			exit(-1);
		}
#else

#endif

		if (perThreadVars[i]->eta.physical_space_data == nullptr)
		{
			std::cout << "ERROR, data == nullptr" << std::endl;
			exit(-1);
		}

		perThreadVars[i]->op.setup(planeDataConfig_local, domain_size);

		// initialize all values to account for first touch policy reason
		perThreadVars[i]->eta.spectral_set_all(0, 0);
		perThreadVars[i]->eta0.spectral_set_all(0, 0);

		perThreadVars[i]->u0.spectral_set_all(0, 0);
		perThreadVars[i]->v0.spectral_set_all(0, 0);

		perThreadVars[i]->h_sum.spectral_set_all(0, 0);
		perThreadVars[i]->u_sum.spectral_set_all(0, 0);
		perThreadVars[i]->v_sum.spectral_set_all(0, 0);

	}


#if SWEET_BENCHMARK_REXI
	stopwatch_preprocessing.reset();
	stopwatch_broadcast.reset();
	stopwatch_reduce.reset();
	stopwatch_solve_rexi_terms.reset();
#endif
}



#if 0
void spectral_zero_nyquist_phase(PlaneData &i_planeData)
{
	/*
	 *
	 */
	for (int j = 0; j < i_planeData.planeDataConfig->spectral_data_size[1]; j++)
	{
		i_planeData.p_spectral_get_nonconstref(j, i_planeData.planeDataConfig->spectral_data_size[0]-1).real(0);
		i_planeData.p_spectral_get_nonconstref(j, i_planeData.planeDataConfig->spectral_data_size[0]-1).imag(0);
	}

	for (int i = 0; i < i_planeData.planeDataConfig->spectral_data_size[0]-1; i++)
	{
		i_planeData.p_spectral_get_nonconstref(i_planeData.planeDataConfig->spectral_data_size[1]/2, i).real(0);
		i_planeData.p_spectral_get_nonconstref(i_planeData.planeDataConfig->spectral_data_size[1]/2, i).imag(0);
	}
}



void spectral_half_nyquist_phase(PlaneData &i_planeData)
{
	/*
	 *
	 */
	for (int j = 0; j < i_planeData.planeDataConfig->spectral_data_size[1]; j++)
		i_planeData.p_spectral_get_nonconstref(j, i_planeData.planeDataConfig->spectral_data_size[0]-1) *= 0.5;

	for (int i = 0; i < i_planeData.planeDataConfig->spectral_data_size[0]-1; i++)
		i_planeData.p_spectral_get_nonconstref(i_planeData.planeDataConfig->spectral_data_size[1]/2, i) *= 0.5;
}
#endif



void SWE_Plane_TS_l_rexi::run_timestep(
		const PlaneData &i_h_pert,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_pert,	///< prognostic variables
		PlaneData &o_u,	///< prognostic variables
		PlaneData &o_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
#if !FOOBAR123

	/// WARNING: i_h_pert might be identical to o_h_pert
	run_timestep_real(i_h_pert, i_u, i_v, o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);

#elif 1

#if 0
	op.diff_c_x.print_physicalData_zeroNumZero();
	std::cout << std::endl;
	op.diff_c_y.print_physicalData_zeroNumZero();
	exit(1);
#endif

	PlaneData tmp_h_pert = i_h_pert;
	PlaneData tmp_u = i_u;
	PlaneData tmp_v = i_v;

	run_timestep_real(i_h_pert, i_u, i_v, o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);


	PlaneData x_o_h_pert(planeDataConfig);
	PlaneData x_o_u(planeDataConfig);
	PlaneData x_o_v(planeDataConfig);

	std::swap(rexi_alpha, x_rexi_alpha);
	std::swap(rexi_beta_re, x_rexi_beta_re);

	run_timestep_real(tmp_h_pert, tmp_u, tmp_v, x_o_h_pert, x_o_u, x_o_v, i_dt, i_simulation_timestamp);

	std::swap(rexi_alpha, x_rexi_alpha);
	std::swap(rexi_beta_re, x_rexi_beta_re);

	PlaneData diff_o_h_pert = x_o_h_pert - o_h_pert;
	PlaneData diff_o_u = x_o_u - o_u;
	PlaneData diff_o_v = x_o_v - o_v;

#if 1
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

//	std::cout << rexi_alpha.size() << ", " << x_rexi_alpha.size() << std::endl;
	std::cout << "tmp_h_pert" << std::endl;
	tmp_h_pert.print_spectralData_zeroNumZero(1e-12);
	std::cout << "prog_h" << std::endl;
	diff_o_h_pert.print_spectralData_zeroNumZero(1e-12);
	std::cout << "u" << std::endl;
	diff_o_u.print_spectralData_zeroNumZero(1e-12);
	std::cout << "v" << std::endl;
	diff_o_v.print_spectralData_zeroNumZero(1e-12);
//	exit(1);
#endif

#else

	PlaneData tmp_h_pert = i_h_pert;
	PlaneData tmp_u = i_u;
	PlaneData tmp_v = i_v;

	run_timestep_real(i_h_pert, i_u, i_v, o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);

	static int i = 0;

	static PlaneData x_i_h_pert(planeDataConfig);
	static PlaneData x_i_u(planeDataConfig);
	static PlaneData x_i_v(planeDataConfig);

	static PlaneData x_o_h_pert(planeDataConfig);
	static PlaneData x_o_u(planeDataConfig);
	static PlaneData x_o_v(planeDataConfig);

	if (i == 0)
	{
		x_i_h_pert = tmp_h_pert;
		x_i_u = tmp_u;
		x_i_v = tmp_v;

//		i = 1;
	}

	std::swap(rexi_alpha, x_rexi_alpha);
	std::swap(rexi_beta_re, x_rexi_beta_re);

	run_timestep_real(x_i_h_pert, x_i_u, x_i_v, x_o_h_pert, x_o_u, x_o_v, i_dt, i_simulation_timestamp);

	std::swap(rexi_alpha, x_rexi_alpha);
	std::swap(rexi_beta_re, x_rexi_beta_re);

	x_i_h_pert = o_h_pert;
	x_i_u = o_u;
	x_i_v = o_v;

	PlaneData diff_o_h_pert = x_o_h_pert - o_h_pert;
	PlaneData diff_o_u = x_o_u - o_u;
	PlaneData diff_o_v = x_o_v - o_v;

#if 0
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	tmp_h_pert.print_spectralData_zeroNumZero();
	tmp_h_pert.request_data_physical();
	std::cout << std::endl;
	tmp_h_pert.print_spectralData_zeroNumZero();
#endif

#if 0
//	std::cout << rexi_alpha.size() << ", " << x_rexi_alpha.size() << std::endl;
	std::cout << "i_prog_h" << std::endl;
	i_h_pert.print_spectralData_zeroNumZero(1e-12);
	std::cout << "prog_h" << std::endl;
	diff_o_h_pert.print_spectralData_zeroNumZero(1e-12);
	std::cout << "u" << std::endl;
	diff_o_u.print_spectralData_zeroNumZero(1e-12);
	std::cout << "v" << std::endl;
	diff_o_v.print_spectralData_zeroNumZero(1e-12);
//	exit(1);
#endif
#endif
}



void SWE_Plane_TS_l_rexi::run_timestep_real(
		const PlaneData &i_h_pert,	///< prognostic variables
		const PlaneData &i_u,		///< prognostic variables
		const PlaneData &i_v,		///< prognostic variables

		PlaneData &o_h_pert,	///< prognostic variables
		PlaneData &o_u,			///< prognostic variables
		PlaneData &o_v,			///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	final_timestep = false;

	if (simVars.rexi.use_direct_solution)
	{
		o_h_pert = i_h_pert;
		o_u = i_u;
		o_v = i_v;
		ts_l_direct.run_timestep(o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);
		return;
	}

	if (i_dt <= 0)
		FatalError("Only constant time step size allowed");


	typedef std::complex<double> complex;

	std::size_t max_N = rexi_alpha.size();


	/*
	 * Request physical or spectral here to avoid parallel race conditions
	 */
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	i_h_pert.request_data_physical();
	i_u.request_data_physical();
	i_v.request_data_physical();
#else
	i_h_pert.request_data_spectral();
	i_u.request_data_spectral();
	i_v.request_data_spectral();
#endif

#if SWEET_MPI

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.start();
#endif

	std::size_t data_size = i_h_pert.planeDataConfig->physical_array_data_number_of_elements;
	MPI_Bcast(i_h_pert.physical_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (std::isnan(i_h_pert.p_physical_get(0,0)))
	{
		final_timestep = true;
		return;
	}


	MPI_Bcast(i_u.physical_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_v.physical_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.stop();
#endif

#endif



#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none) shared(i_dt, i_h_pert, i_u, i_v, max_N, std::cout, std::cerr)
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

		double eta_bar = simVars.sim.h0;
		double g = simVars.sim.gravitation;

		PlaneOperatorsComplex &opc = perThreadVars[i]->op;

		PlaneDataComplex &eta0 = perThreadVars[i]->eta0;
		PlaneDataComplex &u0 = perThreadVars[i]->u0;
		PlaneDataComplex &v0 = perThreadVars[i]->v0;

		PlaneDataComplex &h_sum = perThreadVars[i]->h_sum;
		PlaneDataComplex &u_sum = perThreadVars[i]->u_sum;
		PlaneDataComplex &v_sum = perThreadVars[i]->v_sum;

		PlaneDataComplex &eta = perThreadVars[i]->eta;


		/*
		 * INITIALIZATION - THIS IS THE NON-PARALLELIZABLE PART!
		 */
		h_sum.spectral_set_zero();
		u_sum.spectral_set_zero();
		v_sum.spectral_set_zero();


#if !SWEET_USE_PLANE_SPECTRAL_SPACE
		eta0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(i_h_pert);
		u0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(i_u);
		v0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(i_v);
#else
		if (simVars.rexi.use_half_poles)
		{
			eta0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(i_h_pert);
			u0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(i_u);
			v0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(i_v);
		}
		else
		{
			eta0 = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(i_h_pert);
			u0 = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(i_u);
			v0 = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(i_v);
		}
#endif

		/**
		 * SPECTRAL SOLVER - DO EVERYTHING IN SPECTRAL SPACE
		 *
		 * (alpha+dt L)U = U0
		 *
		 * (alpha/dt+L) (dt U) = U0
		 *
		 * (alpha/dt+L) U = U0/dt
		 */
		// convert to spectral space
		// scale with inverse of tau
		double inv_dt = (1.0/i_dt);
		eta0 = eta0*inv_dt;
		u0 = u0*inv_dt;
		v0 = v0*inv_dt;

#if SWEET_REXI_THREAD_PARALLEL_SUM || SWEET_MPI

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
		int local_thread_id = omp_get_thread_num();
#else
		int local_thread_id = 0;
#endif
		int global_thread_id = local_thread_id + num_local_rexi_par_threads*mpi_rank;

		std::size_t start = std::min(max_N, block_size*global_thread_id);
		std::size_t end = std::min(max_N, start+block_size);
#else
		std::size_t start = 0;
		std::size_t end = max_N;
#endif

		/*
		 * DO SUM IN PARALLEL
		 */
		//
		// precompute a bunch of values
		// this would belong to a serial part according to Amdahl's law
		//
		// (kappa + lhs_a)\eta = kappa/alpha*\eta_0 - (i_parameters.sim.f0*eta_bar/alpha) * rhs_b + rhs_a
		//
		PlaneDataComplex rhs_a = eta_bar*(opc.diff_c_x(u0) + opc.diff_c_y(v0));
		PlaneDataComplex rhs_b = (opc.diff_c_x(v0) - opc.diff_c_y(u0));

		PlaneDataComplex lhs_a = (-g*eta_bar)*(perThreadVars[i]->op.diff2_c_x + perThreadVars[i]->op.diff2_c_y);

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
			complex alpha = rexi_alpha[n]/i_dt;
			complex beta = rexi_beta_re[n];

			if (simVars.sim.f0 == 0)
			{
				/*
				 * TODO: we can even get more performance out of this operations
				 * by partly using the real Fourier transformation
				 */
				PlaneDataComplex rhs =
						eta0*alpha
						+ eta_bar*(opc.diff_c_x(u0) + opc.diff_c_y(v0))
					;

				PlaneDataComplex lhs_a = (-g*eta_bar)*(perThreadVars[i]->op.diff2_c_x + perThreadVars[i]->op.diff2_c_y);
				PlaneDataComplex lhs = lhs_a.spectral_addScalarAll(alpha*alpha);

				eta = rhs.spectral_div_element_wise(lhs);

				PlaneDataComplex u1 = (u0 + g*opc.diff_c_x(eta))*(1.0/alpha);
				PlaneDataComplex v1 = (v0 + g*opc.diff_c_y(eta))*(1.0/alpha);

				h_sum += eta*beta;
				u_sum += u1*beta;
				v_sum += v1*beta;
			}
			else
			{
				// load kappa (k)
				complex kappa = alpha*alpha + simVars.sim.f0*simVars.sim.f0;

				/*
				 * TODO: we can even get more performance out of this operations
				 * by partly using the real Fourier transformation
				 */
				PlaneDataComplex rhs =
						(kappa/alpha) * eta0
						+ (-simVars.sim.f0*eta_bar/alpha) * rhs_b
						+ rhs_a
					;

				PlaneDataComplex lhs = lhs_a.spectral_addScalarAll(kappa);
				eta = rhs.spectral_div_element_wise(lhs);

				PlaneDataComplex uh = u0 + g*opc.diff_c_x(eta);
				PlaneDataComplex vh = v0 + g*opc.diff_c_y(eta);

				PlaneDataComplex u1 = (alpha/kappa) * uh     - (simVars.sim.f0/kappa) * vh;
				PlaneDataComplex v1 = (simVars.sim.f0/kappa) * uh + (alpha/kappa) * vh;

				PlaneData tmp(h_sum.planeDataConfig);

				h_sum += eta*beta;
				u_sum += u1*beta;
				v_sum += v1*beta;
			}
		}

#if SWEET_BENCHMARK_REXI
		if (stopwatch_measure)
			stopwatch_solve_rexi_terms.stop();
#endif
	}

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_reduce.start();
#endif

#if SWEET_REXI_THREAD_PARALLEL_SUM


#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	o_h_pert.physical_set_zero();
	o_u.physical_set_zero();
	o_v.physical_set_zero();

	for (int n = 0; n < num_local_rexi_par_threads; n++)
	{
		perThreadVars[n]->h_sum.request_data_physical();
		perThreadVars[n]->u_sum.request_data_physical();
		perThreadVars[n]->v_sum.request_data_physical();

		// sum real-valued elements
		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < i_h_pert.planeDataConfig->physical_array_data_number_of_elements; i++)
			o_h_pert.physical_space_data[i] += perThreadVars[n]->h_sum.physical_space_data[i].real();

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < i_h_pert.planeDataConfig->physical_array_data_number_of_elements; i++)
			o_u.physical_space_data[i] += perThreadVars[n]->u_sum.physical_space_data[i].real();

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < i_h_pert.planeDataConfig->physical_array_data_number_of_elements; i++)
			o_v.physical_space_data[i] += perThreadVars[n]->v_sum.physical_space_data[i].real();
	}
#else

	o_h_pert.spectral_set_zero();
	o_u.spectral_set_zero();
	o_v.spectral_set_zero();

	for (int n = 0; n < num_local_rexi_par_threads; n++)
	{
		perThreadVars[n]->h_sum.request_data_spectral();
		perThreadVars[n]->u_sum.request_data_spectral();
		perThreadVars[n]->v_sum.request_data_spectral();

		PlaneData tmp(planeDataConfig);

		o_h_pert = o_h_pert + Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(perThreadVars[n]->h_sum);
		o_u = o_u + Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(perThreadVars[n]->u_sum);
		o_v = o_v + Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(perThreadVars[n]->v_sum);
	}
#endif


#else

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	o_h_pert = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->h_sum);
	o_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->u_sum);
	o_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->v_sum);

#else

	if (simVars.rexi.use_half_poles)
	{
		o_h_pert = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->h_sum);
		o_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->u_sum);
		o_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->v_sum);
	}
	else
	{
		o_h_pert = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(perThreadVars[0]->h_sum);
		o_u = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(perThreadVars[0]->u_sum);
		o_v = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(perThreadVars[0]->v_sum);
	}
#endif


#endif


#if SWEET_MPI

#if !SWEET_USE_PLANE_SPECTRAL_SPACE ||1
	PlaneData tmp(o_h_pert.planeDataConfig);

	o_h_pert.request_data_physical();
	int retval = MPI_Reduce(o_h_pert.physical_space_data, tmp.physical_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (retval != MPI_SUCCESS)
	{
		std::cerr << "MPI FAILED!" << std::endl;
		exit(1);
	}

	std::swap(o_h_pert.physical_space_data, tmp.physical_space_data);

	o_u.request_data_physical();
	MPI_Reduce(o_u.physical_space_data, tmp.physical_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(o_u.physical_space_data, tmp.physical_space_data);

	o_v.request_data_physical();
	MPI_Reduce(o_v.physical_space_data, tmp.physical_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(o_v.physical_space_data, tmp.physical_space_data);
#else
#error "TODO: spectral version"
#endif

#endif


#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_reduce.stop();
#endif
}



void SWE_Plane_TS_l_rexi::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	run_timestep(
			io_h, io_u, io_v,
			io_h, io_u, io_v,
			i_dt,
			i_simulation_timestamp
		);
}


void SWE_Plane_TS_l_rexi::cleanup()
{

	for (std::vector<PerThreadVars*>::iterator iter = perThreadVars.begin(); iter != perThreadVars.end(); iter++)
	{
		PerThreadVars* p = *iter;
		delete p;
	}

	perThreadVars.resize(0);
}




SWE_Plane_TS_l_rexi::SWE_Plane_TS_l_rexi(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		planeDataConfig(op.planeDataConfig),
		ts_l_direct(i_simVars, i_op)
{
#if !SWEET_USE_LIBFFT
	std::cerr << "Spectral space required for solvers, use compile option --libfft=enable" << std::endl;
	exit(-1);
#endif


#if SWEET_REXI_THREAD_PARALLEL_SUM

	num_local_rexi_par_threads = omp_get_max_threads();

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



void SWE_Plane_TS_l_rexi::MPI_quitWorkers(
		PlaneDataConfig *i_planeDataConfig
)
{
#if SWEET_MPI
	PlaneData dummyData(i_planeDataConfig);
	dummyData.physical_set_all(NAN);

	MPI_Bcast(dummyData.physical_space_data, dummyData.planeDataConfig->physical_array_data_number_of_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}




SWE_Plane_TS_l_rexi::~SWE_Plane_TS_l_rexi()
{
	cleanup();
}

