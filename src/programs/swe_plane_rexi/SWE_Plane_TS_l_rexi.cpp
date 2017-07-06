/*
 * SWE_Plane_TS_l_rexi.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Plane_TS_l_rexi.hpp"
#include <sweet/plane/PlaneDataComplex.hpp>

#include <sweet/plane/PlaneOperatorsComplex.hpp>

#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>


void SWE_Plane_TS_l_rexi::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	final_timestep = false;

	if (simVars.rexi.use_direct_solution)
	{
		ts_l_direct.run_timestep(io_h, io_u, io_v, o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time);
		return;
	}

	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time - i_simulation_timestamp;


	typedef std::complex<double> complex;

	std::size_t max_N = rexi.alpha.size();

	io_h.request_data_physical();
	io_u.request_data_physical();
	io_v.request_data_physical();

#if SWEET_MPI

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.start();
#endif

	std::size_t data_size = io_h.planeDataConfig->physical_array_data_number_of_elements;
	MPI_Bcast(io_h.physical_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (std::isnan(io_h.physical_get(0,0)))
	{
		final_timestep = true;
		return;
	}


	MPI_Bcast(io_u.physical_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(io_v.physical_space_data, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.stop();
#endif

#endif



#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none) shared(i_fixed_dt, io_h, io_u, io_v, max_N, std::cout, std::cerr)
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
		h_sum.spectral_set_all(0, 0);
		u_sum.spectral_set_all(0, 0);
		v_sum.spectral_set_all(0, 0);

		eta0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_h);
		u0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_u);
		v0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_v);


		/**
		 * SPECTRAL SOLVER - DO EVERYTHING IN SPECTRAL SPACE
		 */
		// convert to spectral space
		// scale with inverse of tau
		eta0 = eta0*(1.0/i_fixed_dt);
		u0 = u0*(1.0/i_fixed_dt);
		v0 = v0*(1.0/i_fixed_dt);

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


		// reuse result from previous computations
		// this significantly speeds up the process
		// initial guess
//		eta.spectral_set_all(0,0);

		/*
		 * DO SUM IN PARALLEL
		 */

		// precompute a bunch of values
		// this would belong to a serial part according to Amdahls law
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
			// we flip the sign to account for the -L used in exp(\tau (-L))
			complex alpha = DQStuff::convertComplex<double>(rexi.alpha[n])/i_fixed_dt;
			complex beta = DQStuff::convertComplex<double>(rexi.beta_re[n]);

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
//			rhs.spectral_div_element_wise(lhs, eta);
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
	io_h.physical_set_all(0);
	io_u.physical_set_all(0);
	io_v.physical_set_all(0);

	for (int n = 0; n < num_local_rexi_par_threads; n++)
	{
		perThreadVars[n]->h_sum.request_data_physical();
		perThreadVars[n]->u_sum.request_data_physical();
		perThreadVars[n]->v_sum.request_data_physical();

		// sum real-valued elements
		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < io_h.planeDataConfig->physical_array_data_number_of_elements; i++)
			io_h.physical_space_data[i] += perThreadVars[n]->h_sum.physical_space_data[i].real();

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < io_h.planeDataConfig->physical_array_data_number_of_elements; i++)
			io_u.physical_space_data[i] += perThreadVars[n]->u_sum.physical_space_data[i].real();

		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < io_h.planeDataConfig->physical_array_data_number_of_elements; i++)
			io_v.physical_space_data[i] += perThreadVars[n]->v_sum.physical_space_data[i].real();
	}

#else

	io_h = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->h_sum);
	io_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->u_sum);
	io_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(perThreadVars[0]->v_sum);

#endif


#if SWEET_MPI
	PlaneData tmp(io_h.planeDataConfig);

	io_h.request_data_physical();
	int retval = MPI_Reduce(io_h.physical_space_data, tmp.physical_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (retval != MPI_SUCCESS)
	{
		std::cerr << "MPI FAILED!" << std::endl;
		exit(1);
	}

	std::swap(io_h.physical_space_data, tmp.physical_space_data);

	io_u.request_data_physical();
	MPI_Reduce(io_u.physical_space_data, tmp.physical_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(io_u.physical_space_data, tmp.physical_space_data);

	io_v.request_data_physical();
	MPI_Reduce(io_v.physical_space_data, tmp.physical_space_data, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(io_v.physical_space_data, tmp.physical_space_data);
#endif


#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_reduce.stop();
#endif

	o_dt = i_fixed_dt;
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


void SWE_Plane_TS_l_rexi::setup(
		double i_h,						///< sampling size
		int i_M,						///< number of sampling points
		int i_L,						///< number of sampling points for Gaussian approximation
										///< set to 0 for auto detection

		bool i_rexi_half,				///< use half-pole reduction
		bool i_rexi_normalization		///< REXI normalization
)
{
	M = i_M;
	h = i_h;

	domain_size[0] = simVars.sim.domain_size[0];
	domain_size[1] = simVars.sim.domain_size[1];

	rexi.setup(0, h, M, i_L, i_rexi_half, i_rexi_normalization);

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

