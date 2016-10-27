/*
 * SWE_Sphere_REXI.cpp
 *
 *  Created on: 25 Oct 2016
 *      Author: martin
 */

#include <cassert>
#include <rexi/REXI.hpp>


#include "SWE_Sphere_REXI.hpp"
#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>


#ifndef SWEET_REXI_THREAD_PARALLEL_SUM
#	define SWEET_REXI_THREAD_PARALLEL_SUM 1
#endif

/*
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



SWE_Sphere_REXI::SWE_Sphere_REXI()	:
	sphereDataConfig(nullptr),
	sphereDataConfigRexi(nullptr)
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



void SWE_Sphere_REXI::cleanup()
{
	for (std::vector<PerThreadVars*>::iterator iter = perThreadVars.begin(); iter != perThreadVars.end(); iter++)
	{
		PerThreadVars* p = *iter;
		delete p;
	}

	perThreadVars.resize(0);

	sphereDataConfigRexi = nullptr;
	sphereDataConfig = nullptr;
}



SWE_Sphere_REXI::~SWE_Sphere_REXI()
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


void SWE_Sphere_REXI::get_workload_start_end(
		std::size_t &o_start,
		std::size_t &o_end
)
{
	std::size_t max_N = rexi.alpha.size();

#if SWEET_REXI_THREAD_PARALLEL_SUM || SWEET_MPI

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
	int local_thread_id = omp_get_thread_num();
#else
	int local_thread_id = 0;
#endif

	int global_thread_id = local_thread_id + num_local_rexi_par_threads*mpi_rank;

	assert(block_size >= 0);
	assert(global_thread_id >= 0);

	o_start = std::min(max_N, block_size*global_thread_id);
	o_end = std::min(max_N, o_start+block_size);

#else

	o_start = 0;
	o_end = max_N;

#endif
}



/**
 * setup the REXI
 */
void SWE_Sphere_REXI::setup(
		double i_h,						///< sampling size
		int i_M,						///< number of sampling points
		int i_L,						///< number of sampling points for Gaussian approx.
										///< set to 0 for auto detection

		SphereDataConfig *i_sphereDataConfig,
		SimulationVariables::Coefficients *i_simCoeffs,
		double i_timestep_size,

		bool i_rexi_half,				///< use half-pole reduction
		bool i_use_robert_functions,	///< use Robert functions
		int i_rexi_use_extended_modes,
		bool i_use_coriolis_rexi_formulation
)
{
	cleanup();

	M = i_M;
	h = i_h;

	sphereDataConfig = i_sphereDataConfig;
	simCoeffs = i_simCoeffs;
	timestep_size = i_timestep_size;

	use_robert_functions = i_use_robert_functions;
	rexi_use_extended_modes = i_rexi_use_extended_modes;
	use_coriolis_rexi_formulation = i_use_coriolis_rexi_formulation;
	use_rexi_preallocation = false;


	if (rexi_use_extended_modes == 0)
	{
		sphereDataConfigRexi = sphereDataConfig;
	}
	else
	{
		// Add modes only along latitude since these are the "problematic" modes
		sphereConfigRexiAddedModes.setupAdditionalModes(
				sphereDataConfig,
				rexi_use_extended_modes,	// TODO: Extend SPH wrapper to also support m != n to set this guy to 0
				rexi_use_extended_modes
		);
		sphereDataConfigRexi = &sphereConfigRexiAddedModes;
	}


	rexi.setup(h, M, i_L, i_rexi_half);

	std::size_t N = rexi.alpha.size();
	block_size = N/num_global_threads;
	if (block_size*num_global_threads != N)
		block_size++;

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
#	pragma omp parallel for schedule(static,1) default(none) shared(std::cout,j)
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

			std::size_t start, end;
			get_workload_start_end(start, end);
			int local_size = (int)end-(int)start;

//#pragma omp critical
//			std::cout << start << "\t" << end << "\t" << local_size << "\t" << block_size << std::endl;

#if SWEET_DEBUG
			if (local_size < 0)
				FatalError("local_size < 0");
#endif

			perThreadVars[i]->alpha.resize(local_size);
			perThreadVars[i]->beta_re.resize(local_size);

			perThreadVars[i]->accum_phi.setup(sphereDataConfigRexi);
			perThreadVars[i]->accum_u.setup(sphereDataConfigRexi);
			perThreadVars[i]->accum_v.setup(sphereDataConfigRexi);


			for (std::size_t n = start; n < end; n++)
			{
				int thread_local_idx = n-start;

				perThreadVars[i]->alpha[thread_local_idx] = rexi.alpha[n];
				perThreadVars[i]->beta_re[thread_local_idx] = rexi.beta_re[n];
			}

			if (use_rexi_preallocation)
			{
				if (use_robert_functions)
					perThreadVars[i]->rexiSPHRobert_vector.resize(local_size);
				else
					perThreadVars[i]->rexiSPH_vector.resize(local_size);

				for (std::size_t n = start; n < end; n++)
				{
					int thread_local_idx = n-start;

					if (use_robert_functions)
					{
						perThreadVars[i]->rexiSPHRobert_vector[thread_local_idx].setup(
								sphereDataConfigRexi,
								perThreadVars[i]->alpha[thread_local_idx],
								perThreadVars[i]->beta_re[thread_local_idx],
								simCoeffs->earth_radius,
								simCoeffs->coriolis_omega,
								simCoeffs->h0 * simCoeffs->gravitation,
								timestep_size,
								use_coriolis_rexi_formulation
						);
					}
					else
					{
						perThreadVars[i]->rexiSPH_vector[thread_local_idx].setup(
								sphereDataConfigRexi,
								perThreadVars[i]->alpha[thread_local_idx],
								perThreadVars[i]->beta_re[thread_local_idx],
								simCoeffs->earth_radius,
								simCoeffs->coriolis_omega,
								simCoeffs->h0*simCoeffs->gravitation,
								timestep_size,
								use_coriolis_rexi_formulation
						);
					}
				}
			}
		}
	}

	if (num_local_rexi_par_threads == 0)
	{
		std::cerr << "FATAL ERROR C: omp_get_max_threads == 0" << std::endl;
		exit(-1);
	}



#if SWEET_BENCHMARK_REXI
	stopwatch_preprocessing.reset();
	stopwatch_broadcast.reset();
	stopwatch_reduce.reset();
	stopwatch_solve_rexi_terms.reset();
#endif
}



/**
 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
 *
 * See
 * 		doc/rexi/understanding_rexi.pdf
 *
 * for further information
 */
bool SWE_Sphere_REXI::run_timestep_rexi(
	SphereData &io_prog_h0,
	SphereData &io_prog_u0,
	SphereData &io_prog_v0,

	double i_timestep_size,	///< timestep size

	const SimulationVariables &i_parameters
)
{
	io_prog_h0.request_data_spectral();
	io_prog_u0.request_data_spectral();
	io_prog_v0.request_data_spectral();


#if SWEET_MPI

	/*
	 * TODO: Maybe we should measure this for the 2nd rank!!!
	 * The reason could be since Bcast might already return before the packages were actually received!
	 */
	#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.start();
	#endif

	std::size_t data_size = io_prog_h0.sphereDataConfigRexi->spectral_array_data_number_of_elements*2;

	MPI_Bcast(io_prog_h0.array_data_spectral_space, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (std::isnan(io_prog_h0.spectral_get(0,0)))
		return false;

	MPI_Bcast(io_prog_u0.array_data_spectral_space, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(io_prog_v0.array_data_spectral_space, data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	#if SWEET_BENCHMARK_REXI
	if (mpi_rank == 0)
		stopwatch_broadcast.stop();
	#endif

#endif


#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for schedule(static,1) default(none) shared(i_parameters, i_timestep_size, io_prog_h0, io_prog_u0, io_prog_v0, std::cout, std::cerr)
#endif
	for (int thread_id = 0; thread_id < num_local_rexi_par_threads; thread_id++)
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

		std::size_t start, end;
		get_workload_start_end(start, end);

		/*
		 * DO SUM IN PARALLEL
		 */
		SphereData thread_prog_phi0(sphereDataConfigRexi);
		SphereData thread_prog_u0(sphereDataConfigRexi);
		SphereData thread_prog_v0(sphereDataConfigRexi);

#if SWEET_DEBUG
		/**
		 * THIS ASSERTION IS VERY IMPORTANT!
		 * OTHERWISE io_prog_*0 will be converted to
		 * spectral space *in parallel* with write
		 * access raceconditions
		 */
		if (	!io_prog_h0.spectral_space_data_valid	||
				!io_prog_u0.spectral_space_data_valid	||
				!io_prog_v0.spectral_space_data_valid
			)
		{
			FatalError("SPECTRAL DATA NOT AVAILABLE, BUT REQUIRED!");
		}
#endif

		if (rexi_use_extended_modes == 0)
		{
			thread_prog_phi0 = io_prog_h0*simCoeffs->gravitation;
			thread_prog_u0 = io_prog_u0;
			thread_prog_v0 = io_prog_v0;
		}
		else
		{
			(io_prog_h0*simCoeffs->gravitation).spectral_copyToDifferentModes(thread_prog_phi0);
			io_prog_u0.spectral_copyToDifferentModes(thread_prog_u0);
			io_prog_v0.spectral_copyToDifferentModes(thread_prog_v0);
		}


		SphereData tmp_prog_phi(sphereDataConfigRexi);
		SphereData tmp_prog_u(sphereDataConfigRexi);
		SphereData tmp_prog_v(sphereDataConfigRexi);

		perThreadVars[thread_id]->accum_phi.spectral_set_zero();
		perThreadVars[thread_id]->accum_u.spectral_set_zero();
		perThreadVars[thread_id]->accum_v.spectral_set_zero();


		for (std::size_t workload_idx = start; workload_idx < end; workload_idx++)
		{
			int local_idx = workload_idx-start;

			std::complex<double> &alpha = perThreadVars[thread_id]->alpha[local_idx];
			std::complex<double> &beta_re = perThreadVars[thread_id]->beta_re[local_idx];


			if (use_robert_functions)
			{
				if (use_rexi_preallocation)
				{
					perThreadVars[thread_id]->rexiSPHRobert_vector[local_idx].solve(
							thread_prog_phi0, thread_prog_u0, thread_prog_v0,
							tmp_prog_phi, tmp_prog_u, tmp_prog_v
						);
				}
				else
				{
					SWERexi_SPHRobert rexiSPHRobert;

					rexiSPHRobert.setup(
							sphereDataConfigRexi,
							alpha,
							beta_re,
							simCoeffs->earth_radius,
							simCoeffs->coriolis_omega,
							simCoeffs->h0*simCoeffs->gravitation,
							i_timestep_size,
							use_coriolis_rexi_formulation
					);

					rexiSPHRobert.solve(
							thread_prog_phi0, thread_prog_u0, thread_prog_v0,
							tmp_prog_phi, tmp_prog_u, tmp_prog_v
						);
				}
			}
			else
			{
				if (use_rexi_preallocation)
				{
					perThreadVars[thread_id]->rexiSPH_vector[local_idx].solve(
							thread_prog_phi0, thread_prog_u0, thread_prog_v0,
							tmp_prog_phi, tmp_prog_u, tmp_prog_v
						);
				}
				else
				{
					SWERexi_SPH rexiSPH;
					rexiSPH.setup(
							sphereDataConfigRexi,
							alpha,
							beta_re,
							simCoeffs->earth_radius,
							simCoeffs->coriolis_omega,
							simCoeffs->h0*simCoeffs->gravitation,
							i_timestep_size,
							use_coriolis_rexi_formulation
					);

					rexiSPH.solve(
							thread_prog_phi0, thread_prog_u0, thread_prog_v0,
							tmp_prog_phi, tmp_prog_u, tmp_prog_v
						);
				}
			}


			perThreadVars[thread_id]->accum_phi += tmp_prog_phi;
			perThreadVars[thread_id]->accum_u += tmp_prog_u;
			perThreadVars[thread_id]->accum_v += tmp_prog_v;
		}

#if SWEET_DEBUG
		if (	!io_prog_h0.spectral_space_data_valid	||
				!io_prog_u0.spectral_space_data_valid	||
				!io_prog_v0.spectral_space_data_valid
			)
		{
			FatalError("SPECTRAL DATA NOT AVAILABLE, BUT REQUIRED!");
		}
#endif

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

	io_prog_h0.physical_set_zero();
	io_prog_u0.physical_set_zero();
	io_prog_v0.physical_set_zero();

	for (int thread_id = 0; thread_id < num_local_rexi_par_threads; thread_id++)
	{
		if (rexi_use_extended_modes == 0)
		{
			assert(io_prog_h0.sphereDataConfig->spectral_array_data_number_of_elements == perThreadVars[0]->accum_phi.sphereDataConfig->spectral_array_data_number_of_elements);


			perThreadVars[thread_id]->accum_phi.request_data_physical();

			#pragma omp parallel for schedule(static) default(none) shared(io_prog_h0, thread_id)
			for (int i = 0; i < io_prog_h0.sphereDataConfig->physical_array_data_number_of_elements; i++)
				io_prog_h0.physical_space_data[i] += perThreadVars[thread_id]->accum_phi.physical_space_data[i];


			perThreadVars[thread_id]->accum_u.request_data_physical();

			#pragma omp parallel for schedule(static) default(none) shared(io_prog_u0, thread_id)
			for (int i = 0; i < io_prog_u0.sphereDataConfig->physical_array_data_number_of_elements; i++)
				io_prog_u0.physical_space_data[i] += perThreadVars[thread_id]->accum_u.physical_space_data[i];


			perThreadVars[thread_id]->accum_v.request_data_physical();

			#pragma omp parallel for schedule(static) default(none) shared(io_prog_v0, thread_id)
			for (int i = 0; i < io_prog_v0.sphereDataConfig->physical_array_data_number_of_elements; i++)
				io_prog_v0.physical_space_data[i] += perThreadVars[thread_id]->accum_v.physical_space_data[i];
		}
		else
		{
			assert(io_prog_h0.sphereDataConfig->spectral_array_data_number_of_elements == sphereDataConfig->spectral_array_data_number_of_elements);

			SphereData tmp(sphereDataConfig);


			perThreadVars[thread_id]->accum_phi.spectral_copyToDifferentModes(tmp);
			tmp.request_data_physical();
			#pragma omp parallel for schedule(static) default(none) shared(io_prog_h0, tmp)
			for (int i = 0; i < io_prog_h0.sphereDataConfig->physical_array_data_number_of_elements; i++)
				io_prog_h0.physical_space_data[i] += tmp.physical_space_data[i]*(1.0/simCoeffs->gravitation);


			perThreadVars[thread_id]->accum_u.spectral_copyToDifferentModes(tmp);
			tmp.request_data_physical();
			#pragma omp parallel for schedule(static) default(none) shared(io_prog_u0, tmp)
			for (int i = 0; i < io_prog_u0.sphereDataConfig->physical_array_data_number_of_elements; i++)
				io_prog_u0.physical_space_data[i] += tmp.physical_space_data[i];


			perThreadVars[thread_id]->accum_v.spectral_copyToDifferentModes(tmp);
			tmp.request_data_physical();
			#pragma omp parallel for schedule(static) default(none) shared(io_prog_v0, tmp)
			for (int i = 0; i < io_prog_v0.sphereDataConfig->physical_array_data_number_of_elements; i++)
				io_prog_v0.physical_space_data[i] += tmp.physical_space_data[i];
		}
	}

#else

	if (rexi_use_extended_modes == 0)
	{
		io_prog_h0 = perThreadVars[0]->accum_phi*(1.0/simCoeffs->gravitation);
		io_prog_u0 = perThreadVars[0]->accum_u;
		io_prog_v0 = perThreadVars[0]->accum_v;
	}
	else
	{
		(perThreadVars[0]->accum_phi*(1.0/simCoeffs->gravitation)).spectral_copyToDifferentModes(io_prog_h0);
		perThreadVars[0]->accum_u.spectral_copyToDifferentModes(io_prog_u0);
		perThreadVars[0]->accum_v.spectral_copyToDifferentModes(io_prog_v0);
	}

#endif

	io_prog_h0.request_data_physical();
	io_prog_u0.request_data_physical();
	io_prog_v0.request_data_physical();

#if SWEET_MPI
	SphereData tmp(io_prog_h0.resolution);

	int retval = MPI_Reduce(io_prog_h0.array_data_physical_space, tmp.array_data_physical_space, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (retval != MPI_SUCCESS)
	{
		std::cerr << "MPI FAILED!" << std::endl;
		exit(1);
	}

	std::swap(io_prog_h0.array_data_physical_space, tmp.array_data_physical_space);

	MPI_Reduce(io_prog_u0.array_data_physical_space, tmp.array_data_physical_space, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(io_prog_u0.array_data_physical_space, tmp.array_data_physical_space);

	MPI_Reduce(io_prog_v0.array_data_physical_space, tmp.array_data_physical_space, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	std::swap(io_prog_v0.array_data_physical_space, tmp.array_data_physical_space);
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


