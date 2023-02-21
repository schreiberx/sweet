/*
 * PlaneDataConfig.hpp
 *
 *  Created on: 17 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PLANE_DATA_CONFIG_HPP_
#define PLANE_DATA_CONFIG_HPP_

/*
 * DOCUMENTATION:
 *
 * The aliasing strategies are described in the folder
 * 		doc/antialiasing/
 */

#include <fftw3.h>
#include <iostream>
#include <stdlib.h>
#include <iostream>
#include <complex>
#include <iomanip>
#include <cmath>
#include <sweet/SWEETError.hpp>
#include <sweet/TransformationPlans.hpp>



#if SWEET_THREADING_SPACE
#	include <omp.h>
#endif



/*
 * This precompiler directive should !!! NEVER !!! be used
 */
#ifndef SWEET_USE_PLANE_SPECTRAL_SPACE
	#define SWEET_USE_PLANE_SPECTRAL_SPACE	1
#endif

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		#error "Dealiasing only available with plane-spectral-space enabled"
	#endif
#endif

/*
 * Activating this option via the precompiler also activates the FFTW library
 */
#ifndef SWEET_USE_LIBFFT
	#define SWEET_USE_LIBFFT	1
#endif

/*
 * Activating the dealiasing creates plans and additional information for dealiasing strategies
 */
#ifndef SWEET_USE_PLANE_SPECTRAL_DEALIASING
	#define SWEET_USE_PLANE_SPECTRAL_DEALIASING 1
#endif


class PlaneDataConfig
{
	friend class PlaneOperators;
	friend class PlaneOperatorsComplex;
	friend class PlaneData;
	friend class PlaneDataComplex;


public:

	/// Physical resolution: Number of cells in each dimension
	std::size_t physical_res[2];

	/// size of physical storage, is identical to physical resolution
	std::size_t physical_data_size[2];

	/// number of real valued data
	std::size_t physical_array_data_number_of_elements;

	/// Different strategies to cope with transformation plans
	TransformationPlans::TRANSFORMATION_PLAN_CACHE reuse_spectral_transformation_plans;


public:

#if SWEET_USE_LIBFFT

	/*
	 * Number of spectral modes
	 *
	 * This is not related to the storage size!
	 *
	 * Also, this includes the modes which cannot be represented according
	 * to the Shannon-Nyquist theorem!
	 */
	std::size_t spectral_modes[2];

	std::size_t spectral_real_modes[2];

	/// allocated size for spectral data for each modes
	/// This storage size is for the real-to-complex transformation
	std::size_t spectral_data_size[2];

	/// iteration ranges for updating data in spectrum
	/// 1st index (left): the id of the range,
	/// 2nd index (middle): dimension the range,
	/// 3rd index (last one): start and end (exclusive) index
	std::size_t spectral_data_iteration_ranges[2][2][2];

	/// total number of complex-valued data elements in spectral space
	std::size_t spectral_array_data_number_of_elements;

private:
	/*
	 * FFTW related stuff
	 */
	fftw_plan	fftw_plan_forward;
	fftw_plan	fftw_plan_backward;

	/// FFTW scaling related stuff for backward transformation
	/// WARNING: FFTW doesn't implement a symmetric FFTW
	/// We only to the rescaling for the backward transformation
	double fftw_backward_scale_factor;


public:
	/// allocated size for spectral data in case of complex data in physical space
	std::size_t spectral_complex_data_size[2];

	/// total number of elements in spectrum
	std::size_t spectral_complex_array_data_number_of_elements;

	/// iteration range for complex valued space
	/// 1st index (left): the id of the range,
	/// 2nd index (middle): dimension the range,
	/// 3rd index (last one): start and end (exclusive) index
	std::size_t spectral_complex_ranges[4][2][2];

	fftw_plan	fftw_plan_complex_forward;
	fftw_plan	fftw_plan_complex_backward;

#endif

	bool fftw_initialized;


	std::string getUniqueIDString()	const
	{
		return getConfigInformationString();
	}

	std::string getConfigInformationString()	const
	{
		std::ostringstream buf;
		buf <<

#if SWEET_USE_LIBFFT
				"M" << spectral_modes[0] << "," << spectral_modes[1] << "_" <<
#endif
				"N" << physical_res[0] << "," << physical_res[1];

		return buf.str();
	}


public:
	void printInformation()	const
	{
		std::cout << std::endl;
		std::cout << "physical_res: " << physical_res[0] << ", " << physical_res[1] << std::endl;
		std::cout << "physical_data_size: " << physical_data_size[0] << ", " << physical_data_size[1] << std::endl;
		std::cout << "physical_array_data_number_of_elements: " << physical_array_data_number_of_elements << std::endl;

#if SWEET_USE_LIBFFT
		std::cout << std::endl;
		std::cout << "spectral_modes: " << spectral_modes[0] << ", " << spectral_modes[1] << std::endl;
		std::cout << "spectral_real_modes: " << spectral_real_modes[0] << ", " << spectral_real_modes[1] << std::endl;
		std::cout << "spectral_data_size: " << spectral_data_size[0] << ", " << spectral_data_size[1] << std::endl;
		std::cout << "spectral_array_data_number_of_elements: " << spectral_array_data_number_of_elements << std::endl;
		std::cout << std::endl;
		std::cout << "spectral_data_iteration_ranges [0][0]: " << spectral_data_iteration_ranges[0][0][0] << ", " << spectral_data_iteration_ranges[0][0][1] << std::endl;
		std::cout << "spectral_data_iteration_ranges [0][1]: " << spectral_data_iteration_ranges[0][1][0] << ", " << spectral_data_iteration_ranges[0][1][1] << std::endl;
		std::cout << "spectral_data_iteration_ranges [1][0]: " << spectral_data_iteration_ranges[1][0][0] << ", " << spectral_data_iteration_ranges[1][0][1] << std::endl;
		std::cout << "spectral_data_iteration_ranges [1][1]: " << spectral_data_iteration_ranges[1][1][0] << ", " << spectral_data_iteration_ranges[1][1][1] << std::endl;
		std::cout << std::endl;
		std::cout << "spectral_complex_data_size: " << spectral_complex_data_size[0] << ", " << spectral_complex_data_size[1] << std::endl;
		std::cout << "spectral_complex_array_data_number_of_elements: " << spectral_complex_array_data_number_of_elements << std::endl;
		std::cout << std::endl;
		std::cout << "spectral_complex_ranges [0][0]: " << spectral_complex_ranges[0][0][0] << ", " << spectral_complex_ranges[0][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [0][1]: " << spectral_complex_ranges[0][1][0] << ", " << spectral_complex_ranges[0][1][1] << std::endl;
		std::cout << "spectral_complex_ranges [1][0]: " << spectral_complex_ranges[1][0][0] << ", " << spectral_complex_ranges[1][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [1][1]: " << spectral_complex_ranges[1][1][0] << ", " << spectral_complex_ranges[1][1][1] << std::endl;
		std::cout << "spectral_complex_ranges [2][0]: " << spectral_complex_ranges[2][0][0] << ", " << spectral_complex_ranges[2][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [2][1]: " << spectral_complex_ranges[2][1][0] << ", " << spectral_complex_ranges[2][1][1] << std::endl;
		std::cout << "spectral_complex_ranges [3][0]: " << spectral_complex_ranges[3][0][0] << ", " << spectral_complex_ranges[3][0][1] << std::endl;
		std::cout << "spectral_complex_ranges [3][1]: " << spectral_complex_ranges[3][1][0] << ", " << spectral_complex_ranges[3][1][1] << std::endl;
		std::cout << std::endl;
#endif
	}


public:
	static
	bool loadWisdom(TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans)
	{
#if 1

		static const char *wisdom_file = "sweet_fftw";

#else
		static const char *wisdom_file = nullptr;

		if (wisdom_file != nullptr)
			return false;

		wisdom_file = getenv("SWEET_FFTW_WISDOM_FROM_FILE");

		if (wisdom_file == nullptr)
		{
			wisdom_file = "sweet_fftw";
			return false;
		}

		std::cout << "Loading SWEET_FFTW_LOAD_WISDOM_FROM_FILE=" << wisdom_file << std::endl;
#endif

		int wisdom_plan_loaded = fftw_import_wisdom_from_filename(wisdom_file);
		if (wisdom_plan_loaded == 0)
		{
			std::cerr << "Failed to load FFTW wisdom from file '" << wisdom_file << "'" << std::endl;
			if (i_reuse_spectral_transformation_plans & TransformationPlans::REQUIRE_LOAD)
			{
				std::cerr << "IMPORTANT: Usage of FFTW wisdoms file enforced, hence exiting here" << std::endl;
				exit(1);
			}
		}

		return true;
	}


public:
	static
	bool storeWisdom()
	{
#if 1

		static const char *wisdom_file = "sweet_fftw";

#else

		static const char *store_wisdom_from_file = nullptr;

		if (wisdom_file != nullptr)
			return false;

		wisdom_file = getenv("SWEET_FFTW_WISDOM_FROM_FILE");

		if (wisdom_file == nullptr)
		{
			wisdom_file = "sweet_fftw";
			return false;
		}

		std::cout << "Loading SWEET_FFTW_LOAD_WISDOM_FROM_FILE=" << wisdom_file << std::endl;
#endif

		int wisdom_plan_loaded = fftw_export_wisdom_to_filename(wisdom_file);
		if (wisdom_plan_loaded == 0)
		{
			std::cerr << "Failed to store FFTW wisdom to file " << wisdom_file << std::endl;
			exit(1);
		}

		return true;
	}


public:
	int& refCounterFftwPlans()
	{
#if SWEET_THREADING_SPACE && SWEET_DEBUG
		if (omp_get_level() != 0)
			SWEETError("PlaneDataConfig is not threadsafe, but called inside parallel region with more than one thread!!!");
#endif

		static int ref_counter = 0;
		return ref_counter;
	}


public:
	PlaneDataConfig()
	{
		physical_res[0] = 0;
		physical_res[1] = 0;

#if SWEET_USE_LIBFFT
		spectral_modes[0] = 0;
		spectral_modes[1] = 0;
#endif

		fftw_initialized = false;
	}



private:
	void setup_internal_data(
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		reuse_spectral_transformation_plans = i_reuse_spectral_transformation_plans;

		if (fftw_initialized)
		{
			//
			// refCounter()--;
			// is inside cleanup!
			cleanup_data();
		}
		else
		{
			// use REF counter in case of multiple plans
			// this allows a clean cleanup of fftw library
			fftw_initialized = true;
		}

		// FFTW PLANS are allocated below
		refCounterFftwPlans()++;

		physical_data_size[0] = physical_res[0];
		physical_data_size[1] = physical_res[1];

		physical_array_data_number_of_elements = physical_res[0]*physical_res[1];

#if SWEET_USE_LIBFFT
		if (	physical_res[0] < spectral_modes[0]	||
				physical_res[1] < spectral_modes[1]
		)
		{
			SWEETError("Lower physical resolution than spectral resolution not supported!");
		}

		assert(physical_res[0] > 0);
		assert(physical_res[1] > 0);

		assert(spectral_modes[0] > 0);
		assert(spectral_modes[1] > 0);


#if SWEET_THREADING_TIME_REXI

		// Threaded time parallel sum
		std::cout << "Using REXI parallel sum, hence using only single FFT thread" << std::endl;

#else

	#if SWEET_THREADING_SPACE
		// Is this the first instance?
		if (refCounterFftwPlans() == 1)
		{
//			std::cout << "FFTW THREADING" << std::endl;
			// initialise FFTW with spatial parallelization
			int retval = fftw_init_threads();
			if (retval == 0)
			{
				std::cerr << "ERROR: fftw_init_threads()" << std::endl;
				exit(1);
			}

			int nthreads = 1;
#if 0
			// Don't use omp_get_max_threads(), since this doesn't reflect the actual number of cores in OMP_NUM_THREADS
#pragma omp parallel
#pragma omp master
			nthreads = omp_get_num_threads();
#else
			nthreads = omp_get_max_threads();
#endif

#if 0
			std::cout << "Using " << nthreads << " for FFTW" << std::endl;
			std::cout << "omp_get_max_threads(): " << omp_get_max_threads() << std::endl;
			std::cout << "omp_get_num_procs() " << omp_get_num_procs() << std::endl;
			std::cout << "PAR MASTER omp_get_num_threads() " << nthreads << std::endl;
#endif

			fftw_plan_with_nthreads(nthreads);

		}
	#endif

#endif

		if (refCounterFftwPlans() == 1)
		{
			// load wisdom the first time
			// this must be done after initializing the threading!
			if (i_reuse_spectral_transformation_plans & TransformationPlans::LOAD)
				loadWisdom(reuse_spectral_transformation_plans);
		}


		unsigned int flags = 0;

		// allow destroying input for faster transformations
		//flags |= FFTW_DESTROY_INPUT;

		if (i_reuse_spectral_transformation_plans & TransformationPlans::QUICK)
		{
			flags |= FFTW_ESTIMATE;
		}
		else
		{
			// estimation base don workload
			unsigned int cells = physical_data_size[0]*physical_data_size[1];

			if (cells < 32*32)
				//flags |= FFTW_EXHAUSTIVE;
				flags |= FFTW_MEASURE;
			else if (cells < 256*256)
				flags |= FFTW_MEASURE;
			else
				flags |= FFTW_MEASURE;
				//flags |= FFTW_PATIENT;

			if (i_reuse_spectral_transformation_plans & TransformationPlans::REQUIRE_LOAD)
			{
				flags |= FFTW_WISDOM_ONLY;	// only allow plans from wisdom
			}
		}


		/*
		 * REAL PHYSICAL SPACE DATA (REAL to COMPLEX FFT)
		 */
		{
			// real-to-complex storage representation
			spectral_data_size[0] = physical_data_size[0]/2+1;
			spectral_data_size[1] = physical_data_size[1];

			if ((physical_data_size[0] & 1) == 1)
				SWEETError("Unsupported odd resolution in x-direction");

			if ((physical_data_size[1] & 1) == 1)
				SWEETError("Unsupported odd resolution in y-direction");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING

			/*
			 * For more information, have a look at
			 * doc/software_development_discussions/antialiasing/implementation_strategy.pdf
			 */

			if (	spectral_modes[0] == physical_res[0] ||
					spectral_modes[1] == physical_res[1]
			)
				SWEETError("Aliasing doesn't make sense since physical resolution is identical to spectral");

			spectral_data_iteration_ranges[0][0][0] = 0;
			spectral_data_iteration_ranges[0][0][1] = (physical_data_size[0]-1)/3;
			spectral_data_iteration_ranges[0][1][0] = 0;
			spectral_data_iteration_ranges[0][1][1] = (physical_data_size[1]-1)/3;


#else

			/*
			 * The central mode is the Shannon-Nyquist one
			 * => Remove this one, since this is not tracked in transformations
			 */
			spectral_data_iteration_ranges[0][0][0] = 0;
			spectral_data_iteration_ranges[0][0][1] = spectral_data_size[0]-1;		// Shannon-Nyquist
			spectral_data_iteration_ranges[0][1][0] = 0;
			spectral_data_iteration_ranges[0][1][1] = spectral_data_size[1]/2;

#endif


			spectral_data_iteration_ranges[1][0][0] = spectral_data_iteration_ranges[0][0][0];
			spectral_data_iteration_ranges[1][0][1] = spectral_data_iteration_ranges[0][0][1];
			spectral_data_iteration_ranges[1][1][0] = spectral_data_size[1] - spectral_data_iteration_ranges[0][1][1] + 1;
			spectral_data_iteration_ranges[1][1][1] = spectral_data_size[1];

			spectral_real_modes[0] = spectral_data_iteration_ranges[0][0][1];
			spectral_real_modes[1] = spectral_data_iteration_ranges[0][1][1];

			spectral_array_data_number_of_elements = spectral_data_size[0]*spectral_data_size[1];



			/*
			 * Physical space data
			 */
			double *data_physical = MemBlockAlloc::alloc<double>(physical_array_data_number_of_elements*sizeof(double));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < physical_array_data_number_of_elements; i++)
				data_physical[i] = 1;	// dummy data

			/*
			 * Spectral space data
			 */
			std::complex<double> *data_spectral = MemBlockAlloc::alloc< std::complex<double> >(spectral_array_data_number_of_elements*sizeof(std::complex<double>));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < spectral_array_data_number_of_elements; i++)
				data_spectral[i] = 1;	// dummy data

			fftw_plan_forward =
				fftw_plan_dft_r2c_2d(
					physical_data_size[1],	// n0 = ny
					physical_data_size[0],	// n1 = nx
					data_physical,
					(fftw_complex*)data_spectral,
					flags | FFTW_PRESERVE_INPUT
				);

			if (fftw_plan_forward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to get forward plan dft_r2c fftw" << std::endl;
				std::cerr << "r2c preverse_input forward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "FFTW-wisdom plan: rf" << physical_res[0] << "x" << physical_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans & TransformationPlans::REQUIRE_LOAD)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			fftw_plan_backward =
					fftw_plan_dft_c2r_2d(
						physical_res[1],	// n0 = ny
						physical_res[0],	// n1 = nx
						(fftw_complex*)data_spectral,
						data_physical,
						flags
					);

			if (fftw_plan_backward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to get backward plan dft_c2r fftw" << std::endl;
				std::cerr << "r2c backward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: rb" << physical_res[0] << "x" << physical_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans & TransformationPlans::REQUIRE_LOAD)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			MemBlockAlloc::free(data_physical, physical_array_data_number_of_elements*sizeof(double));
			MemBlockAlloc::free(data_spectral, spectral_array_data_number_of_elements*sizeof(std::complex<double>));

			// Backward scaling factor
			fftw_backward_scale_factor = 1.0/((double)(physical_data_size[0]*physical_data_size[1]));
		}


		/*
		 * COMPLEX PHYSICAL SPACE DATA
		 */
		{
			spectral_complex_data_size[0] = physical_data_size[0];
			spectral_complex_data_size[1] = physical_data_size[1];

			if ((spectral_complex_data_size[0] & 1) == 1)
				SWEETError("Not supported c");

			if ((spectral_complex_data_size[1] & 1) == 1)
				SWEETError("Not supported d");

			spectral_complex_ranges[0][0][0] = spectral_data_iteration_ranges[0][0][0];
			spectral_complex_ranges[0][0][1] = spectral_data_iteration_ranges[0][0][1];
			spectral_complex_ranges[0][1][0] = spectral_data_iteration_ranges[0][1][0];
			spectral_complex_ranges[0][1][1] = spectral_data_iteration_ranges[0][1][1];

			spectral_complex_ranges[1][0][0] = spectral_complex_ranges[0][0][0];
			spectral_complex_ranges[1][0][1] = spectral_complex_ranges[0][0][1];
			spectral_complex_ranges[1][1][0] = spectral_data_size[1] - spectral_complex_ranges[0][1][1] + 1;
			spectral_complex_ranges[1][1][1] = spectral_data_size[1];

			spectral_complex_ranges[2][0][0] = spectral_complex_data_size[0] - spectral_complex_ranges[0][0][1] + 1;
			spectral_complex_ranges[2][0][1] = spectral_complex_data_size[0];
			spectral_complex_ranges[2][1][0] = 0;
			spectral_complex_ranges[2][1][1] = spectral_complex_ranges[0][1][1];

			spectral_complex_ranges[3][0][0] = spectral_complex_data_size[0] - spectral_complex_ranges[0][0][1] + 1;
			spectral_complex_ranges[3][0][1] = spectral_complex_data_size[0];
			spectral_complex_ranges[3][1][0] = spectral_complex_data_size[1] - spectral_complex_ranges[0][1][1] + 1;
			spectral_complex_ranges[3][1][1] = spectral_complex_data_size[1];

			spectral_complex_array_data_number_of_elements = spectral_complex_data_size[0]*spectral_complex_data_size[1];

			/*
			 * Physical space data
			 */
			std::complex<double> *data_physical = MemBlockAlloc::alloc< std::complex<double> >(physical_array_data_number_of_elements*sizeof(std::complex<double>));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < physical_array_data_number_of_elements; i++)
				data_physical[i] = 1;	// dummy data

			/*
			 * Spectral space data
			 */
			std::complex<double> *data_spectral = MemBlockAlloc::alloc< std::complex<double> >(spectral_complex_array_data_number_of_elements*sizeof(std::complex<double>));

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < spectral_complex_array_data_number_of_elements; i++)
				data_spectral[i] = 1;	// dummy data


			fftw_plan_complex_forward =
					fftw_plan_dft_2d(
						physical_res[1],
						physical_res[0],
						(fftw_complex*)data_physical,
						(fftw_complex*)data_spectral,
						FFTW_FORWARD,
						flags
					);

			if (fftw_plan_complex_forward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to create complex forward plan for fftw" << std::endl;
				std::cerr << "complex forward preverse_input forward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << physical_res[0] << "x" << physical_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans == 2)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			fftw_plan_complex_backward =
					fftw_plan_dft_2d(
						physical_res[1],
						physical_res[0],
						(fftw_complex*)data_spectral,
						(fftw_complex*)data_physical,
						FFTW_BACKWARD,
						flags
					);

			if (fftw_plan_complex_backward == nullptr)
			{
				std::cerr << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
				std::cerr << "Failed to create complex backward plan for fftw" << std::endl;
				std::cerr << "complex backward preverse_input forward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << physical_res[0] << "x" << physical_res[1] << std::endl;
				if (i_reuse_spectral_transformation_plans == 2)
				{
					std::cerr << "********************************************************************************" << std::endl;
					std::cerr << "* IMPORTANT: Usage of FFTW wisdoms file enforced" << std::endl;
					std::cerr << "********************************************************************************" << std::endl;
				}
				exit(-1);
			}

			MemBlockAlloc::free(data_physical, physical_array_data_number_of_elements*sizeof(std::complex<double>));
			MemBlockAlloc::free(data_spectral, spectral_complex_array_data_number_of_elements*sizeof(std::complex<double>));
		}
#endif
	}



#if SWEET_USE_PLANE_SPECTRAL_SPACE

public:
	std::size_t get_spectral_iteration_range_area(int i)	const
	{
		return	(spectral_data_iteration_ranges[i][0][1] - spectral_data_iteration_ranges[i][0][0])*
				(spectral_data_iteration_ranges[i][1][1] - spectral_data_iteration_ranges[i][1][0]);
	}

#endif



#if SWEET_USE_LIBFFT
public:
	void fft_physical_to_spectral(
			double *i_physical_data,
			std::complex<double> *o_spectral_data
	)	const
	{
		fftw_execute_dft_r2c(
				fftw_plan_forward,
				i_physical_data,
				(fftw_complex*)o_spectral_data
			);
	}



	void fft_spectral_to_physical(
			std::complex<double> *i_spectral_data,
			double *o_physical_data
	)	const
	{
		fftw_execute_dft_c2r(
				fftw_plan_backward,
				(fftw_complex*)i_spectral_data,
				o_physical_data
			);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < physical_array_data_number_of_elements; i++)
			o_physical_data[i] *= fftw_backward_scale_factor;
	}



	void fft_complex_physical_to_spectral(
			std::complex<double> *i_physical_data,
			std::complex<double> *o_spectral_data
	)	const
	{
		fftw_execute_dft(
				fftw_plan_complex_forward,
				(fftw_complex*)i_physical_data,
				(fftw_complex*)o_spectral_data
			);
	}



	void fft_complex_spectral_to_physical(
			std::complex<double> *i_spectral_data,
			std::complex<double> *o_physical_data
	)	const
	{
		fftw_execute_dft(
				fftw_plan_complex_backward,
				(fftw_complex*)i_spectral_data,
				(fftw_complex*)o_physical_data
			);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < physical_array_data_number_of_elements; i++)
			o_physical_data[i] *= fftw_backward_scale_factor;
	}
#endif

public:
	void setup(
			int i_physical_res_x,
			int i_physical_res_y,

			int i_spectral_modes_x,
			int i_spectral_modes_y,

			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		physical_res[0] = i_physical_res_x;
		physical_res[1] = i_physical_res_y;

#if SWEET_USE_LIBFFT
		spectral_modes[0] = i_spectral_modes_x;
		spectral_modes[1] = i_spectral_modes_y;
#endif

		setup_internal_data(i_reuse_spectral_transformation_plans);
	}



public:
	void setupAuto(
			int io_physical_res[2],
			int io_spectral_modes[2],
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{

		if (io_physical_res[0] > 0 && io_spectral_modes[0] > 0)
		{
			setup(	io_physical_res[0],
					io_physical_res[1],
					io_spectral_modes[0],
					io_spectral_modes[1],
					i_reuse_spectral_transformation_plans
				);
			return;
		}

		if (io_physical_res[0] > 0)
		{
			setupAutoSpectralSpace(
					io_physical_res[0],
					io_physical_res[1],
					i_reuse_spectral_transformation_plans
				);

#if SWEET_USE_LIBFFT
			io_spectral_modes[0] = spectral_modes[0];
			io_spectral_modes[1] = spectral_modes[1];
#endif
			return;
		}

		if (io_spectral_modes[0] > 0)
		{
#if SWEET_USE_LIBFFT
			setupAutoPhysicalSpace(
					io_spectral_modes[0],
					io_spectral_modes[1],
					i_reuse_spectral_transformation_plans
				);

			io_physical_res[0] = physical_res[0];
			io_physical_res[1] = physical_res[1];
#else
			SWEETError("Setup with spectral modes not enabled");
#endif
			return;
		}

		SWEETError("No resolution/modes selected");
	}



public:
	void setupAutoSpectralSpace(
			int i_physical_res[2],
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		setupAutoSpectralSpace(
				i_physical_res[0],
				i_physical_res[1],
				i_reuse_spectral_transformation_plans
		);
	}



public:
	void setupAutoSpectralSpace(
			int i_physical_res_x,
			int i_physical_res_y,
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		physical_res[0] = i_physical_res_x;
		physical_res[1] = i_physical_res_y;

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING

		// REDUCTION IN EFFECTIVE SPECTRAL MODE RESOLUTION TO CUT OFF ANTI-ALIASED MODES
		spectral_modes[0] = (physical_res[0]*2)/3;
		spectral_modes[1] = (physical_res[1]*2)/3;


#else

	#if SWEET_USE_LIBFFT
		spectral_modes[0] = physical_res[0];
		spectral_modes[1] = physical_res[1];
	#endif

#endif

		setup_internal_data(i_reuse_spectral_transformation_plans);
	}


	void setupAutoSpectralSpace(
			int i_physical_res_x,
			int i_physical_res_y,
			int *o_spectral_res_x,
			int *o_spectral_res_y,
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		setupAutoSpectralSpace(i_physical_res_x, i_physical_res_y, i_reuse_spectral_transformation_plans);


#if SWEET_USE_LIBFFT
		*o_spectral_res_x = spectral_modes[0];
		*o_spectral_res_y = spectral_modes[1];
#endif
	}



#if SWEET_USE_LIBFFT

public:
	void setupAutoPhysicalSpace(
			int i_spectral_res_x,
			int i_spectral_res_y,
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		spectral_modes[0] = i_spectral_res_x;
		spectral_modes[1] = i_spectral_res_y;

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING

		// REDUCTION IN EFFECTIVE SPECTRAL MODE RESOLUTION TO CUT OFF ANTI-ALIASED MODES
		physical_res[0] = (spectral_modes[0]*3+1)/2;
		physical_res[1] = (spectral_modes[1]*3+1)/2;
#else

		physical_res[0] = spectral_modes[0];
		physical_res[1] = spectral_modes[1];

#endif

		setup_internal_data(i_reuse_spectral_transformation_plans);
	}


public:
	void setupAutoPhysicalSpace(
			int i_spectral_res_x,
			int i_spectral_res_y,
			int *o_physical_res_x,
			int *o_physical_res_y,
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		setupAutoPhysicalSpace(
				i_spectral_res_x,
				i_spectral_res_y,
				i_reuse_spectral_transformation_plans
		);

		*o_physical_res_x = physical_res[0];
		*o_physical_res_y = physical_res[1];
	}

	void setupAdditionalModes(
			PlaneDataConfig *i_planeConfig,
			int i_additional_modes_x,
			int i_additional_modes_y,
			TransformationPlans::TRANSFORMATION_PLAN_CACHE i_reuse_spectral_transformation_plans
	)
	{
		setupAutoPhysicalSpace(
				i_planeConfig->spectral_modes[0] + i_additional_modes_x,
				i_planeConfig->spectral_modes[1] + i_additional_modes_y,
				i_reuse_spectral_transformation_plans
		);
	}


#endif


	// TODO: CHECK THIS
	inline
	std::size_t getArrayIndexByModes(
			int n,
			int m
	)	const
	{
		assert(n >= 0);

		return n * spectral_data_size[0] + m;
	}


	// TODO: CHECK THIS
	inline
	std::size_t getArrayIndexByModes_Complex(
			int n,
			int m
	)	const
	{
		assert(n >= 0);

		int idx =  n * spectral_complex_data_size[0] + m;
		return idx;
	}



	void cleanup_data()
	{
		if (fftw_initialized)
		{
#if SWEET_USE_LIBFFT
			fftw_destroy_plan(fftw_plan_forward);
			fftw_destroy_plan(fftw_plan_backward);

			fftw_destroy_plan(fftw_plan_complex_forward);
			fftw_destroy_plan(fftw_plan_complex_backward);

			refCounterFftwPlans()--;
			assert(refCounterFftwPlans() >= 0);

			if (refCounterFftwPlans() == 0)
			{
				// backup wisdom
				if (reuse_spectral_transformation_plans == 1)
					storeWisdom();

#if SWEET_THREADING_SPACE
				fftw_cleanup_threads();
#endif
				fftw_cleanup();
			}
#endif
			fftw_initialized = false;
		}
	}


	void cleanup()
	{
		if (fftw_initialized)
		{
			cleanup_data();

			physical_res[0] = 0;
			physical_res[1] = 0;

	#if SWEET_USE_LIBFFT
			spectral_modes[0] = 0;
			spectral_modes[1] = 0;
	#endif
		}
	}



	~PlaneDataConfig()
	{
		cleanup();
	}
};



#endif
