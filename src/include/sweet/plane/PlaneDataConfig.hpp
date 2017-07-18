/*
 * PlaneDataConfig.hpp
 *
 *  Created on: 17 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <M.Schreiber@exeter.ac.uk>
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
#include <sweet/sweetmath.hpp>
#include <sweet/FatalError.hpp>



#if SWEET_THREADING
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
	std::size_t spectral_complex_data_iteration_ranges[4][2][2];

	fftw_plan	fftw_plan_complex_forward;
	fftw_plan	fftw_plan_complex_backward;

#endif

	bool initialized;


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
		std::cout << "spectral_complex_data_iteration_ranges [0][0]: " << spectral_complex_data_iteration_ranges[0][0][0] << ", " << spectral_complex_data_iteration_ranges[0][0][1] << std::endl;
		std::cout << "spectral_complex_data_iteration_ranges [0][1]: " << spectral_complex_data_iteration_ranges[0][1][0] << ", " << spectral_complex_data_iteration_ranges[0][1][1] << std::endl;
		std::cout << "spectral_complex_data_iteration_ranges [1][0]: " << spectral_complex_data_iteration_ranges[1][0][0] << ", " << spectral_complex_data_iteration_ranges[1][0][1] << std::endl;
		std::cout << "spectral_complex_data_iteration_ranges [1][1]: " << spectral_complex_data_iteration_ranges[1][1][0] << ", " << spectral_complex_data_iteration_ranges[1][1][1] << std::endl;
		std::cout << "spectral_complex_data_iteration_ranges [2][0]: " << spectral_complex_data_iteration_ranges[2][0][0] << ", " << spectral_complex_data_iteration_ranges[2][0][1] << std::endl;
		std::cout << "spectral_complex_data_iteration_ranges [2][1]: " << spectral_complex_data_iteration_ranges[2][1][0] << ", " << spectral_complex_data_iteration_ranges[2][1][1] << std::endl;
		std::cout << "spectral_complex_data_iteration_ranges [3][0]: " << spectral_complex_data_iteration_ranges[3][0][0] << ", " << spectral_complex_data_iteration_ranges[3][0][1] << std::endl;
		std::cout << "spectral_complex_data_iteration_ranges [3][1]: " << spectral_complex_data_iteration_ranges[3][1][0] << ", " << spectral_complex_data_iteration_ranges[3][1][1] << std::endl;
		std::cout << std::endl;
#endif
	}

public:
	static
	bool loadWisdom()
	{
		static const char *load_wisdom_from_file = nullptr;

		if (load_wisdom_from_file != nullptr)
			return false;

		load_wisdom_from_file = getenv("SWEET_FFTW_LOAD_WISDOM_FROM_FILE");

		if (load_wisdom_from_file == nullptr)
			return false;

		std::cout << "Loading SWEET_FFTW_LOAD_WISDOM_FROM_FILE=" << load_wisdom_from_file << std::endl;

		int wisdom_plan_loaded = fftw_import_wisdom_from_filename(load_wisdom_from_file);
		if (wisdom_plan_loaded == 0)
		{
			std::cerr << "Failed to load FFTW wisdom from file " << load_wisdom_from_file << std::endl;
			exit(1);
		}

		return true;
	}


public:
	int& refCounterFftwPlans()
	{
#if SWEET_THREADING
		if (omp_get_level() != 0)
			FatalError("PlaneDataConfig is not threadsafe, but called inside parallel region with more than one thread!!!");
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

		initialized = false;
	}



private:
	void setup_internal_data()
	{
		if (initialized)
		{
			//
			// refCounter()--;
			// is inside cleanup!
			cleanup();
		}
		else
		{
			// use REF counter in case of multiple plans
			// this allows a clean cleanup of fftw library
			initialized = true;
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
			FatalError("Lower physical resolution than spectral resolution not supported!");
		}

		assert(physical_res[0] > 0);
		assert(physical_res[1] > 0);

		assert(spectral_modes[0] > 0);
		assert(spectral_modes[1] > 0);


#if SWEET_REXI_THREAD_PARALLEL_SUM

		// Threaded time parallel sum
		std::cout << "Using REXI parallel sum, hence using only single FFT thread" << std::endl;

#else

	#if SWEET_THREADING
		// Is this the first instance?
		if (refCounterFftwPlans() == 1)
		{
			// initialise FFTW with spatial parallelization
			fftw_init_threads();
			fftw_plan_with_nthreads(omp_get_max_threads());
		}
	#endif

#endif


		/*
		 * Load existing wisdom
		 */
		bool wisdom_loaded = loadWisdom();


		int fftw_estimate_plan = 0;
		const char* fftw_estimate_plan_env = getenv("SWEET_FFTW_ESTIMATE");
		if (fftw_estimate_plan_env != nullptr)
		{
			std::cout << "Using estimated FFTW plan" << std::endl;
			fftw_estimate_plan = FFTW_ESTIMATE;
		}


		/*
		 * REAL PHYSICAL SPACE DATA (REAL to COMPLEX FFT)
		 */
		{
			// real-to-complex storage representation
			spectral_data_size[0] = physical_data_size[0]/2+1;
			spectral_data_size[1] = physical_data_size[1];


#if SWEET_USE_PLANE_SPECTRAL_DEALIASING

			/*
			 * For more information, have a look at
			 * doc/software_development_discussions/antialiasing/implementation_strategy.pdf
			 */

			int M = spectral_data_size[0];
			int N = spectral_data_size[1];

			if (	spectral_modes[0] == physical_res[0] ||
					spectral_modes[1] == physical_res[1]
			)
				FatalError("Aliasing doesn't make sense since physical resolution is identical to spectral");

			spectral_data_iteration_ranges[0][0][0] = 0;
			spectral_data_iteration_ranges[0][0][1] = 2*(M-1)/3+1;
			spectral_data_iteration_ranges[0][1][0] = 0;
			spectral_data_iteration_ranges[0][1][1] = N/3;

			spectral_data_iteration_ranges[1][0][0] = 0;
			spectral_data_iteration_ranges[1][0][1] = 2*(M-1)/3+1;
			spectral_data_iteration_ranges[1][1][0] = N-N/3+1;
			spectral_data_iteration_ranges[1][1][1] = N;

			spectral_data_iteration_ranges[0][0][1]--;
			spectral_data_iteration_ranges[1][0][1]--;

			spectral_real_modes[0] = spectral_data_iteration_ranges[0][0][1];
			spectral_real_modes[1] = spectral_data_iteration_ranges[0][1][1];

#else

			spectral_data_iteration_ranges[0][0][0] = 0;
			spectral_data_iteration_ranges[0][0][1] = spectral_data_size[0];	// padding
			spectral_data_iteration_ranges[0][1][0] = 0;
			spectral_data_iteration_ranges[0][1][1] = spectral_data_size[1]/2;

			spectral_data_iteration_ranges[1][0][0] = 0;
			spectral_data_iteration_ranges[1][0][1] = spectral_data_size[0];	// padding
			spectral_data_iteration_ranges[1][1][0] = spectral_data_size[1]/2;
			spectral_data_iteration_ranges[1][1][1] = spectral_data_size[1];


			spectral_real_modes[0] = spectral_data_iteration_ranges[0][0][1];
			spectral_real_modes[1] = spectral_data_iteration_ranges[0][1][1];
#endif

			spectral_array_data_number_of_elements = spectral_data_size[0]*spectral_data_size[1];



			/*
			 * Physical space data
			 */
			double *data_physical = MemBlockAlloc::alloc<double>(physical_array_data_number_of_elements*sizeof(double));
			//std::cout << physical_array_data_number_of_elements << std::endl;
			//std::cout << spectral_array_data_number_of_elements << std::endl;

	#if SWEET_THREADING
	#pragma omp parallel for OPENMP_PAR_SIMD
	#endif
			for (std::size_t i = 0; i < physical_array_data_number_of_elements; i++)
				data_physical[i] = 1;	// dummy data

			/*
			 * Spectral space data
			 */
			std::complex<double> *data_spectral = MemBlockAlloc::alloc< std::complex<double> >(spectral_array_data_number_of_elements*sizeof(std::complex<double>));

	#if SWEET_THREADING
	#pragma omp parallel for OPENMP_PAR_SIMD
	#endif
			for (std::size_t i = 0; i < spectral_array_data_number_of_elements; i++)
				data_spectral[i] = 1;	// dummy data


			fftw_plan_forward =
					fftw_plan_dft_r2c_2d(
						physical_data_size[1],	// n0 = ny
						physical_data_size[0],	// n1 = nx
						data_physical,
						(fftw_complex*)data_spectral,
						//(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
						(!wisdom_loaded ? fftw_estimate_plan | FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
					);

			if (fftw_plan_forward == nullptr)
			{
				std::cerr << "Failed to create forward plan for fftw" << std::endl;
				std::cerr << "r2c preverse_input forward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "FFTW-wisdom plan: rf" << physical_res[0] << "x" << physical_res[1] << std::endl;
				exit(-1);
			}

			fftw_plan_backward =
					fftw_plan_dft_c2r_2d(
						physical_res[1],	// n0 = ny
						physical_res[0],	// n1 = nx
						(fftw_complex*)data_spectral,
						data_physical,
						(!wisdom_loaded ? fftw_estimate_plan : FFTW_WISDOM_ONLY)
					);

			if (fftw_plan_backward == nullptr)
			{
				std::cerr << "Failed to create backward plan for fftw" << std::endl;
				std::cerr << "r2c backward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: rb" << physical_res[0] << "x" << physical_res[1] << std::endl;
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

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING && 0

			/*
			 * For more information, have a look at
			 * doc/antialiasing/implementation_strategy.pdf
			 */

			spectral_complex_data_iteration_ranges[0][0][0] = 0;
			spectral_complex_data_iteration_ranges[0][0][1] = spectral_complex_data_size[0]*2/3-1;
			spectral_complex_data_iteration_ranges[0][1][0] = 0;
			spectral_complex_data_iteration_ranges[0][1][1] = spectral_complex_data_size[1]/3;

			spectral_complex_data_iteration_ranges[1][0][0] = 0;
			spectral_complex_data_iteration_ranges[1][0][1] = spectral_complex_data_size[0]*2/3-1;
			spectral_complex_data_iteration_ranges[1][1][0] = spectral_complex_data_size[1]-spectral_complex_data_size[1]/3;
			spectral_complex_data_iteration_ranges[1][1][1] = spectral_complex_data_size[1];

#warning "TODO: setup correct ranges"
			spectral_complex_data_iteration_ranges[2][0][0] = spectral_complex_data_size[0] - (spectral_complex_data_size[0]*2/3-1);	// TODO: check this start index
			spectral_complex_data_iteration_ranges[2][0][1] = spectral_complex_data_size[0];
			spectral_complex_data_iteration_ranges[2][1][0] = 0;
			spectral_complex_data_iteration_ranges[2][1][1] = spectral_complex_data_size[1]/3;

			spectral_complex_data_iteration_ranges[3][0][0] = spectral_complex_data_size[0] - (spectral_complex_data_size[0]*2/3-1);	// TODO: check this start index
			spectral_complex_data_iteration_ranges[3][0][1] = spectral_complex_data_size[0];
			spectral_complex_data_iteration_ranges[3][1][0] = spectral_complex_data_size[1]-spectral_complex_data_size[1]/3;
			spectral_complex_data_iteration_ranges[3][1][1] = spectral_complex_data_size[1];

#else

			spectral_complex_data_iteration_ranges[0][0][0] = 0;
			spectral_complex_data_iteration_ranges[0][0][1] = spectral_complex_data_size[0]/2;
			spectral_complex_data_iteration_ranges[0][1][0] = 0;
			spectral_complex_data_iteration_ranges[0][1][1] = spectral_complex_data_size[1]/2;

			spectral_complex_data_iteration_ranges[1][0][0] = 0;
			spectral_complex_data_iteration_ranges[1][0][1] = spectral_complex_data_size[0]/2;
			spectral_complex_data_iteration_ranges[1][1][0] = spectral_complex_data_size[1]/2;
			spectral_complex_data_iteration_ranges[1][1][1] = spectral_complex_data_size[1];

			spectral_complex_data_iteration_ranges[2][0][0] = spectral_complex_data_size[0]/2;
			spectral_complex_data_iteration_ranges[2][0][1] = spectral_complex_data_size[0];
			spectral_complex_data_iteration_ranges[2][1][0] = 0;
			spectral_complex_data_iteration_ranges[2][1][1] = spectral_complex_data_size[1]/2;

			spectral_complex_data_iteration_ranges[3][0][0] = spectral_complex_data_size[0]/2;
			spectral_complex_data_iteration_ranges[3][0][1] = spectral_complex_data_size[0];
			spectral_complex_data_iteration_ranges[3][1][0] = spectral_complex_data_size[1]/2;
			spectral_complex_data_iteration_ranges[3][1][1] = spectral_complex_data_size[1];

#endif

			spectral_complex_array_data_number_of_elements = spectral_complex_data_size[0]*spectral_complex_data_size[1];

			/*
			 * Physical space data
			 */
			std::complex<double> *data_physical = MemBlockAlloc::alloc< std::complex<double> >(physical_array_data_number_of_elements*sizeof(std::complex<double>));

	#if SWEET_THREADING
	#pragma omp parallel for OPENMP_PAR_SIMD
	#endif
			for (std::size_t i = 0; i < physical_array_data_number_of_elements; i++)
				data_physical[i] = 1;	// dummy data

			/*
			 * Spectral space data
			 */
			std::complex<double> *data_spectral = MemBlockAlloc::alloc< std::complex<double> >(spectral_complex_array_data_number_of_elements*sizeof(std::complex<double>));

	#if SWEET_THREADING
	#pragma omp parallel for OPENMP_PAR_SIMD
	#endif
			for (std::size_t i = 0; i < spectral_complex_array_data_number_of_elements; i++)
				data_spectral[i] = 1;	// dummy data


			fftw_plan_complex_forward =
					fftw_plan_dft_2d(
						physical_res[1],
						physical_res[0],
						(fftw_complex*)data_physical,
						(fftw_complex*)data_spectral,
						FFTW_FORWARD,
//						(!wisdom_loaded ? FFTW_PRESERVE_INPUT | fftw_estimate_plan : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
						(!wisdom_loaded ? fftw_estimate_plan : FFTW_WISDOM_ONLY)
					);

			if (fftw_plan_complex_forward == nullptr)
			{
				std::cerr << "Failed to create plan_forward for fftw" << std::endl;
				std::cerr << "complex forward preverse_input forward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << physical_res[0] << "x" << physical_res[1] << std::endl;
				exit(-1);
			}

			fftw_plan_complex_backward =
					fftw_plan_dft_2d(
						physical_res[1],
						physical_res[0],
						(fftw_complex*)data_spectral,
						(fftw_complex*)data_physical,
						FFTW_BACKWARD,
//						(!wisdom_loaded ? FFTW_PRESERVE_INPUT | fftw_estimate_plan : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
						(!wisdom_loaded ? fftw_estimate_plan : FFTW_WISDOM_ONLY)
					);

			if (fftw_plan_complex_backward == nullptr)
			{
				std::cerr << "Failed to create plan_backward for fftw" << std::endl;
				std::cerr << "complex backward preverse_input forward " << physical_res[0] << " x " << physical_res[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << physical_res[0] << "x" << physical_res[1] << std::endl;
				exit(-1);
			}


			MemBlockAlloc::free(data_physical, physical_array_data_number_of_elements*sizeof(std::complex<double>));
			MemBlockAlloc::free(data_spectral, spectral_complex_array_data_number_of_elements*sizeof(std::complex<double>));
		}
#endif


#if 0
		if (((spectral_modes[0] & 1) != 0) || ((spectral_modes[1] & 1) != 0))
		{
//			FatalError("Only even number of spectral modes are supported!");
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

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
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



	void fft_spectral_to_complex_physical(
			std::complex<double> *i_spectral_data,
			std::complex<double> *o_physical_data
	)	const
	{
		fftw_execute_dft(
				fftw_plan_complex_backward,
				(fftw_complex*)i_spectral_data,
				(fftw_complex*)o_physical_data
			);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < physical_array_data_number_of_elements; i++)
			o_physical_data[i] *= fftw_backward_scale_factor;
	}
#endif

public:
	void setup(
			int i_physical_res_x,
			int i_physical_res_y,

			int i_spectral_modes_x,
			int i_spectral_modes_y
	)
	{
		physical_res[0] = i_physical_res_x;
		physical_res[1] = i_physical_res_y;

#if SWEET_USE_LIBFFT
		spectral_modes[0] = i_spectral_modes_x;
		spectral_modes[1] = i_spectral_modes_y;
#endif

		setup_internal_data();
	}



public:
	void setupAuto(
			int io_physical_res[2],
			int io_spectral_modes[2]
	)
	{
//		std::cout << io_physical_res[0] << ", " << io_physical_res[1] << std::endl;
//		std::cout << io_spectral_modes[0] << ", " << io_spectral_modes[1] << std::endl;

		if (io_physical_res[0] > 0 && io_spectral_modes[0] > 0)
		{
			setup(	io_physical_res[0],
					io_physical_res[1],
					io_spectral_modes[0],
					io_spectral_modes[1]
				);
			return;
		}

		if (io_physical_res[0] > 0)
		{
			setupAutoSpectralSpace(
					io_physical_res[0],
					io_physical_res[1]
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
					io_spectral_modes[1]
				);

			io_physical_res[0] = physical_res[0];
			io_physical_res[1] = physical_res[1];
#else
			FatalError("Setup with spectral modes not enabled");
#endif
			return;
		}

		FatalError("No resolution/modes selected");
	}



public:
	void setupAutoSpectralSpace(
			int i_physical_res[2]
	)
	{
		setupAutoSpectralSpace(i_physical_res[0], i_physical_res[1]);
	}



public:
	void setupAutoSpectralSpace(
			int i_physical_res_x,
			int i_physical_res_y
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

		setup_internal_data();
	}


	void setupAutoSpectralSpace(
			int i_physical_res_x,
			int i_physical_res_y,
			int *o_spectral_res_x,
			int *o_spectral_res_y
	)
	{
		setupAutoSpectralSpace(i_physical_res_x, i_physical_res_y);


#if SWEET_USE_LIBFFT
		*o_spectral_res_x = spectral_modes[0];
		*o_spectral_res_y = spectral_modes[1];
#endif
	}



#if SWEET_USE_LIBFFT

public:
	void setupAutoPhysicalSpace(
			int i_spectral_res_x,
			int i_spectral_res_y
	)
	{
		spectral_modes[0] = i_spectral_res_x;
		spectral_modes[1] = i_spectral_res_y;

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING

		// REDUCTION IN EFFECTIVE SPECTRAL MODE RESOLUTION TO CUT OFF ANTI-ALIASED MODES
		// TODO: check for correct anti-aliasing rule
		physical_res[0] = (spectral_modes[0]*3+1)/2;
		physical_res[1] = (spectral_modes[1]*3+1)/2;

#else

	#if SWEET_USE_LIBFFT
		physical_res[0] = spectral_modes[0];
		physical_res[1] = spectral_modes[1];
	#endif

#endif

		setup_internal_data();
	}


public:
	void setupAutoPhysicalSpace(
			int i_spectral_res_x,
			int i_spectral_res_y,
			int *o_physical_res_x,
			int *o_physical_res_y
	)
	{
		setupAutoPhysicalSpace(
				i_spectral_res_x,
				i_spectral_res_y
		);

		*o_physical_res_x = physical_res[0];
		*o_physical_res_y = physical_res[1];
	}

	void setupAdditionalModes(
			PlaneDataConfig *i_planeConfig,
			int i_additional_modes_x,
			int i_additional_modes_y
	)
	{
		setupAutoPhysicalSpace(
				i_planeConfig->spectral_modes[0] + i_additional_modes_x,
				i_planeConfig->spectral_modes[1] + i_additional_modes_y
		);
	}


#endif



	void cleanup()
	{
		if (initialized)
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
#if SWEET_THREADING
				fftw_cleanup_threads();
#endif
				fftw_cleanup();
			}
#endif
		}
	}



	~PlaneDataConfig()
	{
		cleanup();
	}
};



#endif
