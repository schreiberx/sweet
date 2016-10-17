/*
 * PlaneData.hpp
 *
 *  Created on: 28 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_PLANE_DATA_HPP_
#define SRC_PLANE_DATA_HPP_

#include <complex>
#include <cassert>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <memory>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <utility>
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sweet/sweetmath.hpp>
#include <sweet/openmp_helper.hpp>
#include <sweet/MemBlockAlloc.hpp>

#ifndef SWEET_USE_PLANE_SPECTRAL_SPACE
	#define SWEET_USE_PLANE_SPECTRAL_SPACE	1
#endif

#ifndef SWEET_USE_PLANE_SPECTRAL_DEALIASING
	#define SWEET_USE_PLANE_SPECTRAL_DEALIASING 1
#endif

#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_LIBFFT
#	include <fftw3.h>
#endif

#if SWEET_THREADING
#	include <omp.h>
#endif



/**
 * Plane data and operator support.
 *
 * Here, we assume the Cartesian coordinate system given similar to the following sketch:
 *
 *  Y ^
 *    |
 *    |
 *    |
 *    +-------->
 *           X
 *
 * Also the arrays are stored in this way:
 * 		A[Y0...YN-1][X0...XN-1]
 *
 */
class PlaneData
{
public:
	/**
	 * global size of allocated array
	 * (x,y[,z])
	 */
	std::size_t resolution[2];

	/**
	 * local data in cartesian space
	 */
	std::size_t array_data_cartesian_length;
	double *array_data_cartesian_space;

#if SWEET_USE_LIBFFT || SWEET_USE_PLANE_SPECTRAL_SPACE
	std::size_t resolution_spec[2];
	bool array_data_cartesian_space_valid;

	/**
	 * local data in spectral space
	 */
	std::size_t array_data_spectral_length;
	double *array_data_spectral_space;
	bool array_data_spectral_space_valid;

	bool aliasing_scaled;
#endif

	/**
	 * local ranges
	 */
	std::size_t range_start[2];
	std::size_t range_end[2];
	std::size_t range_size[2];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	std::size_t range_spec_start[2];
	std::size_t range_spec_end[2];
	std::size_t range_spec_size[2];
#else
	int kernel_size = -1;
	double *kernel_data = nullptr;
	int kernel_id = -1;
#endif

	/**
	 * temporary data?
	 *
	 * Temporary data can be created if e.g. operators are evaluated:
	 * h = hu/h
	 *
	 * This first creates a temporary PlaneData to compute hu/h.
	 *
	 * This is then followed by an assignment of this data to h.
	 */
//	bool temporary_data;


#if 0
private:
	// http://stackoverflow.com/questions/124856/how-do-i-prevent-a-class-from-being-allocated-via-the-new-operator-id-like
	// Prevent heap allocation
	void *operator new(std::size_t);
	void *operator new[](std::size_t);
	void operator delete(void*);
	void operator delete[](void*);
#endif

#if 0
	/**
	 * allocator which allocated memory blocks aligned at 128 byte boundaries
	 */
public:
	template <typename T=void>
	static
	T *alloc_aligned_mem(
			std::size_t i_size
	)
	{
		T *data;
		int retval = posix_memalign((void**)&data, 128, i_size);
		if (retval != 0)
		{
			std::cerr << "Unable to allocate memory" << std::endl;
			assert(false);
			exit(-1);
		}
		return data;
	}


	/**
	 * Free memory which was previously allocated
	 */
public:
	template <typename T>
	void free(T *i_ptr)
	{
		::free(i_ptr);
	}
#endif



	/**
	 * allocate buffers
	 *
	 * The size is given in array_data_cartesian_length and array_data_spectral_length
	 */
private:
	inline
	void p_allocate_buffers(
			bool i_first_touch_initialize = true	///< true: initialize the data buffers with dummy data for first touch policy of page allocation on shared-memory systems
	)
	{
		MemBlockAlloc::free(array_data_cartesian_space, array_data_cartesian_length*sizeof(double));
		array_data_cartesian_space = MemBlockAlloc::alloc<double>(array_data_cartesian_length*sizeof(double));

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		MemBlockAlloc::free(array_data_spectral_space, array_data_cartesian_length*sizeof(double));
		array_data_spectral_space = MemBlockAlloc::alloc<double>(array_data_spectral_length*sizeof(double));
#endif

		if (i_first_touch_initialize)
		{
			// use parallel setup for first touch policy!
#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				array_data_cartesian_space[i] = -12345;	// dummy data

#if SWEET_USE_PLANE_SPECTRAL_SPACE
#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				array_data_spectral_space[i] = -12345;	// dummy data
#endif
		}
	}



	/**
	 * request that the buffers are allocated with the specified resolution
	 */
private:
	void p_request_buffers_with_resolution(
		const std::size_t i_resolution[2]		///< requested resolution
	)
	{
		// check if the resolution is maybe already the same as requested
		int i = 0;
		for (; i < 2; i++)
			if (resolution[i] != i_resolution[i])
				break;

		if (i == 2)
			return;

		array_data_cartesian_length = 1;
		for (int i = 0; i < 2; i++)
		{
			array_data_cartesian_length *= i_resolution[i];

			resolution[i] = i_resolution[i];
			range_start[i] = 0;
			range_end[i] = i_resolution[i];
			range_size[i] = range_end[i]-range_start[i];
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = false;

		array_data_spectral_length = array_data_cartesian_length/i_resolution[0];	/// see FFTW documentation for allocation of memory buffers
		array_data_spectral_length *= 2*(i_resolution[0]/2+1);

		for (int i = 0; i < 2; i++)
		{
			if (i == 0)
				resolution_spec[0] = (resolution[0]/2+1);
			else
				resolution_spec[i] = resolution[i];

			range_spec_start[i] = 0;
			range_spec_end[i] = resolution_spec[i];
			range_spec_size[i] = resolution_spec[i];
		}

		array_data_spectral_space_valid = false;
#endif

		p_allocate_buffers(false);
	}


#if SWEET_USE_LIBFFT || SWEET_USE_PLANE_SPECTRAL_SPACE

public:
	class FFTWSingletonClass
	{
	public:
		int ref_counter;

	private:
		fftw_plan	plan_forward;
		std::size_t plan_forward_output_length;

		fftw_plan	plan_backward;
		std::size_t plan_backward_output_length;

	public:
		static
		bool loadWisdom()
		{

			// only initialize it, if this class is not called to initialize aliasing
#if SWEET_REXI_THREAD_PARALLEL_SUM
			std::cout << "Using REXI parallel sum, hence using only single FFT thread" << std::endl;
			// only use serial FFT in case of REXI parallel sum
			// automatically use single-threaded fftw only
//#if SWEET_THREADING
//#	error "Incompatible/unsupported compile flags detected"
//#endif
			//fftw_init_threads();
			//fftw_plan_with_nthreads(1);
#else

	#if SWEET_THREADING
			// support threading
			fftw_init_threads();
			fftw_plan_with_nthreads(omp_get_max_threads());
	#endif

#endif
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
		FFTWSingletonClass(
				PlaneData &i_dataArray
		)	:
			ref_counter(0)
		{
			plan_backward_output_length = i_dataArray.array_data_cartesian_length;
			plan_forward_output_length = i_dataArray.array_data_spectral_length;

			double *data_cartesian = MemBlockAlloc::alloc<double>(i_dataArray.array_data_cartesian_length*sizeof(double));
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < i_dataArray.array_data_cartesian_length; i++)
				data_cartesian[i] = -123;	// dummy data

			double *data_spectral = MemBlockAlloc::alloc<double>(i_dataArray.array_data_spectral_length*sizeof(double));
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < i_dataArray.array_data_spectral_length; i++)
				data_spectral[i] = -123;	// dummy data

			// allow to search for a good plan only for 60 seconds
//			fftw_set_timelimit(60);


			/*
			 * FFT WISDOM INIT
			 */
			bool wisdom_loaded = loadWisdom();

			plan_forward =
					fftw_plan_dft_r2c_2d(
						i_dataArray.resolution[1],	// n0 = ny
						i_dataArray.resolution[0],	// n1 = nx
						data_cartesian,
						(fftw_complex*)data_spectral,
//						FFTW_PRESERVE_INPUT
						(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
					);

//			fftw_print_plan(plan_forward);
			if (plan_forward == nullptr)
			{
				std::cerr << "Failed to create forward plan for fftw" << std::endl;
				std::cerr << "r2c preverse_input forward " << i_dataArray.resolution[0] << " x " << i_dataArray.resolution[1] << std::endl;
				std::cerr << "fftw-wisdom plan: rf" << i_dataArray.resolution[0] << "x" << i_dataArray.resolution[1] << std::endl;
				exit(-1);
			}

			plan_backward =
					fftw_plan_dft_c2r_2d(
						i_dataArray.resolution[1],	// n0 = ny
						i_dataArray.resolution[0],	// n1 = nx
						(fftw_complex*)data_spectral,
						data_cartesian,
						(!wisdom_loaded ? 0 : FFTW_WISDOM_ONLY)
					);

			if (plan_backward == nullptr)
			{
				std::cerr << "Failed to create backward plan for fftw" << std::endl;
				std::cerr << "r2c backward " << i_dataArray.resolution[0] << " x " << i_dataArray.resolution[1] << std::endl;
				std::cerr << "fftw-wisdom plan: rb" << i_dataArray.resolution[0] << "x" << i_dataArray.resolution[1] << std::endl;
				exit(-1);
			}

			// always store plans - maybe they got extended with another one
//			if (wisdom_plan_loaded == 0)
#if 0
			fftw_export_wisdom_to_filename(fftw_load_wisdom_filename);
#endif

			MemBlockAlloc::free(data_cartesian, i_dataArray.array_data_cartesian_length*sizeof(double));
			MemBlockAlloc::free(data_spectral, i_dataArray.array_data_spectral_length*sizeof(double));
		}

		void fft_forward(
				PlaneData &io_dataArray
		)
		{
			fftw_execute_dft_r2c(plan_forward, io_dataArray.array_data_cartesian_space, (fftw_complex*)io_dataArray.array_data_spectral_space);
		}

		void fft_backward(
				PlaneData &io_dataArray
		)
		{
			fftw_execute_dft_c2r(plan_backward, (fftw_complex*)io_dataArray.array_data_spectral_space, io_dataArray.array_data_cartesian_space);
			// spectral data is not valid anymore, since c2r is destructive!
			io_dataArray.array_data_spectral_space_valid = false;

			double scale = (1.0/(double)plan_backward_output_length);
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < plan_backward_output_length; i++)
				io_dataArray.array_data_cartesian_space[i] *= scale;
		}

		~FFTWSingletonClass()
		{
			fftw_destroy_plan(plan_forward);
			fftw_destroy_plan(plan_backward);

#if SWEET_THREADING
			fftw_cleanup_threads();
#endif

			fftw_cleanup();
		}
	};
#endif



	/**
	 * prohibit empty initialization by making this method private
	 */
private:
	PlaneData()
#if 1
	{}
#else
	{
		resolution[0] = 0;
		resolution[1] = 0;

		temporary_data = false;

		array_data_cartesian_space = nullptr;
		array_data_cartesian_length = 0;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = false;

		array_data_spectral_space = nullptr;
		array_data_spectral_space_valid = false;
		array_data_spectral_length = 0;
#endif
	}
#endif




public:
	/**
	 * copy constructor, used e.g. in
	 * 	PlaneData tmp_h = h;
	 * 	PlaneData tmp_h2(h);
	 *
	 * Duplicate all data
	 */
	PlaneData(
			const PlaneData &i_dataArray
	)	:
		array_data_cartesian_space(nullptr)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		,array_data_spectral_space(nullptr),
		aliasing_scaled(i_dataArray.aliasing_scaled)
#endif

//		,temporary_data(false)
	{
		for (int i = 0; i < 2; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			resolution_spec[i] = i_dataArray.resolution_spec[i];

			range_spec_start[i] = i_dataArray.range_spec_start[i];
			range_spec_end[i] = i_dataArray.range_spec_end[i];
			range_spec_size[i] = i_dataArray.range_spec_size[i];
#endif
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_spectral_length = i_dataArray.array_data_spectral_length;
#endif

		p_allocate_buffers(false);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = i_dataArray.array_data_cartesian_space_valid;
		if (array_data_cartesian_space_valid)
#endif
		{
			// use parallel copy for first touch policy!
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				array_data_cartesian_space[i] = i_dataArray.array_data_cartesian_space[i];
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_spectral_length = i_dataArray.array_data_spectral_length;
		array_data_spectral_space_valid = i_dataArray.array_data_spectral_space_valid;

		if (array_data_spectral_space_valid)
		{
			// use parallel copy for first touch policy!
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				array_data_spectral_space[i] = i_dataArray.array_data_spectral_space[i];
		}

		if (aliasing_scaled)
		{
			FFTWSingletonClass* fft_ptr = *fftAliasingGetSingletonPtr();
			fft_ptr->ref_counter++;
		}

		{
			FFTWSingletonClass* fft_ptr = *fftGetSingletonPtr();
			fft_ptr->ref_counter++;
		}
#endif
	}



public:
	/**
	 * Move constructor
	 */
	PlaneData(
			PlaneData &&i_dataArray
	)	:
		array_data_cartesian_space(nullptr)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		,array_data_spectral_space(nullptr),
		aliasing_scaled(i_dataArray.aliasing_scaled)
#endif
//		,temporary_data(false)
	{
		for (int i = 0; i < 2; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			resolution_spec[i] = i_dataArray.resolution_spec[i];

			range_spec_start[i] = i_dataArray.range_spec_start[i];
			range_spec_end[i] = i_dataArray.range_spec_end[i];
			range_spec_size[i] = i_dataArray.range_spec_size[i];
#endif
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;
		array_data_cartesian_space = i_dataArray.array_data_cartesian_space;
		i_dataArray.array_data_cartesian_space = nullptr;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = i_dataArray.array_data_cartesian_space_valid;
		i_dataArray.array_data_cartesian_space_valid = false;

		array_data_spectral_length = i_dataArray.array_data_spectral_length;
		array_data_spectral_space = i_dataArray.array_data_spectral_space;
		i_dataArray.array_data_spectral_space = nullptr;
		array_data_spectral_space_valid = i_dataArray.array_data_spectral_space_valid;
		i_dataArray.array_data_spectral_space_valid = false;

		if (aliasing_scaled)
		{
			FFTWSingletonClass* fft_ptr = *fftAliasingGetSingletonPtr();
			fft_ptr->ref_counter++;
		}

		{
			FFTWSingletonClass* fft_ptr = *fftGetSingletonPtr();
			fft_ptr->ref_counter++;
		}
#endif
	}


#if 0
	/**
	 * special empty constructor
	 *
	 * note, that the class may only be used after calling the setup() method.
	 */
public:
	PlaneData(
		int i_int
	)	:
		array_data_cartesian_space(nullptr),
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_spectral_space(nullptr),
		aliasing_scaled(false),
#endif
		temporary_data(false)
	{
	}
#endif


public:
	inline
	void checkConsistency(bool debug = false)	const
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE && SWEET_DEBUG_MODE==1 && 0

		static bool volatile inCheck = false;

		bool shouldReturn = false;

#pragma omp critical
		{
			if (inCheck)
				shouldReturn = true;
			else
				inCheck = true;
		}


		if (shouldReturn)
			return;

		if (array_data_cartesian_space_valid == true && array_data_spectral_space_valid == true)
		{
			{
				PlaneData a(*this);
				a.array_data_cartesian_space_valid = false;
				a.requestDataInCartesianSpace();

				if (debug)
				{
					std::cout << "DEBUG: PRINT CART DATA" << std::endl;
					a.printArrayData();
				}

				double sum = 0;
				for (std::size_t i = 0; i < array_data_cartesian_length; i++)
					sum += std::abs(a.array_data_cartesian_space[i] - array_data_cartesian_space[i]);

				sum /= (double)array_data_cartesian_length;

				if (sum > 1e-6)
				{
					std::cout << std::endl;
					std::cout << "CARTESIAN COMPARISON: Inconsistent data detected" << std::endl;
					std::cout << std::endl;
					std::cout << "************************************" << std::endl;
					std::cout << "Cart data:" << std::endl;
					printArrayData();
					std::cout << std::endl;
					std::cout << "************************************" << std::endl;

					std::cout << "************************************" << std::endl;
					std::cout << "Spec data:" << std::endl;
					printSpectrum();
					std::cout << "************************************" << std::endl;

					assert(false);
					exit(1);
				}
			}


			{
				PlaneData a(*this);
				a.array_data_spectral_space_valid = false;
				a.requestDataInSpectralSpace();


				double sum_re = 0;
				double sum_im = 0;
				for (std::size_t i = 0; i < array_data_cartesian_length/2; i++)
				{
					sum_re += std::abs(a.array_data_spectral_space[2*i+0] - array_data_spectral_space[2*i+0]);
					sum_im += std::abs(a.array_data_spectral_space[2*i+1] - array_data_spectral_space[2*i+1]);
				}

				sum_re /= (double)(array_data_spectral_length/2);
				sum_im /= (double)(array_data_spectral_length/2);

				if (sum_re > 1e-6|| sum_im > 1e-6)
				{
					std::cout << "************************************" << std::endl;
					std::cout << "SPECTRUM COMPARISON: Inconsistent data detected" << std::endl;
					std::cout << "sumerr: " << sum_re << ", " << sum_im << std::endl;
					a.printSpectrum();
					std::cout << std::endl;
					std::cout << "************************************" << std::endl;
					std::cout << "Cart data:" << std::endl;
					printArrayData();
					std::cout << std::endl;
					std::cout << "************************************" << std::endl;

					std::cout << "************************************" << std::endl;
					std::cout << "Spec data:" << std::endl;
					printSpectrum();
					std::cout << "************************************" << std::endl;

					assert(false);
					exit(1);
				}
			}
		}

		inCheck = false;

//		assert(!(array_data_cartesian_space_valid == true && array_data_spectral_space_valid == true));
#endif
	}


public:
	/**
	 * setup the PlaneData in case that the special
	 * empty constructor with int as a parameter was used.
	 *
	 * Calling this setup function should be in general avoided.
	 */
public:
	void setup(
			const std::size_t i_resolution[2],	///< size of array
			bool i_anti_aliasing = false		///< set to true if this should be used as an anti-aliasing dataArray
	)
	{
		assert(array_data_cartesian_space == nullptr);

		for (int i = 0; i < 2; i++)
		{
			resolution[i] = 0;
#if SWEET_USE_PLANE_SPECTRAL_SPACE
			if (i_resolution[i] & 1)
			{
				std::cerr << "Sorry - only even resolution supported for spectral space, makes the life much easier!" << std::endl;
				exit(1);
			}
#endif
		}

		p_request_buffers_with_resolution(i_resolution);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		// initialize fft if not yet done
		fftTestAndInit(*this);

		aliasing_scaled = i_anti_aliasing;
		if (i_anti_aliasing)
			fftAliasingTestAndInit(*this);
#endif

		checkConsistency();
	}


	/**
	 * default constructor
	 */
public:
	PlaneData(
		const std::size_t i_resolution[2],	///< size of array for each dimension,
		bool i_anti_aliasing = false		///< set to true if this should be used as an anti-aliasing dataArray
	)	:
		array_data_cartesian_space(nullptr)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		,array_data_spectral_space(nullptr),
		aliasing_scaled(false)
#endif
//		,temporary_data(false)
	{
		setup(i_resolution, i_anti_aliasing);
	}



	~PlaneData()
	{
		MemBlockAlloc::free(array_data_cartesian_space, array_data_cartesian_length*sizeof(double));

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		MemBlockAlloc::free(array_data_spectral_space, array_data_spectral_length*sizeof(double));

		/*
		 * If this is an aliasing PlaneData, reduce the reference counter of this array
		 */
		if (aliasing_scaled)
		{
			FFTWSingletonClass* fft_anti_ptr = *fftAliasingGetSingletonPtr();
			if (fft_anti_ptr == nullptr)
			{
				std::cerr << "FFTS for aliasing PlaneDatas still existing, but no FFTW handlers!" << std::endl;
				exit(1);
			}
			fft_anti_ptr->ref_counter--;

			/*
			 * NOTE: The aliasing plans are not free'd here. See below for the reason
			 */
		}

		{
			FFTWSingletonClass* fft_ptr = *fftGetSingletonPtr();
			assert(fft_ptr != nullptr);
			fft_ptr->ref_counter--;

			assert(fft_ptr->ref_counter >= 0);
			if (fft_ptr->ref_counter == 0)
			{
				delete *fftGetSingletonPtr();
				*fftGetSingletonPtr() = nullptr;

				{
					/*
					 * Handle anti-aliasing stuff here.
					 * We don't free the Anti-aliasing FFTW plans if their reference counter is zero
					 * since these this reference counter could get 0 during each time step.
					 * This frequent plan creation would significantly slow down the program.
					 */

					// also free the aliasing stuff
					FFTWSingletonClass* fft_anti_ptr = *fftAliasingGetSingletonPtr();
					if (fft_anti_ptr != nullptr)
					{
						if (fft_anti_ptr->ref_counter != 0)
						{
							std::cerr << "SWEET requires all PlaneDatas which were used for Antialiasing to be freed before the standard PlaneDatas!" << std::endl;
							assert(false);
							exit(1);
						}

						delete fft_anti_ptr;
						*fftAliasingGetSingletonPtr() = nullptr;
					}
				}
			}
		}


#else
		MemBlockAlloc::free(kernel_data, sizeof(double)*kernel_size*kernel_size);
#endif
	}


	inline
	void set(
			std::size_t j,
			std::size_t i,
			double i_value
	)
	{
		assert(i >= range_start[0] && i < range_end[0]);
		assert(j >= range_start[1] && j < range_end[1]);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;
#endif

		array_data_cartesian_space[
							(j-range_start[1])*range_size[0]+
							(i-range_start[0])
						] = i_value;


		checkConsistency();
	}



#if SWEET_USE_PLANE_SPECTRAL_SPACE
	inline
	void set_spec(
			std::size_t j,
			std::size_t i,
			std::complex<double> &i_value
	)
	{
		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

//		requestDataInSpectralSpace();

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

		std::size_t idx = ((j-range_spec_start[1])*range_spec_size[0]+(i-range_spec_start[0]))*2;
		array_data_spectral_space[idx] = i_value.real();
		array_data_spectral_space[idx+1] = i_value.imag();

		checkConsistency();
	}
#endif



	inline
	double get(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInCartesianSpace();

		assert(i >= range_start[0] && i < range_end[0]);
		assert(j >= range_start[1] && j < range_end[1]);

		return array_data_cartesian_space[
							(j-range_start[1])*range_size[0]+
							(i-range_start[0])
						];

		checkConsistency();
	}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	inline
	std::complex<double> get_spec(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInSpectralSpace();

		((PlaneData*)this)->array_data_spectral_space_valid = true;

		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		std::size_t idx = ((j-range_spec_start[1])*range_spec_size[0]+(i-range_spec_start[0]))*2;
		return std::complex<double>(array_data_spectral_space[idx], array_data_spectral_space[idx+1]);

		checkConsistency();
	}
#endif


	inline
	void set_all(
			double i_value
	)
	{
		// TODO: implement this in spectral space!

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] = i_value;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;
#endif

		checkConsistency();
	}


	/**
	 * Set the values in the specified row
	 */
	inline
	void set_row(
			int i_row,
			double i_value
	)
	{
		if (i_row < 0)
			i_row += resolution[1];

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]; i++)
			array_data_cartesian_space[i_row*resolution[0]+i] = i_value;

		checkConsistency();
	}

	/**
	 * Copy values from another row and flip sign
	 */
	inline
	void copy_row_inv_sign(
			int i_src_row,
			int i_dst_row
	)
	{
		if (i_src_row < 0)
			i_src_row += resolution[1];

		if (i_dst_row < 0)
			i_dst_row += resolution[1];

		requestDataInCartesianSpace();

		std::size_t src_idx = i_src_row*resolution[0];
		std::size_t dst_idx = i_dst_row*resolution[0];

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]; i++)
			array_data_cartesian_space[dst_idx+i] = -array_data_cartesian_space[src_idx+i];

		checkConsistency();
	}


	/**
	 * Copy values from another row
	 */
	inline
	void copy_row(
			int i_src_row,
			int i_dst_row
	)
	{
		if (i_src_row < 0)
			i_src_row = resolution[1]+i_src_row;

		if (i_dst_row < 0)
			i_dst_row = resolution[1]+i_dst_row;

		requestDataInCartesianSpace();

		std::size_t src_idx = i_src_row*resolution[0];
		std::size_t dst_idx = i_dst_row*resolution[0];

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]; i++)
			array_data_cartesian_space[dst_idx+i] = array_data_cartesian_space[src_idx+i];

		checkConsistency();
	}






#if SWEET_USE_PLANE_SPECTRAL_SPACE==1

	inline
	double spec_getRe(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInSpectralSpace();

		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		std::size_t idx =	(j-range_spec_start[1])*range_spec_size[0]+
							(i-range_spec_start[0]);

		checkConsistency();
		return array_data_spectral_space[idx*2+0];
	}



	inline
	double spec_getIm(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInSpectralSpace();

		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		std::size_t idx =	(j-range_spec_start[1])*range_spec_size[0]+
							(i-range_spec_start[0]);

		checkConsistency();
		return array_data_spectral_space[idx*2+1];
	}


	inline
	void set_spec(
			std::size_t j,
			std::size_t i,
			double i_value_re,
			double i_value_im
	)
	{
		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		std::size_t idx =	(j-range_spec_start[1])*range_spec_size[0]+
							(i-range_spec_start[0]);

		array_data_spectral_space[idx*2+0] = i_value_re;
		array_data_spectral_space[idx*2+1] = i_value_im;

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

		checkConsistency();
	}



	inline
	void set_spec_all(
			double i_value_re,
			double i_value_im
	)
	{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			array_data_spectral_space[i+0] = i_value_re;
			array_data_spectral_space[i+1] = i_value_im;
		}

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

		checkConsistency();
	}



	/**
	 * Set the spectrum of a frequency.
	 *
	 * This is different to the default spec_set function, since
	 * it directly sets up the y-mirrored frequencies as well.
	 *
	 * Note, that the x-frequencies are already mirrored.
	 */
	inline
	void set_spec_diff(
			std::size_t j,		///< j in [0, res[1]/2-1]
			std::size_t i,		///< i in [0, res[0]/2-1]
			double i_value_re,
			double i_value_im
	)
	{
		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		/*
		 * Note the padding in the x-direction:
		 *
		 * res_spec_x = 2 * (nx/2 + 1)
		 *
		 * SWEET (so far) only supports even resolution.
		 * Hence we have a padding of 1 complex value in the end.
		 *
		 * The frequencies in the x direction are then arranged in the following way:
		 *     [0, 1, 2, 3, ..., N/2-1, "padding=0"]
		 *
		 * Note, that setting a frequency here also sets the frequency on the "virtually mirrored" side.
		 *
		 * For the y-axis, the frequencies are given as follows (vector is transposed):
		 *
		 *     [0, 1, 2, 3, ..., N/2-1, N/2, N/2-1..., 3, 2, 1]
		 *
		 * Setting the frequency for N/2 to zero is a kind of obvious
		 */
		assert(i >= 0 && i < resolution[0]/2);
		assert(j >= 0 && j < resolution[1]/2);

		{	// lower part of y
			std::size_t idx =	(j-range_spec_start[1])*range_spec_size[0]+
								(i-range_spec_start[0]);

			array_data_spectral_space[idx*2+0] = i_value_re;
			array_data_spectral_space[idx*2+1] = i_value_im;
		}

		if (j != 0)
		{	// upper part of y
			std::size_t idx =	((resolution[1]-j)-range_spec_start[1])*range_spec_size[0]+
								(i-range_spec_start[0]);

			array_data_spectral_space[idx*2+0] = i_value_re;
			// IMPORTANT! the imaginary component is mirrored!!!
			array_data_spectral_space[idx*2+1] = i_value_im;
		}

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

		checkConsistency();
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void set_spec_spectrum_with_ampl_and_phase(
			std::size_t j,		///< j in [0, res[1]/2-1]
			std::size_t i,		///< i in [0, res[0]/2-1]
			double i_value_amplitude,	///< amplitude in |R
			double i_value_phase_shift	///< phase shift in [0;1[
	)
	{
		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		double c = cos(2.0*M_PIl*i_value_phase_shift)*i_value_amplitude;
		double s = sin(2.0*M_PIl*i_value_phase_shift)*i_value_amplitude;

		set_spec_spectrum(j, i, c, s);

		checkConsistency();
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void set_spec_spectrum(
			std::size_t j,		///< j in [0, res[1]/2-1]
			std::size_t i,		///< i in [0, res[0]/2-1]
			double i_value_re,	///< amplitude in |R
			double i_value_im	///< phase shift in [0;1[
	)
	{
		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		/*
		 * Note the padding in the x-direction:
		 *
		 * res_spec_x = 2 * (nx/2 + 1)
		 *
		 * SWEET (so far) only supports even resolution.
		 * Hence we have a padding of 1 complex value in the end.
		 *
		 * The frequencies in the x direction are then arranged in the following way:
		 *     [0, 1, 2, 3, ..., N/2-1, "padding=0"]
		 *
		 * Note, that setting a frequency here also sets the frequency on the "virtually mirrored" side.
		 *
		 * For the y-axis, the frequencies are given as follows (vector is transposed):
		 *
		 *     [0, 1, 2, 3, ..., N/2-1, N/2, N/2-1..., 3, 2, 1]
		 *
		 * Setting the frequency for N/2 to zero is a kind of obvious
		 */
		assert(i >= 0 && i < resolution[0]/2);
		assert(j >= 0 && j < resolution[1]/2);

		{	// lower part of y
			std::size_t idx =	(j-range_spec_start[1])*range_spec_size[0]+
								(i-range_spec_start[0]);

			array_data_spectral_space[idx*2+0] = i_value_re;
			array_data_spectral_space[idx*2+1] = i_value_im;
		}

		if (j != 0)
		{	// upper part of y
			std::size_t idx =	((resolution[1]-j)-range_spec_start[1])*range_spec_size[0]+
								(i-range_spec_start[0]);

			array_data_spectral_space[idx*2+0] = i_value_re;
			// IMPORTANT! the imaginary component is mirrored!!!
			array_data_spectral_space[idx*2+1] = -i_value_im;
		}

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

		checkConsistency();
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void set_spec_spectrum_A(
			std::size_t j,		///< j in [0, res[1]/2-1]
			std::size_t i,		///< i in [0, res[0]/2-1]
			double i_value_re,	///< amplitude in |R
			double i_value_im	///< phase shift in [0;1[
	)
	{
		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		/*
		 * Note the padding in the x-direction:
		 *
		 * res_spec_x = 2 * (nx/2 + 1)
		 *
		 * SWEET (so far) only supports even resolution.
		 * Hence we have a padding of 1 complex value in the end.
		 *
		 * The frequencies in the x direction are then arranged in the following way:
		 *     [0, 1, 2, 3, ..., N/2-1, "padding=0"]
		 *
		 * Note, that setting a frequency here also sets the frequency on the "virtually mirrored" side.
		 *
		 * For the y-axis, the frequencies are given as follows (vector is transposed):
		 *
		 *     [0, 1, 2, 3, ..., N/2-1, N/2, N/2-1..., 3, 2, 1]
		 *
		 * Setting the frequency for N/2 to zero is a kind of obvious
		 */
		assert(i >= 0 && i < resolution[0]/2);
		assert(j >= 0 && j < resolution[1]/2);

		{	// lower part of y
			std::size_t idx =	(j-range_spec_start[1])*range_spec_size[0]+
								(i-range_spec_start[0]);

			array_data_spectral_space[idx*2+0] = i_value_re;
			array_data_spectral_space[idx*2+1] = i_value_im;
		}

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

		checkConsistency();
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void set_spec_spectrum_B(
			std::size_t j,		///< j in [0, res[1]/2-1]
			std::size_t i,		///< i in [0, res[0]/2-1]
			double i_value_re,	///< amplitude in |R
			double i_value_im	///< phase shift in [0;1[
	)
	{
		assert(i >= range_spec_start[0] && i < range_spec_end[0]);
		assert(j >= range_spec_start[1] && j < range_spec_end[1]);

		/*
		 * Note the padding in the x-direction:
		 *
		 * res_spec_x = 2 * (nx/2 + 1)
		 *
		 * SWEET (so far) only supports even resolution.
		 * Hence we have a padding of 1 complex value in the end.
		 *
		 * The frequencies in the x direction are then arranged in the following way:
		 *     [0, 1, 2, 3, ..., N/2-1, "padding=0"]
		 *
		 * Note, that setting a frequency here also sets the frequency on the "virtually mirrored" side.
		 *
		 * For the y-axis, the frequencies are given as follows (vector is transposed):
		 *
		 *     [0, 1, 2, 3, ..., N/2-1, N/2, N/2-1..., 3, 2, 1]
		 *
		 * Setting the frequency for N/2 to zero is a kind of obvious
		 */
		assert(i >= 0 && i < resolution[0]/2);
		assert(j >= 0 && j < resolution[1]/2);

		if (j != 0)
		{	// upper part of y
			std::size_t idx =	((resolution[1]-j)-range_spec_start[1])*range_spec_size[0]+
								(i-range_spec_start[0]);

			array_data_spectral_space[idx*2+0] = i_value_re;
			// IMPORTANT! the imaginary component is mirrored!!!
			array_data_spectral_space[idx*2+1] = -i_value_im;
		}

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

		checkConsistency();
	}




private:
	static FFTWSingletonClass** fftGetSingletonPtr()
	{
		static FFTWSingletonClass *fftw_singleton_data = nullptr;
		return &fftw_singleton_data;
	}


private:
	static FFTWSingletonClass* fftTestAndInit(
		PlaneData &i_dataArray
	)
	{
		FFTWSingletonClass **fftw_singleton_data = fftGetSingletonPtr();

		if (*fftw_singleton_data != nullptr)
		{
			(*fftw_singleton_data)->ref_counter++;
			return *fftw_singleton_data;
		}


		*fftw_singleton_data = new FFTWSingletonClass(i_dataArray);
		(*fftw_singleton_data)->ref_counter++;

		return *fftw_singleton_data;
	}

private:
	/**
	 * TODO: Replace this with returning a reference
	 */
	static FFTWSingletonClass** fftAliasingGetSingletonPtr()
	{
		static FFTWSingletonClass *fftw_singleton_data = nullptr;
		return &fftw_singleton_data;
	}

private:
	FFTWSingletonClass* fftAliasingTestAndInit(
		PlaneData &i_dataArray
	)	const
	{
		assert(i_dataArray.aliasing_scaled);

		FFTWSingletonClass **fftw_singleton_data = fftAliasingGetSingletonPtr();

		if (*fftw_singleton_data != nullptr)
		{
			(*fftw_singleton_data)->ref_counter++;
			return *fftw_singleton_data;
		}

		*fftw_singleton_data = new FFTWSingletonClass(i_dataArray);
		(*fftw_singleton_data)->ref_counter++;

		return *fftw_singleton_data;
	}

public:
	static void checkRefCounters()
	{
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING || SWEET_USE_PLANE_SPECTRAL_SPACE
		if (*PlaneData::fftGetSingletonPtr() != nullptr)
		{
			std::cerr << "FFT plans not yet released." << std::endl;
			std::cerr << "This typically means a memory leak if this is the end of the program and all classes deconstructors have been called" << std::endl;
			std::cerr << " Reference counter: " << (*PlaneData::fftGetSingletonPtr())->ref_counter << std::endl;
		}

		if (*PlaneData::fftAliasingGetSingletonPtr() != nullptr)
		{
			std::cerr << "FFT plans for ANTI ALIASING not yet released." << std::endl;
			std::cerr << "This typically means a memory leak if this is the end of the program and all classes deconstructors have been called" << std::endl;
			std::cerr << " Reference counter: " << (*PlaneData::fftAliasingGetSingletonPtr())->ref_counter << std::endl;
		}
#endif
	}

#else


public:
	static void checkRefCounters()
	{
	}

#endif

public:
	inline
	const PlaneData& requestDataInSpectralSpace() const
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE==0
		std::cerr << "requestDataInSpectralSpace: spectral space is disabled" << std::endl;
		exit(-1);
#else


#if SWEET_DEBUG
	#if SWEET_THREADING
		if (omp_get_num_threads() > 0)
		{
			std::cerr << "Threading race conditions likely" << std::endl;
			assert(false);
		}
	#endif
#endif

		if (array_data_spectral_space_valid)
			return *this;		// nothing to do

		if (!array_data_cartesian_space_valid)
		{
			std::cerr << "Spectral data not available! Is this maybe a non-initialized operator?" << std::endl;
			assert(false);
			exit(1);
		}

		PlaneData *rw_array_data = (PlaneData*)this;

		if (aliasing_scaled)
			(*fftAliasingGetSingletonPtr())->fft_forward(*rw_array_data);
		else
			(*fftGetSingletonPtr())->fft_forward(*rw_array_data);

		rw_array_data->array_data_spectral_space_valid = true;

		/*
		 * TODO: this is only a debugging helper
		 */
//		rw_array_data->array_data_cartesian_space_valid = false;
		checkConsistency();
#endif
		return *this;
	}


	inline
	const PlaneData& requestDataInCartesianSpace() const
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE==1
		checkConsistency();

		if (array_data_cartesian_space_valid)
			return *this;		// nothing to do

		assert(array_data_spectral_space_valid == true);

		PlaneData *rw_array_data = (PlaneData*)this;

		if (!aliasing_scaled)
			(*fftGetSingletonPtr())->fft_backward(*rw_array_data);
		else
			(*fftAliasingGetSingletonPtr())->fft_backward(*rw_array_data);

		rw_array_data->array_data_cartesian_space_valid = true;

		/*
		 * TODO: this is only a debugging helper
		 */
//		rw_array_data->array_data_spectral_space_valid = false;

		checkConsistency();
#endif
	return *this;
	}


	inline
	PlaneData return_one_if_positive()
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] > 0 ? 1 : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif

		out.checkConsistency();
		return out;
	}



	inline
	PlaneData return_value_if_positive()	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] > 0 ? array_data_cartesian_space[i] : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif

		out.checkConsistency();
		return out;
	}


	inline
	PlaneData return_one_if_negative()	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] < 0 ? 1 : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif

		out.checkConsistency();
		return out;
	}


	inline
	PlaneData return_value_if_negative()	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] < 0 ? array_data_cartesian_space[i] : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif

		out.checkConsistency();
		return out;
	}



	/**
	 * return true, if any value is infinity
	 */
	bool reduce_all_finite() const
	{
		requestDataInCartesianSpace();

		bool isallfinite = true;
#if SWEET_THREADING
#pragma omp parallel for reduction(&&:isallfinite)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			isallfinite = isallfinite && std::isfinite(array_data_cartesian_space[i]);
//			isallfinite = isallfinite && (array_data_cartesian_space[i]<1000);


		return isallfinite;
	}



	/**
	 * return the maximum of all absolute values
	 */
	double reduce_maxAbs()	const
	{
		requestDataInCartesianSpace();

		double maxabs = -1;
#if SWEET_THREADING
#pragma omp parallel for reduction(max:maxabs)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			maxabs = std::max(maxabs, std::abs(array_data_cartesian_space[i]));


		checkConsistency();
		return maxabs;
	}



	/**
	 * reduce to root mean square
	 */
	double reduce_rms()
	{
		requestDataInCartesianSpace();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += array_data_cartesian_space[i]*array_data_cartesian_space[i];

		sum = std::sqrt(sum/double(array_data_cartesian_length));

		checkConsistency();
		return sum;
	}


	/**
	 * reduce to root mean square
	 */
	double reduce_rms_quad()
	{
		requestDataInCartesianSpace();

		double sum = 0;
		double c = 0;

#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
		{
			double value = array_data_cartesian_space[i]*array_data_cartesian_space[i];

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		sum = std::sqrt(sum/double(array_data_cartesian_length));

		checkConsistency();
		return sum;
	}



	/**
	 * return the maximum of all absolute values
	 */
	double reduce_max()	const
	{
		requestDataInCartesianSpace();

		double maxvalue = -std::numeric_limits<double>::max();
#if SWEET_THREADING
#pragma omp parallel for reduction(max:maxvalue)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			maxvalue = std::max(maxvalue, array_data_cartesian_space[i]);


		checkConsistency();
		return maxvalue;
	}


	/**
	 * return the maximum of all absolute values
	 */
	double reduce_min()	const
	{
		requestDataInCartesianSpace();

		double minvalue = std::numeric_limits<double>::max();
#if SWEET_THREADING
#pragma omp parallel for reduction(min:minvalue)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			minvalue = std::min(minvalue, array_data_cartesian_space[i]);


		checkConsistency();
		return minvalue;
	}


	/**
	 * return the maximum of all absolute values
	 */
	double reduce_sum()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += array_data_cartesian_space[i];


		checkConsistency();
		return sum;
	}


	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double reduce_sum_quad()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
		double c = 0;
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
		{
			double value = array_data_cartesian_space[i];

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;


		checkConsistency();
		return sum;
	}

	/**
	 * return the maximum of all absolute values
	 */
	double reduce_norm1()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += std::abs(array_data_cartesian_space[i]);


		checkConsistency();
		return sum;
	}

	/**
	 * return the sum of the absolute values.
	 */
	double reduce_norm1_quad()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
		double c = 0;
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
		{

			double value = std::abs(array_data_cartesian_space[i]);
			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;


		checkConsistency();
		return sum;
	}


	/**
	 * return the sqrt of the sum of the squared values
	 */
	double reduce_norm2()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += array_data_cartesian_space[i]*array_data_cartesian_space[i];


		checkConsistency();
		return std::sqrt(sum);
	}


	/**
	 * return the sqrt of the sum of the squared values, use quad precision for reduction
	 */
	double reduce_norm2_quad()	const
	{
		requestDataInCartesianSpace();

		double sum = 0.0;
		double c = 0.0;

#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
		{
			double value = array_data_cartesian_space[i]*array_data_cartesian_space[i];

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;


		checkConsistency();
		return std::sqrt(sum);
	}


#if SWEET_USE_PLANE_SPECTRAL_SPACE

	/**
	 * return the maximum of all absolute values
	 */
	double reduce_spec_maxAbs()	const
	{
		requestDataInSpectralSpace();

		double maxabs = -1;
#if SWEET_THREADING
#pragma omp parallel for reduction(max:maxabs)
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			double re = array_data_spectral_space[i];
			double im = array_data_spectral_space[i+1];
			maxabs = std::max(maxabs, std::sqrt(re*re + im*im));
		}


		checkConsistency();
		return maxabs;
	}


	/**
	 * return centroid of frequency:
	 *
	 * Note, that the centroid is given in logarithmic space
	 */
	double reduce_spec_getPolvaniCentroid()	const
	{
		requestDataInSpectralSpace();

		double nom = 0;
		double denom = 0;

		for (std::size_t j = 0; j < resolution_spec[1]; j++)
		{
			for (std::size_t i = 0; i < resolution_spec[0]; i++)
			{
				double re = spec_getRe(j, i);
				double im = spec_getIm(j, i);

				std::size_t ka = i;
//				std::size_t ka = (i < resolution_spec[0] ? i : resolution_spec[0]-i);

				// consider mirroring in y direction
				std::size_t kb = (j < resolution_spec[1]/2 ? j : resolution_spec[1]-j);

				double k = std::sqrt((double)(ka*ka)+(double)(kb*kb));
				double value = std::sqrt(re*re+im*im);

				nom += k*value;
				denom += value;
			}
		}


		checkConsistency();
		return nom/denom;
	}
#endif

	constexpr
	static
	int get_kernel_mask3x3(
			int i_0,
			int i_1,
			int i_2,
			int i_3,
			int i_4,
			int i_5,
			int i_6,
			int i_7,
			int i_8
	)
	{
		return
				(i_0 << 0) |
				(i_1 << 1) |
				(i_2 << 2) |
				(i_3 << 3) |
				(i_4 << 4) |
				(i_5 << 5) |
				(i_6 << 6) |
				(i_7 << 7) |
				(i_8 << 8);
	}


public:
	template <int S>
	void kernel_stencil_setup(
			const double i_kernel_array[S][S],
			double i_scale = 1.0
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE == 0

		kernel_size = S;
		kernel_data = MemBlockAlloc::alloc<double>(sizeof(double)*S*S);
		for (int y = 0; y < S; y++)
			for (int x = 0; x < S; x++)
				kernel_data[y*S+x] = i_kernel_array[S-1-y][x];

		for (int i = 0; i < S*S; i++)
			kernel_data[i] *= i_scale;


		if (S == 3)
		{
			kernel_id = get_kernel_mask3x3(
					kernel_data[0] != 0,
					kernel_data[1] != 0,
					kernel_data[2] != 0,
					kernel_data[3] != 0,
					kernel_data[4] != 0,
					kernel_data[5] != 0,
					kernel_data[6] != 0,
					kernel_data[7] != 0,
					kernel_data[8] != 0
				);
		}
		else
		{
			kernel_id = -1;
		}

#else

		double inv_kernel_array[S][S];

		for (int j = 0; j < S; j++)
			for (int i = 0; i < S; i++)
				inv_kernel_array[j][i] = i_kernel_array[j][S-i-1]*i_scale;

		// assure symmetric kernel
		assert((S & 1) == 1);

		// radius of kernel (half size)
		std::size_t R = S>>1;

		set_all(0);

		// left lower corner
		//kernel_cart[    0:R[0]+1,       0:R[0]+1    ] = conv_kernel[    R[0]:,  R[0]:   ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = ky-R;
				std::size_t cx = kx-R;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// right bottom corner
		//kernel_cart[    0:R[0]+1,       res-R[0]:   ] = conv_kernel[    R[0]:,  0:R[0]  ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = ky-R;
				std::size_t cx = resolution[0] - R + kx;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// left top corner
		//kernel_cart[    res-R[0]:,      0:R[0]+1   ] = conv_kernel[    0:R[0],  R[0]:  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = resolution[1] - R + ky;
				std::size_t cx = kx-R;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// right top corner
		//kernel_cart[    res-R[0]:,      res-R[0]:   ] = conv_kernel[    0:R[0], 0:R[0]  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = resolution[1] - R + ky;
				std::size_t cx = resolution[0] - R + kx;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				set(cy, cx, inv_kernel_array[ky][kx]);
			}

		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;

		requestDataInSpectralSpace();	/// convert kernel_data; to spectral space
#endif

		checkConsistency();
	}


	/**
	 * apply a 3x3 stencil
	 */
	PlaneData op_stencil_Re_3x3(
			const double *i_kernel_data
	)
	{
		PlaneData out = *this;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		bool was_spectral = false;
		if (this->array_data_spectral_space_valid)
		{
			this->requestDataInCartesianSpace();
			was_spectral = true;
			this->array_data_spectral_space_valid=false;
		}else if(!this->array_data_cartesian_space_valid){
			std::cout << "Uninitialized PlaneData in op_stencil_Re_3x3" << std::endl;
			exit(-1);
		}
		out.array_data_spectral_space_valid=false;
#endif


		int res_x = resolution[0];
		int res_y = resolution[1];


#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y, out, i_kernel_data)
#endif
		for (int y = 0; y < res_y; y++)
		{
			for (int x = 0; x < res_x; x++)
			{
				double *data_out = &out.array_data_cartesian_space[(y*res_x+x)];
				data_out[0] = 0;

				for (int j = -1; j <= 1; j++)
				{
					int pos_y = y+j;

					pos_y -= (pos_y >= res_y ? res_y : 0);
					pos_y += (pos_y < 0 ? res_y : 0);

					assert(pos_y >= 0 && pos_y < res_y);

					for (int i = -1; i <= 1; i++)
					{
						int pos_x = x+i;

						pos_x -= (pos_x >= res_x ? res_x : 0);
						pos_x += (pos_x < 0 ? res_x : 0);

						assert(pos_x >= 0 && pos_x < res_x);

						int idx = (j+1)*3+(i+1);
						assert(idx >= 0 && idx < 9);

						double kre = i_kernel_data[idx];

						double dre = array_data_cartesian_space[(pos_y*res_x+pos_x)];

						data_out[0] += dre*kre;
					}
				}
			}
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (was_spectral)
		{
			this->requestDataInSpectralSpace();
			out.requestDataInSpectralSpace();
		}
#endif

		return out;
	}


	/**
	 * apply a 3x3 stencil with updates
	 */
	PlaneData op_stencil_Re_3x3_update(
			const double *i_kernel_data
	)
	{
		PlaneData out(resolution);
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		bool was_spectral = false;
		if (this->array_data_spectral_space_valid)
		{
			this->requestDataInCartesianSpace();
			was_spectral = true;
			this->array_data_spectral_space_valid=false;
		}else if(!this->array_data_cartesian_space_valid){
			std::cout << "Uninitialized PlaneData in op_stencil_Re_3x3" << std::endl;
			exit(-1);
		}
		out.array_data_spectral_space_valid=false;
#endif

		int res_x = resolution[0];
		int res_y = resolution[1];


#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y, out, i_kernel_data)
#endif
		for (int y = 0; y < res_y; y++)
		{
			for (int x = 0; x < res_x; x++)
			{
				double *data_out = &out.array_data_cartesian_space[(y*res_x+x)];
				data_out[0] = 0;

				for (int j = -1; j <= 1; j++)
				{
					int pos_y = y+j;

					pos_y -= (pos_y >= res_y ? res_y : 0);
					pos_y += (pos_y < 0 ? res_y : 0);

					assert(pos_y >= 0 && pos_y < res_y);

					for (int i = -1; i <= 1; i++)
					{
						int pos_x = x+i;

						pos_x -= (pos_x >= res_x ? res_x : 0);
						pos_x += (pos_x < 0 ? res_x : 0);

						assert(pos_x >= 0 && pos_x < res_x);

						int idx = (j+1)*3+(i+1);
						assert(idx >= 0 && idx < 9);

						double kre = i_kernel_data[idx];

						double dre = array_data_cartesian_space[(pos_y*res_x+pos_x)];

						data_out[0] += dre*kre;
					}
				}
			}
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (was_spectral)
		{
			this->requestDataInSpectralSpace();
			out.requestDataInSpectralSpace();
		}
#endif

		return out;
	}


#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Invert the application of a linear operator in spectral space.
	 * The operator is given in i_array_data
	 */
	inline
	PlaneData spec_div_element_wise(
			const PlaneData &i_array_data,	///< operator
			double i_denom_zeros_scalar = 0.0,
			double i_tolerance = 0.0			///< set tolerance to 0, since we setup the values in the spectral operator directly
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		// only makes sense, if this is an operator created in spectral space
		assert(i_array_data.array_data_spectral_space_valid == true);

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

		// determine maximum value for tolerance
		double max_value = i_array_data.reduce_spec_maxAbs();
		i_tolerance *= max_value;
		i_tolerance *= (resolution[0]+resolution[1]);	// the larger the matrix, the less the accuracy

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			double ar = array_data_spectral_space[i];
			double ai = array_data_spectral_space[i+1];
			double br = i_array_data.array_data_spectral_space[i];
			double bi = i_array_data.array_data_spectral_space[i+1];

			double den = (br*br+bi*bi);

			if (std::abs(den) <= i_tolerance)
			{
				// For inverting differential operators, this is the integration constant C
				out.array_data_spectral_space[i] = ar*i_denom_zeros_scalar;
				out.array_data_spectral_space[i+1] = ai*i_denom_zeros_scalar;
			}
			else
			{
				out.array_data_spectral_space[i] = (ar*br + ai*bi)/den;
				out.array_data_spectral_space[i+1] = (ai*br - ar*bi)/den;
			}
		}

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;


		checkConsistency();
		return out;
	}
#endif


#if SWEET_USE_PLANE_SPECTRAL_DEALIASING && SWEET_USE_PLANE_SPECTRAL_SPACE

	/**
	 * scale the solution to avoid aliasing effects for handling non-linear terms
	 */
	inline
	PlaneData aliasing_scaleUp(
			std::size_t *i_new_resolution = nullptr
	)	const
	{
		std::size_t new_resolution[2];

		//Resolution for augmented cartesian space arrays (3N/2)
		if (i_new_resolution == nullptr)
		{
			for (int i = 0; i < 2; i++)
			{
				// TODO
//				new_resolution[i] = (this->resolution[i]*4)/2;
				new_resolution[i] = (this->resolution[i]*3)/2;
//				new_resolution[i] = (this->resolution[i]*3)/2-1;
//				new_resolution[i] = (this->resolution[i]*2);

				if ((new_resolution[i] & 1) != 0)
				{
					std::cerr << "Odd resolution for dealiasing not supported, please change your resolution" << std::endl;
					exit(-1);
				}
			}
		}
		else
		{
			for (int i = 0; i < 2; i++)
				new_resolution[i] = i_new_resolution[i];
		}

		// Augmented array (temporary)
		PlaneData out(new_resolution);
//		out.temporary_data = false;
		out.aliasing_scaled = true;

		/*
		 * Test if the FFTW plans for aliasing were initialized before.
		 * If not, then initialize them.
		 * These FFTW plans are of higher resolution than the standard ones.
		 */
		fftAliasingTestAndInit(out);

		//Get spectral data
		requestDataInSpectralSpace();

		//Zero temporary array
		out.set_spec_all(0, 0);

		// TODO: this does not work once distributed memory is available
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution_spec[1]/2; j++)
		{
			/*
			 * Copy the spectrum of data_array to the temporary array
			 *    --> this will leave blank the high modes of the tmp array, since it has 3N/2 modes instead of N
			 * We use memcpy, since there are typically optimized routines available.
			 * Writing this as a for loop can result in similar or even better performance.
			 * lower quadrant
			 */
			memcpy(
					out.array_data_spectral_space+(j*out.range_spec_size[0])*2,
					array_data_spectral_space+(j*range_spec_size[0])*2,
					sizeof(double)*(range_spec_size[0]*2-2)
				);
		}

		std::size_t disp = out.range_spec_size[1] - range_spec_size[1];
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = resolution_spec[1]/2+1; j < resolution_spec[1]; j++)
		{
			// top quadrant
			memcpy(
					out.array_data_spectral_space+((j+disp)*out.range_spec_size[0])*2,
					array_data_spectral_space+(j*range_spec_size[0])*2,
					sizeof(double)*(range_spec_size[0]*2-2)
				);
		}

		//Scale for temporary array for the correct Fourier constant relative to resolution
		double scale = ((double)new_resolution[0]*(double)new_resolution[1])/((double)resolution[0]*(double)resolution[1]);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < out.array_data_spectral_length; i++)
			out.array_data_spectral_space[i] *= scale;

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;


		checkConsistency();

		/*
		 *The temporary array becomes the main array, but this means that it will have 3N/2 modes and ...
		 * The Cartesian space data is lost.
		 * In general, both data sets have to match when converted forward or backward or only one of them is valid.
		 * If they are not consistently matching (converting forward/backward), this would result in an exception.
		 */

		return out;
	}

	PlaneData aliasing_scaleDown(
			std::size_t *i_new_resolution
	)
	{
		assert(aliasing_scaled == true);
//		aliasing_scaled = true;

		PlaneData out(i_new_resolution);
//		out.temporary_data = false;
		out.aliasing_scaled = false;

//		fftAliasingTestAndInit(out);

		requestDataInSpectralSpace();

		out.set_spec_all(0, 0);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < out.resolution_spec[1]/2; j++)
		{
			// lower quadrant
			memcpy(
					out.array_data_spectral_space+(j*out.range_spec_size[0])*2,
					array_data_spectral_space+(j*range_spec_size[0])*2,
					sizeof(double)*(out.range_spec_size[0]*2-2)
				);
		}

		std::size_t disp = range_spec_size[1] - out.range_spec_size[1];
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = out.resolution_spec[1]/2+1; j < out.resolution_spec[1]; j++)
		{
			// top quadrant
			memcpy(
					out.array_data_spectral_space+(j*out.range_spec_size[0])*2,
					array_data_spectral_space+((j+disp)*range_spec_size[0])*2,
					sizeof(double)*(out.range_spec_size[0]*2-2)
				);
		}
		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;


		double scale = ((double)i_new_resolution[0]*(double)i_new_resolution[1])/((double)resolution[0]*(double)resolution[1]);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < out.array_data_spectral_length; i++)
			out.array_data_spectral_space[i] *= scale;



		checkConsistency();
		return out;
	}
#endif


#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Zero high frequency modes (beyond 2N/3)
	 *
	 *Example of Spectrum with N=16: all high modes will be set to zero
	 * 	(0, 1, Low )	(1, 1, Low )	(2, 1, Low )	(3, 1, Low )	(4, 1, Low )	(5, 1, Low )	(6, 1, High)	(7, 1, High)	(8, 1, High)
	 * 	(0, 2, Low )	(1, 2, Low )	(2, 2, Low )	(3, 2, Low )	(4, 2, Low )	(5, 2, Low )	(6, 2, High)	(7, 2, High)	(8, 2, High)
	 * 	(0, 3, Low )	(1, 3, Low )	(2, 3, Low )	(3, 3, Low )	(4, 3, Low )	(5, 3, Low )	(6, 3, High)	(7, 3, High)	(8, 3, High)
	 *	(0, 4, Low )	(1, 4, Low )	(2, 4, Low )	(3, 4, Low )	(4, 4, Low )	(5, 4, Low )	(6, 4, High)	(7, 4, High)	(8, 4, High)
	 *	(0, 5, Low )	(1, 5, Low )	(2, 5, Low )	(3, 5, Low )	(4, 5, Low )	(5, 5, Low )	(6, 5, High)	(7, 5, High)	(8, 5, High)
	 *	(0, 6, High)	(1, 6, High)	(2, 6, High)	(3, 6, High)	(4, 6, High)	(5, 6, High)	(6, 6, High)	(7, 6, High)	(8, 6, High)
	 *	(0, 7, High)	(1, 7, High)	(2, 7, High)	(3, 7, High)	(4, 7, High)	(5, 7, High)	(6, 7, High)	(7, 7, High)	(8, 7, High)
	 *	(0, 8, High)	(1, 8, High)	(2, 8, High)	(3, 8, High)	(4, 8, High)	(5, 8, High)	(6, 8, High)	(7, 8, High)	(8, 8, High)
	 *	(0, 7, High)	(1, 7, High)	(2, 7, High)	(3, 7, High)	(4, 7, High)	(5, 7, High)	(6, 7, High)	(7, 7, High)	(8, 7, High)
	 *	(0, 6, High)	(1, 6, High)	(2, 6, High)	(3, 6, High)	(4, 6, High)	(5, 6, High)	(6, 6, High)	(7, 6, High)	(8, 6, High)
	 *	(0, 5, Low)		(1, 5, Low)		(2, 5, Low)		(3, 5, Low)		(4, 5, Low)		(5, 5, Low)		(6, 5, High)	(7, 5, High)	(8, 5, High)
	 *	(0, 4, Low)		(1, 4, Low)		(2, 4, Low)		(3, 4, Low)		(4, 4, Low)		(5, 4, Low)		(6, 4, High)	(7, 4, High)	(8, 4, High)
	 *	(0, 3, Low)		(1, 3, Low)		(2, 3, Low)		(3, 3, Low)		(4, 3, Low)		(5, 3, Low)		(6, 3, High)	(7, 3, High)	(8, 3, High)
	 *	(0, 2, Low)		(1, 2, Low)		(2, 2, Low)		(3, 2, Low)		(4, 2, Low)		(5, 2, Low)		(6, 2, High)	(7, 2, High)	(8, 2, High)
	 *	(0, 1, Low)		(1, 1, Low)		(2, 1, Low)		(3, 1, Low)		(4, 1, Low)		(5, 1, Low)		(6, 1, High)	(7, 1, High)	(8, 1, High)
	 *	(0, 0, Low)		(1, 0, Low)		(2, 0, Low)		(3, 0, Low)		(4, 0, Low)		(5, 0, Low)		(6, 0, High)	(7, 0, High)	(8, 0, High)
	 *
	 *
	 */
	inline
	PlaneData& aliasing_zero_high_modes()
	{
		//std::cout<<"Cartesian"<<std::endl;
		//printArrayData();
		//Get spectral data
		requestDataInSpectralSpace();
		//std::cout<<"Spectral input data"<<std::endl;
		//printSpectrum();

		//Upper part
		std::size_t i=1; //modenumber in x
		std::size_t j=1; //modenumber in y

		// TODO: this does not work once distributed memory is available
#if SWEET_THREADING
//#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t y = resolution_spec[1]-1; y > resolution_spec[1]/2; y--)
		{
			i=0;
			for (std::size_t x = 0; x < resolution_spec[0]; x++)
			{

				//double value_re = spec_getRe(y, x);
				//double value_im = spec_getIm(y, x);
				if( x > 2*(resolution_spec[0]-1)/3 || j > 2*(resolution_spec[1]/2)/3 )
				//if( x > 2*(resolution_spec[0]-1)/3-1 || j > 2*(resolution_spec[1]/2)/3-1 )
				{
					set_spec( y, x, 0.0, 0.0);
					//std::cout << "(" << i << ", " << j << ", High)\t";
				}
				else
				{
					//std::cout << "(" << i << ", " << j << ", Low )\t";
				}
				i++;
			}
			j++;
			//std::cout << std::endl;
		}


		//Lower part
		i=0; //modenumber in x
		j=resolution_spec[1]/2; //modenumber in y
		// TODO: this does not work once distributed memory is available
#if SWEET_THREADING
//#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (int y = (int) resolution_spec[1]/2; y >= 0; y--)
		{
			i=0;
			for (std::size_t x = 0; x < resolution_spec[0]; x++)
			{
				//double value_re = spec_getRe(y, x);
				//double value_im = spec_getIm(y, x);
				if( x > 2*(resolution_spec[0]-1)/3 ||  y > 2*((int) resolution_spec[1]/2)/3 )
				//if( x > 2*(resolution_spec[0]-1)/3-1 ||  y > 2*((int) resolution_spec[1]/2)/3-1 )
				{
					set_spec( y, x, 0.0, 0.0);
					//std::cout << "(" << i << ", " << j << ", High)\t";
				}
				else
				{
					//std::cout << "(" << i << ", " << j << ", Low)\t";
				}
				i++;
			}
			j--;
			//std::cout << std::endl;
		}

		requestDataInCartesianSpace();

		checkConsistency();

		return *this;
	}
#endif


public:
	/**
	 * assignment operator
	 */
	PlaneData &operator=(double i_value)
	{
		set_all(i_value);

		checkConsistency();
		return *this;
	}


public:
	/**
	 * assignment operator
	 */
	PlaneData &operator=(int i_value)
	{
		set_all(i_value);

		checkConsistency();
		return *this;
	}


public:
	/**
	 * assignment operator
	 *
	 * hasdfasdf = h*hasdf;
	 */
	PlaneData &operator=(
			const PlaneData &i_dataArray
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		aliasing_scaled = i_dataArray.aliasing_scaled;
#endif

		for (int i = 0; i < 2; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (i_dataArray.array_data_cartesian_space_valid)
		{
			array_data_cartesian_space_valid = true;
#endif
			/**
			 * If this data was generated based on temporary data sets (e.g. via h = hu/u), then only swap pointers.
			 */
//			if (i_dataArray.temporary_data)
//			{
//				std::swap(array_data_cartesian_space, ((PlaneData &)i_dataArray).array_data_cartesian_space);
//			}
//			else
			{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
				for (std::size_t i = 0; i < array_data_cartesian_length; i++)
					array_data_cartesian_space[i] = i_dataArray.array_data_cartesian_space[i];
			}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		}
		else
		{
			array_data_cartesian_space_valid = false;
		}

		array_data_spectral_length = i_dataArray.array_data_spectral_length;

		if (i_dataArray.array_data_spectral_space_valid)
		{
			array_data_spectral_space_valid = true;
//			if (i_dataArray.temporary_data)
//			{
//				std::swap(array_data_spectral_space, ((PlaneData &)i_dataArray).array_data_spectral_space);
//			}
//			else
			{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
				for (std::size_t i = 0; i < array_data_spectral_length; i++)
					array_data_spectral_space[i] = i_dataArray.array_data_spectral_space[i];
			}
		}
		else
		{
			array_data_spectral_space_valid = false;
		}
#endif

		checkConsistency();
		return *this;
	}


	/**
	 * Apply a linear operator given by this class to the input data array.
	 */
	inline
	PlaneData operator()(
			const PlaneData &i_array_data
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			double ar = array_data_spectral_space[i];
			double ai = array_data_spectral_space[i+1];
			double br = i_array_data.array_data_spectral_space[i];
			double bi = i_array_data.array_data_spectral_space[i+1];

			out.array_data_spectral_space[i] = ar*br - ai*bi;
			out.array_data_spectral_space[i+1] = ar*bi + ai*br;
		}

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;

#else

		/**
		 * TODO: optimize this!!!
		 *
		 *  - cache blocking
		 *  - if branching elimination
		 *  - etc.....
		 */
		int res_x = resolution[0];
		int res_y = resolution[1];

		if (kernel_size == 3)
		{
			switch (kernel_id)
			{
			case get_kernel_mask3x3(0, 0, 0, 1, 0, 1, 0, 0, 0):	// (X, 0, X)
#if SWEET_THREADING
#pragma omp parallel for schedule(static)
#endif
				for (int y = 0; y < res_y; y++)
				{
					for (int x = 0; x < res_x; x++)
					{
						double &data_out = out.array_data_cartesian_space[y*res_x+x];
						data_out = 0;

						int pos_y = y;
						assert(pos_y >= 0 && pos_y < res_y);

						if (x > 0 && x < res_x-1)
						{
							double *kernel_scalar_ptr = &kernel_data[3];
							double *data_scalar_ptr = &i_array_data.array_data_cartesian_space[pos_y*res_x+x-1];

							data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
							data_out += kernel_scalar_ptr[2]*data_scalar_ptr[2];
						}
						else
						{
							for (int i = -1; i <= 1; i+=2)
							{
								int pos_x = x+i;
								pos_x -= (pos_x >= res_x ? res_x : 0);
								pos_x += (pos_x < 0 ? res_x : 0);
								int idx = i+4;
								double kernel_scalar = kernel_data[idx];
								double data_scalar = i_array_data.array_data_cartesian_space[pos_y*res_x+pos_x];

								data_out += kernel_scalar*data_scalar;
							}
						}
					}
				}
				break;



			case get_kernel_mask3x3(0, 0, 0, 1, 1, 1, 0, 0, 0):	// (X, X, X)
#if SWEET_THREADING
#pragma omp parallel for schedule(static)
#endif
					for (int y = 0; y < res_y; y++)
					{
						for (int x = 0; x < res_x; x++)
						{
							double &data_out = out.array_data_cartesian_space[y*res_x+x];
							data_out = 0;

							int pos_y = y+res_y;
							pos_y -= (pos_y >= res_y ? res_y : 0);
							pos_y += (pos_y < 0 ? res_y : 0);

							assert(pos_y >= 0 && pos_y < res_y);

							if (x > 0 && x < res_x-1)
							{
								double *kernel_scalar_ptr = &kernel_data[3];
								double *data_scalar_ptr = &i_array_data.array_data_cartesian_space[pos_y*res_x+x-1];

								data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
								data_out += kernel_scalar_ptr[1]*data_scalar_ptr[1];
								data_out += kernel_scalar_ptr[2]*data_scalar_ptr[2];
							}
							else
							{
								for (int i = -1; i <= 1; i++)
								{
									int pos_x = (x+i);
									pos_x -= (pos_x >= res_x ? res_x : 0);
									pos_x += (pos_x < 0 ? res_x : 0);
									int idx = i+4;
									double kernel_scalar = kernel_data[idx];
									double data_scalar = i_array_data.array_data_cartesian_space[pos_y*res_x+pos_x];

									data_out += kernel_scalar*data_scalar;
								}
							}
						}
					}
					break;

			case get_kernel_mask3x3(0, 1, 0, 0, 0, 0, 0, 1, 0):	// (X, 0, X)^T
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
					for (int y = 0; y < res_y; y++)
					{
						for (int x = 0; x < res_x; x++)
						{
							double &data_out = out.array_data_cartesian_space[y*res_x+x];
							data_out = 0;

							if (y > 0 && y < res_y-1)
							{
								double *kernel_scalar_ptr = &kernel_data[1];
								double *data_scalar_ptr = &i_array_data.array_data_cartesian_space[(y-1)*res_x+x];

								data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
								data_out += kernel_scalar_ptr[6]*data_scalar_ptr[2*res_x];
							}
							else
							{
								for (int j = -1; j <= 1; j+=2)
								{
									int pos_y = y+j;

									pos_y -= (pos_y >= res_y ? res_y : 0);
									pos_y += (pos_y < 0 ? res_y : 0);

									int idx = (j+1)*3+1;

									double kernel_scalar = kernel_data[idx];
									double data_scalar = i_array_data.array_data_cartesian_space[pos_y*res_x+x];

									data_out += kernel_scalar*data_scalar;
								}
							}
						}
					}
					break;

			case get_kernel_mask3x3(0, 1, 0, 0, 1, 0, 0, 1, 0):	// (X, 0, X)^T
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
					for (int y = 0; y < res_y; y++)
					{
						for (int x = 0; x < res_x; x++)
						{
							double &data_out = out.array_data_cartesian_space[y*res_x+x];
							data_out = 0;

							if (y > 0 && y < res_y-1)
							{
								double *kernel_scalar_ptr = &kernel_data[1];
								double *data_scalar_ptr = &i_array_data.array_data_cartesian_space[(y-1)*res_x+x];

								data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
								data_out += kernel_scalar_ptr[3]*data_scalar_ptr[res_x];
								data_out += kernel_scalar_ptr[6]*data_scalar_ptr[2*res_x];
							}
							else
							{
								for (int j = -1; j <= 1; j++)
								{
									int pos_y = y+j;

									pos_y -= (pos_y >= res_y ? res_y : 0);
									pos_y += (pos_y < 0 ? res_y : 0);

									int idx = (j+1)*3+1;

									double kernel_scalar = kernel_data[idx];
									double data_scalar = i_array_data.array_data_cartesian_space[pos_y*res_x+x];

									data_out += kernel_scalar*data_scalar;
								}
							}
						}
					}
					break;

			default:
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
				for (int y = 0; y < res_y; y++)
				{
					for (int x = 0; x < res_x; x++)
					{
						double &data_out = out.array_data_cartesian_space[y*res_x+x];
						data_out = 0;

						for (int j = -1; j <= 1; j++)
						{
							int pos_y = y+j;

							pos_y -= (pos_y >= res_y ? res_y : 0);
							pos_y += (pos_y < 0 ? res_y : 0);

							assert(pos_y >= 0 && pos_y < res_y);

							for (int i = -1; i <= 1; i++)
							{
								int pos_x = x+i;

								pos_x -= (pos_x >= res_x ? res_x : 0);
								pos_x += (pos_x < 0 ? res_x : 0);

								assert(pos_x >= 0 && pos_x < res_x);

								int idx = (j+1)*3+(i+1);
								assert(idx >= 0 && idx < 9);

								double kernel_scalar = kernel_data[idx];
								double data_scalar = i_array_data.array_data_cartesian_space[pos_y*res_x+pos_x];

								data_out += kernel_scalar*data_scalar;
							}
						}
					}
				}
			}
		}
		else
		{
			std::cerr << "Not yet implemented" << std::endl;
		}
#endif


		out.checkConsistency();
		return out;
	}



	/**
	 * Compute element-wise addition
	 */
	inline
	PlaneData operator+(
			const PlaneData &i_array_data
	)	const
	{
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] =
					array_data_spectral_space[i]+
					i_array_data.array_data_spectral_space[i];

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = array_data_cartesian_space[i] + i_array_data.array_data_cartesian_space[i];

#endif

		out.checkConsistency();

		return out;
	}



	/**
	 * Compute element-wise addition
	 */
	inline
	PlaneData operator+(
			const double i_value
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		{
			requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				out.array_data_spectral_space[i] = array_data_spectral_space[i];

			double scale = resolution[0]*resolution[1];
			out.array_data_spectral_space[0] += i_value*scale;

			out.array_data_cartesian_space_valid = false;
			out.array_data_spectral_space_valid = true;
		}

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = array_data_cartesian_space[i]+i_value;

#endif

		out.checkConsistency();
		return out;
	}




	/**
	 * Compute element-wise addition
	 */
	inline
	PlaneData& operator+=(
			const PlaneData &i_array_data
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] +=
					i_array_data.array_data_spectral_space[i];

		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] += i_array_data.array_data_cartesian_space[i];
#endif


		checkConsistency();
		return *this;
	}



	/**
	 * Compute element-wise addition
	 */
	inline
	PlaneData& operator+=(
			const double i_value
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		double scale = resolution[0]*resolution[1];
		array_data_spectral_space[0] += i_value*scale;

		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] += i_value;

#endif

		checkConsistency();
		return *this;
	}



	/**
	 * Compute multiplication with scalar
	 */
	inline
	PlaneData& operator*=(
			const double i_value
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] *= i_value;

		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] *= i_value;

#endif

		checkConsistency();
		return *this;
	}


	/**
	 * Compute division with scalar
	 */
	inline
	PlaneData& operator/=(
			const double i_value
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] /= i_value;


		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] /= i_value;

#endif

		checkConsistency();
		return *this;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	PlaneData& operator-=(
			const PlaneData &i_array_data
	)
	{
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] -=
					i_array_data.array_data_spectral_space[i];

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;
#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] -= i_array_data.array_data_cartesian_space[i];

#endif


		checkConsistency();
		return *this;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	PlaneData operator-(
			const PlaneData &i_array_data
	)	const
	{
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		PlaneData out(resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] =
					array_data_spectral_space[i]-
					i_array_data.array_data_spectral_space[i];

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]-
					i_array_data.array_data_cartesian_space[i];

#endif


		checkConsistency();
		return out;
	}



	/**
	 * Compute element-wise subtraction
	 */
	inline
	PlaneData operator-(
			const double i_value
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE

//		if (array_data_spectral_space_valid)
//		{
			requestDataInSpectralSpace();

//#pragma error "TOOD: make this depending on rexi_par_sum!"
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				out.array_data_spectral_space[i] = array_data_spectral_space[i];

			double scale = resolution[0]*resolution[1];
			out.array_data_spectral_space[0] -= i_value*scale;

			out.array_data_cartesian_space_valid = false;
			out.array_data_spectral_space_valid = true;
//		}
//		else
//		{
//#if SWEET_THREADING
//#pragma omp parallel for OPENMP_PAR_SIMD
//#endif
//			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
//				out.array_data_cartesian_space[i] =
//						array_data_cartesian_space[i]-i_value;
//
//			out.array_data_cartesian_space_valid = true;
//			out.array_data_spectral_space_valid = false;
//		}
#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]-i_value;

#endif

		out.checkConsistency();
		return out;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	PlaneData valueMinusThis(
			const double i_value
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] = -array_data_spectral_space[i];

		double scale = resolution[0]*resolution[1];
		out.array_data_spectral_space[0] = i_value*scale + out.array_data_spectral_space[0];

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = i_value - array_data_cartesian_space[i];

#endif

		out.checkConsistency();
		return out;
	}

	/**
	 * Invert sign
	 */
	inline
	PlaneData operator-()
	{
		PlaneData out(this->resolution);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] = -array_data_spectral_space[i];

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();

		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = -array_data_cartesian_space[i];

#endif


		out.checkConsistency();
		return out;
	}



	/**
	 * Compute element-wise multiplication
	 */
	inline
	PlaneData operator*(
			const PlaneData &i_array_data	///< this class times i_array_data
	)	const
	{

		// Call as
		// 		data_array.operator*(i_array_data)
		// which is identical to
		//		data_array * i_array_data
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		//This is the actual product result, with the correct N resolution
		PlaneData out(i_array_data.resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE && SWEET_USE_PLANE_SPECTRAL_DEALIASING

		// Array on the left of *, augmented to 3N/2 with zeros on high end spectrum
		// Any previous data in cartesian space will be lost?
		PlaneData u = aliasing_scaleUp();
		assert(u.aliasing_scaled);

		// The input array (right of *) augmented to 3N/2 with zeros on high end spectrum
		PlaneData v = rw_array_data.aliasing_scaleUp();
		assert(v.aliasing_scaled);

		//Convert to cartesian
		u.requestDataInCartesianSpace();
		v.requestDataInCartesianSpace();

		// This array is augmented to 3N/2, since u is
		PlaneData scaled_output(u.resolution, true);

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		//Calculate the product element wise in cartesian space
		for (std::size_t i = 0; i < scaled_output.array_data_cartesian_length; i++)
			scaled_output.array_data_cartesian_space[i] =
					u.array_data_cartesian_space[i]*
					v.array_data_cartesian_space[i];

		scaled_output.array_data_cartesian_space_valid = true;
		scaled_output.array_data_spectral_space_valid = false;

		//Copies the spectrum of the product to the output data_array, which has the correct resolution N
		// As a consequence, all high modes are ignored (beyond N to 3N/2)
		out = scaled_output.aliasing_scaleDown(out.resolution);

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else
		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]*
					i_array_data.array_data_cartesian_space[i];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
#endif

#endif

		out.checkConsistency();
		return out;
	}

	/**
	 * Compute element-wise multiplication in cartesian space
	 * if de-aliasing activated, do a 2/3 truncation in spectrum before and after multiplication
	 *
	 *  *** This is not equivalent to operator* with dealiasing !! It kill more modes than necessary.
	 */
	inline
	PlaneData mult(
			const PlaneData &i_array_data	///< this class times i_array_data
	)	const
	{

		// Call as
		// 		data_array.mult(i_array_data)
		// to represent
		//      data_array*i_array_data

		//This is the actual product result, with the correct de-aliasing
//		PlaneData &rw_array_data = (PlaneData&)i_array_data;
		const PlaneData &data_in1_const= *this;
		const PlaneData &data_in2_const= i_array_data;

		PlaneData out(i_array_data.resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE && SWEET_USE_PLANE_SPECTRAL_DEALIASING
		// Truncate arrays to 2N/3 high end spectrum

		PlaneData data_in1 = data_in1_const;
		data_in1.aliasing_zero_high_modes();

		PlaneData data_in2 = data_in2_const;
		data_in2.aliasing_zero_high_modes();

		out=data_in1*data_in2;
		//Truncate the product, since the high modes could contain alias
		out=out.aliasing_zero_high_modes();
#else
		out=data_in1_const*data_in2_const;
#endif

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
#endif

		out.checkConsistency();
		return out;
	}


	/**
	 * Compute multiplication with a scalar
	 */
	inline
	PlaneData operator*(
			const double i_value
	)	const
	{
		PlaneData out(resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		if (array_data_spectral_space_valid)
		{
#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				out.array_data_spectral_space[i] =
						array_data_spectral_space[i]*i_value;

			out.array_data_cartesian_space_valid = false;
			out.array_data_spectral_space_valid = true;
		}
		else
		{
			assert(array_data_cartesian_space_valid);

#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				out.array_data_cartesian_space[i] =
						array_data_cartesian_space[i]*i_value;

			out.array_data_cartesian_space_valid = true;
			out.array_data_spectral_space_valid = false;
		}

#else
#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]*i_value;
#endif

		out.checkConsistency();
		return out;
	}


	/**
	 * Compute element-wise division
	 */
	inline
	PlaneData operator/(
			const double &i_value
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (array_data_cartesian_space_valid)
		{
#endif
#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				out.array_data_cartesian_space[i] = array_data_cartesian_space[i] / i_value;
#if SWEET_USE_PLANE_SPECTRAL_SPACE
			out.array_data_cartesian_space_valid = true;
			out.array_data_spectral_space_valid = false;
		}
		else
		{
			assert(array_data_spectral_space_valid);

#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				out.array_data_spectral_space[i] = array_data_spectral_space[i] / i_value;

			out.array_data_cartesian_space_valid = false;
			out.array_data_spectral_space_valid = true;
		}
#endif


		out.checkConsistency();
		return out;
	}


	/**
	 * Compute element-wise division
	 */
	inline
	PlaneData operator/(
			const PlaneData &i_array_data
	)	const
	{
		checkConsistency();
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		PlaneData out(this->resolution);
//		out.temporary_data = true;

#if SWEET_USE_PLANE_SPECTRAL_SPACE && SWEET_USE_PLANE_SPECTRAL_DEALIASING

		PlaneData u = aliasing_scaleUp();
		PlaneData v = rw_array_data.aliasing_scaleUp();

		u.requestDataInCartesianSpace();
		v.requestDataInCartesianSpace();

		PlaneData scaled_output(u.resolution, true);

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < scaled_output.array_data_cartesian_length; i++)
			scaled_output.array_data_cartesian_space[i] =
					u.array_data_cartesian_space[i]/
					v.array_data_cartesian_space[i];

		scaled_output.array_data_cartesian_space_valid = true;
		scaled_output.array_data_spectral_space_valid = false;

		out = scaled_output.aliasing_scaleDown(out.resolution);

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else
		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]/
					i_array_data.array_data_cartesian_space[i];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
#endif

#endif

		out.checkConsistency();
		return out;
	}


	friend
	inline
	std::ostream& operator<<(
			std::ostream &o_ostream,
			const PlaneData &i_dataArray
	)
	{
		PlaneData &rw_array_data = (PlaneData&)i_dataArray;

		rw_array_data.requestDataInCartesianSpace();

		for (int y = rw_array_data.resolution[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < rw_array_data.resolution[0]; x++)
			{
				double value = rw_array_data.get(y, x);
//					if (std::abs(value) < 1e-13)
//						value = 0;
				std::cout << value << "\t";
			}
			std::cout << std::endl;
		}

		i_dataArray.checkConsistency();
		return o_ostream;
	}


#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Add scalar to all spectral modes
	 */
	inline
	PlaneData spec_addScalarAll(
			const double &i_value
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			out.array_data_spectral_space[i] = array_data_spectral_space[i] + i_value;
			out.array_data_spectral_space[i+1] = array_data_spectral_space[i+1];
		}

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

		out.checkConsistency();
		return out;
	}
#endif

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Invert all spectral coefficients a+bi --> 1/(a+bi)
	 */
	inline
	PlaneData spec_invert(
	)	const
	{
		PlaneData out(this->resolution);
//		out.temporary_data = true;

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			//get spectral coefficient ar+i*ai
			double ar = array_data_spectral_space[i];
			double ai = array_data_spectral_space[i+1];
			double norm=ar*ar+ai*ai;
			// Calculate 1/(ar+i*ai) and split into real and imag parts
			out.array_data_spectral_space[i] = ar/norm;
			out.array_data_spectral_space[i+1] = -ai/norm;
		}

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

		out.checkConsistency();
		return out;
	}


	inline
	void printSpectrum()	const
	{

		checkConsistency();
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.requestDataInSpectralSpace();

		for (int y = rw_array_data.resolution_spec[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < rw_array_data.resolution_spec[0]; x++)
			{
				double value_re = rw_array_data.spec_getRe(y, x);
				double value_im = rw_array_data.spec_getIm(y, x);
				std::cout << "(" << value_re << ", " << value_im << ")\t";
			}
			std::cout << std::endl;
		}
	}

	inline
	void printSpectrumIndex()	const
	{

		checkConsistency();
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.requestDataInSpectralSpace();

		for (int y = rw_array_data.resolution_spec[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < rw_array_data.resolution_spec[0]; x++)
			{
				double value_re = rw_array_data.spec_getRe(y, x);
				double value_im = rw_array_data.spec_getIm(y, x);
				//if(std::abs(value_re)<1.0e-13)
					//value_re=0.0;
				//if(std::abs(value_im)<1.0e-13)
					//value_im=0.0;

				std::cout << "(" << x << ", "<< y << ", "<< value_re << ", " << value_im << ")\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

	}

	inline
	void printSpectrumNonZero()	const
	{

		checkConsistency();
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.requestDataInSpectralSpace();

		for (int y = rw_array_data.resolution_spec[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < rw_array_data.resolution_spec[0]; x++)
			{
				double value_re = rw_array_data.spec_getRe(y, x);
				double value_im = rw_array_data.spec_getIm(y, x);
				if(value_re*value_re+value_im*value_im>1.0e-13)
					std::cout << "(" << x << ", "<< y << ", "<< value_re << ", " << value_im << ")" <<std::endl;;
			}
			//std::cout << std::endl;
		}
		//std::cout << std::endl;

	}


	inline //TODO - print for each K^2+l^2 the energy in the spectrum
		void printSpectrumEnergy_y()	const
		{

			checkConsistency();
			PlaneData &rw_array_data = (PlaneData&)*this;

			rw_array_data.requestDataInSpectralSpace();
			std::cout << "Energy (sqr)" <<std::endl;

			for (int y = rw_array_data.resolution_spec[1]/2; y >= 0; y--)
			{
			//	for (std::size_t x = 0; x < rw_array_data.resolution_spec[0]; x++)
				std::size_t x=0;
				{
					double value_re = rw_array_data.spec_getRe(y, x);
					double value_im = rw_array_data.spec_getIm(y, x);
					//std::cout << "(" << x << ", " << y << ", "<< value_re << ", "<< value_im*value_im << ")\t";
					double energy=value_re*value_re + value_im*value_im;
					if(energy>1e-14)
						std::cout << "(" << x << ", " << y << ", "<< energy << ")" <<std::endl;
				}
				//std::cout << std::endl;
			}
		}
#endif


	/**
	 * Print data
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool printArrayData(
			int i_precision = 8		///< number of floating point digits
			)	const
	{

		checkConsistency();
		requestDataInCartesianSpace();

		std::ostream &o_ostream = std::cout;

		o_ostream << std::setprecision(i_precision);

		for (int y = resolution[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < resolution[0]; x++)
			{
				o_ostream << get(y, x);

				if (x < resolution[0]-1)
					o_ostream << '\t';
				else
					o_ostream << std::endl;
			}
		}

		checkConsistency();
		return true;
	}



	/**
	 * Write data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_saveData_ascii(
			const char *i_filename,		///< Name of file to store data to
			char i_separator = '\t',	///< separator to use for each line
			int i_precision = 12		///< number of floating point digits
	)
	{

		checkConsistency();
		requestDataInCartesianSpace();

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		for (int y = resolution[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < resolution[0]; x++)
			{
				file << get(y, x);

				if (x < resolution[0]-1)
					file << i_separator;
				else
					file << std::endl;
			}
		}


		checkConsistency();
		return true;
	}


#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Write data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_saveSpectralData_ascii(
			const char *i_filename,		///< Name of file to store data to
			char i_separator = '\t',	///< separator to use for each line
			int i_precision = 12		///< number of floating point digits
	)
	{

		checkConsistency();
		requestDataInSpectralSpace();

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		for (int y = resolution_spec[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < resolution_spec[0]; x++)
			{
				double value_re = spec_getRe(y, x);
				double value_im = spec_getIm(y, x);
				//file << "(" << value_re << ", " << value_im << ")";
				file << sqrt(value_re*value_re+value_im*value_im);
				if (x < resolution_spec[0]-1)
					file << i_separator;
				else
					file << std::endl;
			}
		}

		checkConsistency();
		return true;
	}
#endif


	/**
	 * Write data to VTK file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_saveData_vtk(
			const char *i_filename,		///< Name of file to store data to
			const char *i_title,		///< Title of scalars
			int i_precision = 12		///< number of floating point digits
	)
	{

		checkConsistency();
		requestDataInCartesianSpace();

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		file << "# vtk DataFile Version 2.0" << std::endl;
		file << "Rectangular solid example" << std::endl;
		file << "ASCII" << std::endl;
		file << "DATASET RECTILINEAR_GRID" << std::endl;
		file << "DIMENSIONS " << resolution[0]+1 << " " << resolution[1]+1 << " 1" << std::endl;

		file << "X_COORDINATES " << resolution[0]+1 << " float" << std::endl;
		for (std::size_t x = 0; x < resolution[0]+1; x++)
			file << (double)x/((double)resolution[0]+1) << std::endl;

		file << "Y_COORDINATES " << resolution[1]+1 << " float" << std::endl;
		for (std::size_t y = 0; y < resolution[1]+1; y++)
			file << (double)y/((double)resolution[1]+1) << std::endl;

		file << "Z_COORDINATES 1 float" << std::endl;
		file << "0" << std::endl;


		std::string title = i_title;
		std::replace(title.begin(), title.end(), ' ', '_');
		file << "CELL_DATA " << resolution[0]*resolution[1] << std::endl;
		file << "SCALARS " << title << " float 1" << std::endl;
		file << "LOOKUP_TABLE default" << std::endl;

		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			file << array_data_cartesian_space[i] << std::endl;


		checkConsistency();
		return true;
	}


	/**
	 * Load data from ASCII file.
	 * This is a non-bullet proof implementation, so be careful for invalid file formats.
	 *
	 * New array rows are initialized with a newline.
	 * Each line then has the floating point values stored separated with space ' ' or tabs '\t'
	 *
	 * Note, that the number of values in the ASCII file have to match the resolution of the PlaneData.
	 *
	 * \return true if data was successfully read
	 */
	bool file_loadData(
			const char *i_filename,		///< Name of file to load data from
			bool i_binary_data = false	///< load as binary data (disabled per default)
	)
	{
		checkConsistency();
		if (i_binary_data)
		{
			std::ifstream file(i_filename, std::ios::binary);

			if (!file)
			{
				std::cerr << "Failed to open file " << i_filename << std::endl;
				exit(-1);
			}
			file.seekg(0, std::ios::end);
			std::size_t size = file.tellg();
			file.seekg(0, std::ios::beg);


			std::size_t expected_size = sizeof(double)*resolution[0]*resolution[1];

			if (size != expected_size)
			{
				std::cerr << "Error while loading data from file " << i_filename << ":" << std::endl;
				std::cerr << "Size of file " << size << " does not match expected size of " << expected_size << std::endl;
				exit(-1);
			}

			if (!file.read((char*)array_data_cartesian_space, expected_size))
			{
				std::cerr << "Error while loading data from file " << i_filename << std::endl;
				exit(1);
			}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			array_data_cartesian_space_valid = true;
			array_data_spectral_space_valid = false;
#endif
			return true;
		}
		std::ifstream file(i_filename);

		for (std::size_t row = 0; row < resolution[1]; row++)
		{
			std::string line;
			std::getline(file, line);
			if (!file.good())
			{
				std::cerr << "Failed to read data from file " << i_filename << " in line " << row << std::endl;
				return false;
			}

			std::size_t last_pos = 0;
			std::size_t col = 0;
			for (std::size_t pos = 0; pos < line.size()+1; pos++)
			{
				if (pos < line.size())
					if (line[pos] != '\t' && line[pos] != ' ')
						continue;

				std::string strvalue = line.substr(last_pos, pos-last_pos);

				double i_value = atof(strvalue.c_str());

				set(resolution[1]-row-1, col, i_value);

				col++;
				last_pos = pos+1;
		    }

			if (col < resolution[0])
			{
				std::cerr << "Failed to read data from file " << i_filename << " in line " << row << ", column " << col << std::endl;
				return false;
			}
		}

		checkConsistency();
		return true;
	}
};


/**
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
PlaneData operator*(
		const double i_value,
		const PlaneData &i_array_data
)
{
	i_array_data.checkConsistency();
	return ((PlaneData&)i_array_data)*i_value;
}

/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 *
 * Otherwise, we'd have to write it as arrayData-1.5
 *
 */
inline
static
PlaneData operator-(
		const double i_value,
		const PlaneData &i_array_data
)
{
	i_array_data.checkConsistency();
	return ((PlaneData&)i_array_data).valueMinusThis(i_value);
//	return -(((PlaneData&)i_array_data).operator-(i_value));
}
/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */
inline
static
PlaneData operator+(
		const double i_value,
		const PlaneData &i_array_data
)
{
	i_array_data.checkConsistency();
	return ((PlaneData&)i_array_data)+i_value;
}

#endif /* SRC_DATAARRAY_HPP_ */
