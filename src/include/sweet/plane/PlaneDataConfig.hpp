/*
 * PlaneDataConfig.hpp
 *
 *  Created on: 17 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHSETUP_HPP_
#define SPHSETUP_HPP_


#include <libmath/shtns_inc.hpp>
#include <fftw3.h>
#include <iostream>
#include <sweet/sweetmath.hpp>
#include <stdlib.h>
#include <iostream>

#if SWEET_THREADING
#	include <omp.h>
#endif



/*
 * This precompiler directive should !!! NEVER !!! be used
 */
#ifndef SWEET_USE_PLANE_SPECTRAL_SPACE
	#define SWEET_USE_PLANE_SPECTRAL_SPACE	1
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

private:

#if SWEET_USE_LIBFFT

	/// number of spectral modes
	std::size_t spectral_modes[2];

	/// allocated size for spectral data
	std::size_t spectral_data_size[2];

	/// iteration end for updating data in spectrum:
	/// This is different in case of anti-aliasing
	std::size_t spectral_data_iteration_end[2];

	/// number of complex-valued data
	std::size_t spectral_array_data_number_of_elements;


	/*
	 * FFTW related stuff
	 */
	fftw_plan	fftw_plan_forward;
	fftw_plan	fftw_plan_backward;

	double fftw_backward_scale_factor;



	/// allocated size for spectral data in case of complex data in physical space
	std::size_t spectral_complex_data_size[2];
	std::size_t spectral_complex_data_iteration_end[2];
	std::size_t spectral_complex_array_data_number_of_elements;

	fftw_plan	fftw_plan_complex_forward;
	fftw_plan	fftw_plan_complex_backward;

#endif

	bool initialized;



public:
	static
	bool loadWisdom()
	{
#if SWEET_REXI_THREAD_PARALLEL_SUM

		// Threaded time parallel sum
		std::cout << "Using REXI parallel sum, hence using only single FFT thread" << std::endl;

#else

#if SWEET_THREADING
		// initialise FFTW with spatial parallelization
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
	PlaneDataConfig()
	{
		physical_res[0] = 0;
		physical_res[1] = 0;
#if SWEET_USE_LIBFFT
		spectral_modes[0] = 0;
		spectral_modes[1] = 0;
#endif

		initialized = false;

//		refCounter()++;
	}



private:
	void setup_internal_data()
	{
		if (initialized)
		{
			cleanup();
		}

		physical_data_size[0] = physical_res[0];
		physical_data_size[1] = physical_res[1];

		physical_array_data_number_of_elements = physical_res[0]*physical_res[1];

#if SWEET_USE_LIBFFT
		if (	physical_res[0] < spectral_modes[0]	||
				physical_res[1] < spectral_modes[1]
		)
		{
			std::cerr << "Lower physical resolution than spectral resolution not supported!" << std::endl;
			exit(-1);
		}

		assert(physical_res[0] > 0);
		assert(physical_res[1] > 0);

		assert(spectral_modes[0] > 0);
		assert(spectral_modes[1] > 0);

		// real-to-complex storage representation
		spectral_data_size[0] = spectral_modes[0]/2+1;
		spectral_data_size[1] = spectral_modes[1];

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
//#error "TODO: check all that stuff"
		// TODO: setup dealiasing limitations correctly
		spectral_data_iteration_end[0] = spectral_data_size[0]*2/3-1;
		spectral_data_iteration_end[1] = spectral_data_size[1]*2/3-1;
#else
		spectral_data_iteration_end[0] = spectral_data_size[0];
		spectral_data_iteration_end[1] = spectral_data_size[1];
#endif

		spectral_array_data_number_of_elements = spectral_data_size[0]*spectral_data_size[1];


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
			fftw_backward_scale_factor = 1.0/((double)(spectral_modes[0]*spectral_modes[1]));
		}


		/*
		 * COMPLEX PHYSICAL SPACE DATA
		 */
		{
			spectral_complex_array_data_number_of_elements = spectral_modes[0]*spectral_modes[1];
			spectral_complex_data_iteration_end[0] = spectral_modes[0];
			spectral_complex_data_iteration_end[1] = spectral_modes[1];

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
						(!wisdom_loaded ? FFTW_PRESERVE_INPUT | fftw_estimate_plan : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
//						(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
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
						(!wisdom_loaded ? FFTW_PRESERVE_INPUT | fftw_estimate_plan : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
//						(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
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
	}



#if SWEET_USE_LIBFFT
	void fft_physical_to_spectral(
			double *i_physical_data,
			std::complex<double> *o_spectral_data
	)
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
	)
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
	void setup(
			std::size_t i_physical_res[2],
			std::size_t i_spectral_modes[2]
	)
	{
		if (i_spectral_modes[0] <= 0)
		{
			setupAutoSpectralSpace(
					i_physical_res[0],
					i_physical_res[1]
				);
		}
		else if (i_physical_res[0] > 0)
		{
			setup(	i_physical_res[0],
					i_physical_res[1],
					i_spectral_modes[0],
					i_spectral_modes[1]
				);
		}
		else
		{
			std::cerr << "No resolution selected" << std::endl;
			exit(-1);
		}
	}

#if 0
public:
	int& refCounter()
	{
		static int ref_counter = 0;
		return ref_counter;
	}

#endif


public:
	void setupAutoSpectralSpace(
			std::size_t i_physical_res[2]
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
		spectral_modes[0] = physical_res[0]/2+1;
		spectral_modes[1] = physical_res[1];

#else

	#if SWEET_USE_LIBFFT
		spectral_modes[0] = physical_res[0];
		spectral_modes[1] = physical_res[1];
	#endif

#endif

		setup_internal_data();
	}



	void cleanup()
	{

#if SWEET_USE_LIBFFT
//		assert(refCounter() >= 0);

//		if (refCounter() == 0)
		{
			fftw_destroy_plan(fftw_plan_forward);
			fftw_destroy_plan(fftw_plan_backward);

			fftw_destroy_plan(fftw_plan_complex_forward);
			fftw_destroy_plan(fftw_plan_complex_backward);

#if SWEET_THREADING
			fftw_cleanup_threads();
#endif
			fftw_cleanup();
		}
#endif
	}



	~PlaneDataConfig()
	{
//		refCounter()--;
		cleanup();
	}
};



#endif /* SPHSETUP_HPP_ */
