/*
 * DataArray.hpp
 *
 *  Created on: 28 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_DATAARRAY_HPP_
#define SRC_DATAARRAY_HPP_

#include <cassert>
#include <cstddef>
#include <cmath>
#include <memory>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <utility>
#include <limits>
#include <sweet/openmp_helper.hpp>

#ifndef SWEET_USE_SPECTRAL_SPACE
	#define SWEET_USE_SPECTRAL_SPACE	1
#endif

#ifndef SWEET_USE_SPECTRAL_DEALIASING
	#define SWEET_USE_SPECTRAL_DEALIASING 1
#endif

#if SWEET_USE_SPECTRAL_SPACE
#	include <fftw3.h>
#endif

#include <omp.h>



/**
 * Data array and operator support.
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
template <int D>
class DataArray
{
public:
	/**
	 * global size of allocated array
	 * (x,y[,z])
	 */
	std::size_t resolution[D];

	/**
	 * local data in cartesian space
	 */
	std::size_t array_data_cartesian_length;
	double *array_data_cartesian_space;

#if SWEET_USE_SPECTRAL_SPACE
	std::size_t resolution_spec[D];
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
	std::size_t range_start[D];
	std::size_t range_end[D];
	std::size_t range_size[D];

#if SWEET_USE_SPECTRAL_SPACE
	std::size_t range_spec_start[D];
	std::size_t range_spec_end[D];
	std::size_t range_spec_size[D];
#else
	int kernel_size = -1;
	double *kernel_data = nullptr;
#endif

	/**
	 * temporary data?
	 *
	 * Temporary data can be created if e.g. operators are evaluated:
	 * h = hu/h
	 *
	 * This first creates a temporary DataArray to compute hu/h.
	 *
	 * This is then followed by an assignment of this data to h.
	 */
	bool temporary_data;

/*

private:
	// http://stackoverflow.com/questions/124856/how-do-i-prevent-a-class-from-being-allocated-via-the-new-operator-id-like
	// Prevent heap allocation
	void *operator new(std::size_t);
	void *operator new[](std::size_t);
	void operator delete(void*);
	void operator delete[](void*);
*/

	/**
	 * allocator which allocated memory blocks aligned at 128 byte boundaries
	 */
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
	 * allocate buffers
	 *
	 * The size is given in array_data_cartesian_length and array_data_spectral_length
	 */
private:
	void p_allocate_buffers(
			bool i_first_touch_initialize = true	///< true: initialize the data buffers with dummy data for first touch policy of page allocation on shared-memory systems
	)
	{
		free(array_data_cartesian_space);
		array_data_cartesian_space = alloc_aligned_mem<double>(array_data_cartesian_length*sizeof(double));

#if SWEET_USE_SPECTRAL_SPACE
		free(array_data_spectral_space);
		array_data_spectral_space = alloc_aligned_mem<double>(array_data_spectral_length*sizeof(double));
#endif

		if (i_first_touch_initialize)
		{
			// use parallel setup for first touch policy!
			#pragma omp parallel for
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				array_data_cartesian_space[i] = -12345;	// dummy data

#if SWEET_USE_SPECTRAL_SPACE
			#pragma omp parallel for
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
		const std::size_t i_resolution[D]		///< requested resolution
	)
	{
		// check if the resolution is maybe already the same as requested
		int i = 0;
		for (; i < D; i++)
			if (resolution[i] != i_resolution[i])
				break;

		if (i == D)
			return;

		array_data_cartesian_length = 1;
		for (int i = 0; i < D; i++)
		{
			array_data_cartesian_length *= i_resolution[i];

			resolution[i] = i_resolution[i];
			range_start[i] = 0;
			range_end[i] = i_resolution[i];
			range_size[i] = range_end[i]-range_start[i];
		}

#if SWEET_USE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = false;

		array_data_spectral_length = array_data_cartesian_length/i_resolution[0];	/// see FFTW documentation for allocation of memory buffers
		array_data_spectral_length *= 2*(i_resolution[0]/2+1);

		for (int i = 0; i < D; i++)
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

		p_allocate_buffers();
	}



	/**
	 * prohobit empty initialization by making this method private
	 */
private:
	DataArray()
#if 1
	{}
#else
	{
		resolution[0] = 0;
		resolution[1] = 0;

		temporary_data = false;

		array_data_cartesian_space = nullptr;
		array_data_cartesian_length = 0;

#if SWEET_USE_SPECTRAL_SPACE
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
	 * 	DataArray<2> tmp_h = h;
	 * 	DataArray<2> tmp_h2(h);
	 *
	 * Duplicate all data
	 */
	DataArray(
			const DataArray<D> &i_dataArray
	)	:
		array_data_cartesian_space(nullptr),
#if SWEET_USE_SPECTRAL_SPACE
		array_data_spectral_space(nullptr),
		aliasing_scaled(i_dataArray.aliasing_scaled),
#endif

		temporary_data(false)
	{
		for (int i = 0; i < D; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];

#if SWEET_USE_SPECTRAL_SPACE
			resolution_spec[i] = i_dataArray.resolution_spec[i];

			range_spec_start[i] = i_dataArray.range_spec_start[i];
			range_spec_end[i] = i_dataArray.range_spec_end[i];
			range_spec_size[i] = i_dataArray.range_spec_size[i];
#endif
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;
#if SWEET_USE_SPECTRAL_SPACE
		array_data_spectral_length = i_dataArray.array_data_spectral_length;
#endif

		p_allocate_buffers(false);

#if SWEET_USE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = i_dataArray.array_data_cartesian_space_valid;
		if (array_data_cartesian_space_valid)
#endif
		{
			// use parallel copy for first touch policy!
#pragma omp parallel for
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				array_data_cartesian_space[i] = i_dataArray.array_data_cartesian_space[i];
		}

#if SWEET_USE_SPECTRAL_SPACE
		array_data_spectral_length = i_dataArray.array_data_spectral_length;
		array_data_spectral_space = alloc_aligned_mem<double>(array_data_spectral_length*sizeof(double));
		array_data_spectral_space_valid = i_dataArray.array_data_spectral_space_valid;

		if (array_data_spectral_space_valid)
		{
			// use parallel copy for first touch policy!
#pragma omp parallel for
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				array_data_spectral_space[i] = i_dataArray.array_data_spectral_space[i];

		}

		auto fft_ptr = *fftGetSingletonPtr();
		fft_ptr->ref_counter++;
#endif
	}



public:
	/**
	 * Move constructor
	 */
	DataArray(
			DataArray<D> &&i_dataArray
	)	:
		array_data_cartesian_space(nullptr),
#if SWEET_USE_SPECTRAL_SPACE
		array_data_spectral_space(nullptr),
		aliasing_scaled(i_dataArray.aliasing_scaled),
#endif
		temporary_data(false)
	{
		for (int i = 0; i < D; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];

#if SWEET_USE_SPECTRAL_SPACE
			resolution_spec[i] = i_dataArray.resolution_spec[i];

			range_spec_start[i] = i_dataArray.range_spec_start[i];
			range_spec_end[i] = i_dataArray.range_spec_end[i];
			range_spec_size[i] = i_dataArray.range_spec_size[i];
#endif
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;
		array_data_cartesian_space = i_dataArray.array_data_cartesian_space;
		i_dataArray.array_data_cartesian_space = nullptr;

#if SWEET_USE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = i_dataArray.array_data_cartesian_space_valid;
		i_dataArray.array_data_cartesian_space_valid = false;

		array_data_spectral_length = i_dataArray.array_data_spectral_length;
		array_data_spectral_space = i_dataArray.array_data_spectral_space;
		i_dataArray.array_data_spectral_space = nullptr;
		array_data_spectral_space_valid = i_dataArray.array_data_spectral_space_valid;
		i_dataArray.array_data_spectral_space_valid = false;

		auto fft_ptr = *fftGetSingletonPtr();
		fft_ptr->ref_counter++;
#endif
	}



	/**
	 * default constructor
	 */
public:
	DataArray(
		const std::size_t i_resolution[D]	///< size of array
	)	:
		array_data_cartesian_space(nullptr),
#if SWEET_USE_SPECTRAL_SPACE
		array_data_spectral_space(nullptr),
		aliasing_scaled(false),
#endif
		temporary_data(false)
	{
		for (int i = 0; i < D; i++)
		{
			resolution[i] = 0;
#if SWEET_USE_SPECTRAL_SPACE
			if (i_resolution[i] & 1)
			{
				std::cerr << "Sorry - only even resolution supported for spectral space, makes the life much easier!" << std::endl;
				exit(1);
			}
#endif
		}

		p_request_buffers_with_resolution(i_resolution);

#if SWEET_USE_SPECTRAL_SPACE
		// initialize fft if not yet done
		fftTestAndInit(*this);
#endif
	}



	~DataArray()
	{
		free(array_data_cartesian_space);

#if SWEET_USE_SPECTRAL_SPACE
		free(array_data_spectral_space);

		{
			auto fft_ptr = *fftGetSingletonPtr();
			fft_ptr->ref_counter--;

			assert(fft_ptr->ref_counter >= 0);
			if (fft_ptr->ref_counter == 0)
			{
				delete *fftGetSingletonPtr();
				*fftGetSingletonPtr() = nullptr;

				// also free the aliasing stuff
//				if (fft_ptr->ref_counter == 0)
				{
					delete *fftAliasingGetSingletonPtr();
					*fftAliasingGetSingletonPtr() = nullptr;
				}
			}
		}

		if (aliasing_scaled)
		{
			auto fft_ptr = *fftAliasingGetSingletonPtr();
			fft_ptr->ref_counter--;

#if 0
			assert(fft_ptr->ref_counter >= 0);
			if (fft_ptr->ref_counter == 0)
			{
				delete *fftAliasingGetSingletonPtr();
				*fftAliasingGetSingletonPtr() = nullptr;
			}
#endif
		}
#else
		free(kernel_data);
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

#if SWEET_USE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;
#endif

		array_data_cartesian_space[
							(j-range_start[1])*range_size[0]+
							(i-range_start[0])
						] = i_value;
	}



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
	}


	inline
	void setAll(
			double i_value
	)
	{
		// TODO: implement this in spectral space!

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] = i_value;

#if SWEET_USE_SPECTRAL_SPACE
		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;
#endif
	}



#if SWEET_USE_SPECTRAL_SPACE==1

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
		return array_data_spectral_space[idx*2+1];
	}


	inline
	void spec_set(
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
	}



	inline
	void spec_setAll(
			double i_value_re,
			double i_value_im
	)
	{
#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			array_data_spectral_space[i+0] = i_value_re;
			array_data_spectral_space[i+1] = i_value_im;
		}

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;
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
	void spec_diff_set(
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
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void spec_set_spectrum_with_ampl_and_phase(
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

		spec_set_spectrum(j, i, c, s);
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void spec_set_spectrum(
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
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void spec_set_spectrum_A(
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
	}


	/**
	 * Set the spectrum of a frequency with amplitude and phase.
	 *
	 * The amplitude specifies the magnitude in real space.
	 * The phase specifies the shift from one amplitude to another one.
	 */
	inline
	void spec_set_spectrum_B(
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
	}



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
		FFTWSingletonClass(
				DataArray<D> &i_dataArray
		)	:
			ref_counter(0)
		{
#if SWEET_THREADING
			// support threading
			fftw_init_threads();

		    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

			plan_backward_output_length = i_dataArray.array_data_cartesian_length;
			plan_forward_output_length = i_dataArray.array_data_spectral_length;

			double *data_cartesian = alloc_aligned_mem<double>(i_dataArray.array_data_cartesian_length*sizeof(double));
#pragma omp parallel for
			for (std::size_t i = 0; i < i_dataArray.array_data_cartesian_length; i++)
				data_cartesian[i] = -123;	// dummy data

			double *data_spectral = alloc_aligned_mem<double>(i_dataArray.array_data_spectral_length*sizeof(double));
#pragma omp parallel for
			for (std::size_t i = 0; i < i_dataArray.array_data_spectral_length; i++)
				data_spectral[i] = -123;	// dummy data

			// allow to search for a good plan only for 60 seconds
//			fftw_set_timelimit(60);

			static const char *wisdom_filename = "FFTW_WISDOM.cached";

			int wisdom_plan_loaded = fftw_import_wisdom_from_filename(wisdom_filename);
			if (wisdom_plan_loaded > 0)
			{
				std::cout << "Successfully loaded FFTW wisdom from file " << wisdom_filename << std::endl;
			}

			plan_forward =
					fftw_plan_dft_r2c_2d(
						i_dataArray.resolution[1],	// n0 = ny
						i_dataArray.resolution[0],	// n1 = nx
						data_cartesian,
						(fftw_complex*)data_spectral,
						FFTW_PRESERVE_INPUT
					);

			if (plan_forward == nullptr)
			{
				std::cerr << "Failed to create forward plan for fftw" << std::endl;
				exit(-1);
			}

			plan_backward =
					fftw_plan_dft_c2r_2d(
						i_dataArray.resolution[1],	// n0 = ny
						i_dataArray.resolution[0],	// n1 = nx
						(fftw_complex*)data_spectral,
						data_cartesian,
						0
					);

			if (plan_backward == nullptr)
			{
				std::cerr << "Failed to create backward plan for fftw" << std::endl;
				exit(-1);
			}

			// always store plans - maybe they got extended with another one
//			if (wisdom_plan_loaded == 0)
			fftw_export_wisdom_to_filename(wisdom_filename);


			free(data_cartesian);
			free(data_spectral);
		}

		void fft_forward(
				DataArray<D> &io_dataArray
		)
		{
			fftw_execute_dft_r2c(plan_forward, io_dataArray.array_data_cartesian_space, (fftw_complex*)io_dataArray.array_data_spectral_space);
		}

		void fft_backward(
				DataArray<D> &io_dataArray
		)
		{
			fftw_execute_dft_c2r(plan_backward, (fftw_complex*)io_dataArray.array_data_spectral_space, io_dataArray.array_data_cartesian_space);
			// spectral data is not valid anymore, since c2r is destructive!
			io_dataArray.array_data_spectral_space_valid = false;

			double scale = (1.0/(double)plan_backward_output_length);
#pragma omp parallel for OPENMP_SIMD
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


private:
	FFTWSingletonClass** fftGetSingletonPtr()	const
	{
		static FFTWSingletonClass *fftw_singleton_data = nullptr;
		return &fftw_singleton_data;
	}


private:
	FFTWSingletonClass* fftTestAndInit(
		DataArray<D> &i_dataArray
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
	FFTWSingletonClass** fftAliasingGetSingletonPtr()	const
	{
		static FFTWSingletonClass *fftw_singleton_data = nullptr;
		return &fftw_singleton_data;
	}

private:
	FFTWSingletonClass* fftAliasingTestAndInit(
		DataArray<D> &i_dataArray
	)	const
	{
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
#endif

public:
	inline
	void requestDataInSpectralSpace()	const
	{
#if SWEET_USE_SPECTRAL_SPACE==0
		std::cerr << "requestDataInSpectralSpace: spectral space is disabled" << std::endl;
		exit(-1);
#else

		if (array_data_spectral_space_valid)
			return;		// nothing to do

		if (!array_data_cartesian_space_valid)
		{
			std::cerr << "Spectral data not available! Is this maybe a non-initialized operator?" << std::endl;
			assert(false);
			exit(1);
		}

		DataArray<D> *rw_array_data = (DataArray<D>*)this;

		if (aliasing_scaled)
			(*fftAliasingGetSingletonPtr())->fft_forward(*rw_array_data);
		else
			(*fftGetSingletonPtr())->fft_forward(*rw_array_data);

		rw_array_data->array_data_spectral_space_valid = true;
#endif
	}


	inline
	void requestDataInCartesianSpace()	const
	{
#if SWEET_USE_SPECTRAL_SPACE==1
		if (array_data_cartesian_space_valid)
			return;		// nothing to do

		assert(array_data_spectral_space_valid == true);

		DataArray<D> *rw_array_data = (DataArray<D>*)this;

		if (aliasing_scaled)
		{
			(*fftAliasingGetSingletonPtr())->fft_backward(*rw_array_data);
		}
		else
			(*fftGetSingletonPtr())->fft_backward(*rw_array_data);

		rw_array_data->array_data_cartesian_space_valid = true;
#endif
	}


	inline
	DataArray<D> return_one_if_positive()
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] > 0 ? 1 : 0);

#if SWEET_USE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif
		return out;
	}



	inline
	DataArray<D> return_value_if_positive()	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] > 0 ? array_data_cartesian_space[i] : 0);

#if SWEET_USE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif
		return out;
	}


	inline
	DataArray<D> return_one_if_negative()	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] < 0 ? 1 : 0);

#if SWEET_USE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif
		return out;
	}


	inline
	DataArray<D> return_value_if_negative()	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = (array_data_cartesian_space[i] < 0 ? array_data_cartesian_space[i] : 0);

#if SWEET_USE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
#endif
		return out;
	}



	/**
	 * return true, if any value is infinity
	 */
	bool reduce_all_finite() const
	{
		requestDataInCartesianSpace();

		bool isallfinite = true;
#pragma omp parallel for OPENMP_SIMD reduction(&&:isallfinite)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			isallfinite = isallfinite && std::isfinite(array_data_cartesian_space[i]);

		return isallfinite;
	}



	/**
	 * return the maximum of all absolute values
	 */
	double reduce_maxAbs()	const
	{
		requestDataInCartesianSpace();

		double maxabs = -1;
#pragma omp parallel for OPENMP_SIMD reduction(max:maxabs)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			maxabs = std::max(maxabs, std::abs(array_data_cartesian_space[i]));

		return maxabs;
	}


	/**
	 * reduce to root mean square
	 */
	double reduce_rms()
	{
		requestDataInCartesianSpace();

		double sum = 0;
#pragma omp parallel for OPENMP_SIMD reduction(+:sum)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += array_data_cartesian_space[i]*array_data_cartesian_space[i];

		sum = std::sqrt(sum/double(array_data_cartesian_length));
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

#pragma omp parallel for reduction(+:sum,c)
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
		return sum;
	}



	/**
	 * return the maximum of all absolute values
	 */
	double reduce_max()	const
	{
		requestDataInCartesianSpace();

		double maxvalue = -std::numeric_limits<double>::max();
#pragma omp parallel for OPENMP_SIMD reduction(max:maxvalue)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			maxvalue = std::max(maxvalue, array_data_cartesian_space[i]);

		return maxvalue;
	}


	/**
	 * return the maximum of all absolute values
	 */
	double reduce_min()	const
	{
		requestDataInCartesianSpace();

		double minvalue = std::numeric_limits<double>::max();
#pragma omp parallel for OPENMP_SIMD reduction(min:minvalue)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			minvalue = std::min(minvalue, array_data_cartesian_space[i]);

		return minvalue;
	}


	/**
	 * return the maximum of all absolute values
	 */
	double reduce_sum()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
#pragma omp parallel for OPENMP_SIMD reduction(+:sum)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += array_data_cartesian_space[i];

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
#pragma omp parallel for reduction(+:sum,c)
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

		return sum;
	}

	/**
	 * return the maximum of all absolute values
	 */
	double reduce_norm1()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
#pragma omp parallel for OPENMP_SIMD reduction(+:sum)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += std::abs(array_data_cartesian_space[i]);

		return sum;
	}

	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double reduce_norm1_quad()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
		double c = 0;
#pragma omp parallel for reduction(+:sum,c)
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

		return sum;
	}


	/**
	 * return the maximum of all absolute values
	 */
	double reduce_norm2()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
#pragma omp parallel for OPENMP_SIMD reduction(+:sum)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			sum += array_data_cartesian_space[i]*array_data_cartesian_space[i];

		return std::sqrt(sum);
	}


	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double reduce_norm2_quad()	const
	{
		requestDataInCartesianSpace();

		double sum = 0.0;
		double c = 0.0;

#pragma omp parallel for reduction(+:sum,c)
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

		return std::sqrt(sum);
	}


#if SWEET_USE_SPECTRAL_SPACE

	/**
	 * return the maximum of all absolute values
	 */
	double reduce_spec_maxAbs()	const
	{
		requestDataInSpectralSpace();

		double maxabs = -1;
#pragma omp parallel for OPENMP_SIMD reduction(max:maxabs)
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			double re = array_data_spectral_space[i];
			double im = array_data_spectral_space[i+1];
			maxabs = std::max(maxabs, std::sqrt(re*re + im*im));
		}

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

		return nom/denom;
	}
#endif


public:
	template <int S>
	void setup_kernel(
			const double i_kernel_array[S][S],
			double i_scale = 1.0
	)
	{
#if SWEET_USE_SPECTRAL_SPACE == 0

		kernel_size = S;
		kernel_data = alloc_aligned_mem<double>(sizeof(double)*S*S);
		for (int y = 0; y < S; y++)
			for (int x = 0; x < S; x++)
				kernel_data[y*S+x] = i_kernel_array[S-1-y][x];

		for (int i = 0; i < S*S; i++)
			kernel_data[i] *= i_scale;

#else

		double inv_kernel_array[S][S];

		for (int j = 0; j < S; j++)
			for (int i = 0; i < S; i++)
				inv_kernel_array[j][i] = i_kernel_array[j][S-i-1]*i_scale;

		assert(D == 2);

		// assure symmetric kernel
		assert((S & 1) == 1);

		// radius of kernel (half size)
		std::size_t R = S>>1;

		setAll(0);

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
	}


public:
	template <int S>
	void setup_kernel(
			const double i_kernel_array[S][S][S]
	)
	{
		assert(D == 3);
		assert(false);//TODO
	}



public:
	/**
	 * assignment operator
	 */
	DataArray<D> &operator=(double i_value)
	{
		setAll(i_value);
		return *this;
	}


public:
	/**
	 * assignment operator
	 */
	DataArray<D> &operator=(int i_value)
	{
		setAll(i_value);
		return *this;
	}


public:
	/**
	 * assignment operator
	 *
	 * hasdfasdf = h*hasdf;
	 */
	DataArray<D> &operator=(
			const DataArray<D> &i_dataArray
	)
	{
#if SWEET_USE_SPECTRAL_SPACE
		aliasing_scaled = i_dataArray.aliasing_scaled;
#endif

		for (int i = 0; i < D; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;

#if SWEET_USE_SPECTRAL_SPACE
		if (i_dataArray.array_data_cartesian_space_valid)
		{
			array_data_cartesian_space_valid = true;
#endif
			/**
			 * If this data was generated based on temporary data sets (e.g. via h = hu/u), then only swap pointers.
			 */
			if (i_dataArray.temporary_data)
			{
				std::swap(array_data_cartesian_space, ((DataArray<D> &)i_dataArray).array_data_cartesian_space);
			}
			else
			{
				#pragma omp parallel for OPENMP_SIMD
				for (std::size_t i = 0; i < array_data_cartesian_length; i++)
					array_data_cartesian_space[i] = i_dataArray.array_data_cartesian_space[i];
			}

#if SWEET_USE_SPECTRAL_SPACE
		}
		else
		{
			array_data_cartesian_space_valid = false;
		}

		array_data_spectral_length = i_dataArray.array_data_spectral_length;

		if (i_dataArray.array_data_spectral_space_valid)
		{
			array_data_spectral_space_valid = true;
			if (i_dataArray.temporary_data)
			{
				std::swap(array_data_spectral_space, ((DataArray<D> &)i_dataArray).array_data_spectral_space);
			}
			else
			{
				#pragma omp parallel for OPENMP_SIMD
				for (std::size_t i = 0; i < array_data_spectral_length; i++)
					array_data_spectral_space[i] = i_dataArray.array_data_spectral_space[i];
			}
		}
		else
		{
			array_data_spectral_space_valid = false;
		}
#endif

		return *this;
	}


	/**
	 * Apply a linear operator given by this class to the input data array.
	 */
	inline
	DataArray<D> operator()(
			const DataArray<D> &i_array_data
	)	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#pragma omp parallel for OPENMP_SIMD
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

		int res_x = resolution[0];
		int res_y = resolution[1];

		if (kernel_size == 3)
		{
#pragma omp parallel for
			for (int y = 0; y < res_y; y++)
			{
				for (int x = 0; x < res_x; x++)
				{
					double &data_out = out.array_data_cartesian_space[y*res_x+x];
					data_out = 0;

					for (int j = -1; j <= 1; j++)
					{
						int pos_y = (y+j+res_y) % res_y;
						assert(pos_y >= 0 && pos_y < res_y);

						for (int i = -1; i <= 1; i++)
						{
							int pos_x = (x+i+res_x) % res_x;
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
		else
		{
			std::cerr << "Not yet implemented" << std::endl;
		}
#endif

		return out;
	}


#if SWEET_USE_SPECTRAL_SPACE
	/**
	 * Invert the application of a linear operator in spectral space.
	 * The operator is given in i_array_data
	 */
	inline
	DataArray<D> spec_div_element_wise(
			const DataArray<D> &i_array_data,	///< operator
			double i_denom_zeros_scalar = 0.0,
			double i_tolerance = 0.0			///< set tolerance to 0, since we setup the values in the spectral operator directly
	)	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		// only makes sense, if this is an operator created in spectral space
		assert(i_array_data.array_data_spectral_space_valid == true);

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

		// determine maximum value for tolerance
		double max_value = i_array_data.reduce_spec_maxAbs();
		i_tolerance *= max_value;
		i_tolerance *= (resolution[0]+resolution[1]);	// the larger the matrix, the less the accuracy

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			double ar = array_data_spectral_space[i];
			double ai = array_data_spectral_space[i+1];
			double br = i_array_data.array_data_spectral_space[i];
			double bi = i_array_data.array_data_spectral_space[i+1];

			double den = (br*br+bi*bi);

			if (std::abs(den) <= i_tolerance)
			{
				// For Laplace solution, this is the integration constant C
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

		return out;
	}
#endif


#if SWEET_USE_SPECTRAL_DEALIASING && SWEET_USE_SPECTRAL_SPACE

	/**
	 * scale the solution to avoid aliasing effects for handling non-linear terms
	 */
	inline
	DataArray<D> aliasing_scaleUp(
			std::size_t *i_new_resolution = nullptr
	)	const
	{
		std::size_t new_resolution[D];

		if (i_new_resolution == nullptr)
		{
			for (int i = 0; i < D; i++)
			{
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
			for (int i = 0; i < D; i++)
				new_resolution[i] = i_new_resolution[i];
		}

		DataArray<D> out(new_resolution);
//		out.temporary_data = true;
		out.aliasing_scaled = true;

		fftAliasingTestAndInit(out);

		requestDataInSpectralSpace();

		if (D != 2)
		{
			std::cerr << "TODO: Only 2D so far supported" << std::endl;
			exit(-1);
		}

		out.spec_setAll(0, 0);

		// TODO: this does not work once distributed memory is activated
		#pragma omp parallel for
		for (std::size_t j = 0; j < resolution_spec[1]/2; j++)
		{
			// lower quadrant
			memcpy(
					out.array_data_spectral_space+(j*out.range_spec_size[0])*2,
					array_data_spectral_space+(j*range_spec_size[0])*2,
					sizeof(double)*range_spec_size[0]*2
				);

			// top quadrant
			memcpy(
					out.array_data_spectral_space+((out.resolution_spec[1]-resolution_spec[1]/2+j)*out.range_spec_size[0])*2,
					array_data_spectral_space+((resolution_spec[1]-resolution_spec[1]/2+j)*range_spec_size[0])*2,
					sizeof(double)*range_spec_size[0]*2
				);
		}

		double scale = ((double)new_resolution[0]*(double)new_resolution[1])/((double)resolution[0]*(double)resolution[1]);

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < out.array_data_spectral_length; i++)
			out.array_data_spectral_space[i] *= scale;

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;

		return out;
	}

	DataArray<D> aliasing_scaleDown(
			std::size_t *i_new_resolution
	)
	{
		aliasing_scaled = true;

		DataArray<D> out(i_new_resolution);
//		out.temporary_data = true;
		out.aliasing_scaled = false;

		fftAliasingTestAndInit(out);

		requestDataInSpectralSpace();

		if (D != 2)
		{
			std::cerr << "TODO: Only 2D so far supported" << std::endl;
			exit(-1);
		}

		#pragma omp parallel for
		for (std::size_t j = 0; j < out.resolution_spec[1]/2; j++)
		{
			// lower quadrant
			memcpy(
					out.array_data_spectral_space+(j*out.range_spec_size[0])*2,
					array_data_spectral_space+(j*range_spec_size[0])*2,
					sizeof(double)*out.range_spec_size[0]*2
				);

			// top quadrant
			memcpy(
					out.array_data_spectral_space+((out.resolution_spec[1]-out.resolution_spec[1]/2+j)*out.range_spec_size[0])*2,
					array_data_spectral_space+((resolution_spec[1]-out.resolution_spec[1]/2+j)*range_spec_size[0])*2,
					sizeof(double)*out.range_spec_size[0]*2
				);
		}

		double scale = ((double)i_new_resolution[0]*(double)i_new_resolution[1])/((double)resolution[0]*(double)resolution[1]);

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < out.array_data_spectral_length; i++)
			out.array_data_spectral_space[i] *= scale;

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;

		return out;
	}

#endif


	/**
	 * Compute element-wise addition
	 */
	inline
	DataArray<D> operator+(
			const DataArray<D> &i_array_data
	)	const
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] =
					array_data_spectral_space[i]+
					i_array_data.array_data_spectral_space[i];

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = array_data_cartesian_space[i] + i_array_data.array_data_cartesian_space[i];

#endif

		return out;
	}




	/**
	 * Compute element-wise addition
	 */
	inline
	DataArray<D> operator+(
			const double i_value
	)	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] = array_data_spectral_space[i];

		double scale = resolution[0]*resolution[1];
		out.array_data_spectral_space[0] += i_value*scale;

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = array_data_cartesian_space[i]+i_value;

#endif
		return out;
	}




	/**
	 * Compute element-wise addition
	 */
	inline
	DataArray<D>& operator+=(
			const DataArray<D> &i_array_data
	)
	{
#if SWEET_USE_SPECTRAL_SPACE
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] +=
					i_array_data.array_data_spectral_space[i];

		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] += i_array_data.array_data_cartesian_space[i];
#endif

		return *this;
	}



	/**
	 * Compute element-wise addition
	 */
	inline
	DataArray<D>& operator+=(
			const double i_value
	)
	{
#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		double scale = resolution[0]*resolution[1];
		array_data_spectral_space[0] += i_value*scale;

		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] += i_value;

#endif
		return *this;
	}



	/**
	 * Compute multiplication with scalar
	 */
	inline
	DataArray<D>& operator*=(
			const double i_value
	)
	{
#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] *= i_value;

		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] *= i_value;

#endif
		return *this;
	}


	/**
	 * Compute division with scalar
	 */
	inline
	DataArray<D>& operator/=(
			const double i_value
	)
	{
#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] /= i_value;


		array_data_spectral_space_valid = true;
		array_data_cartesian_space_valid = false;

#else

		requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] /= i_value;

#endif
		return *this;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	DataArray<D>& operator-=(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

#if SWEET_USE_SPECTRAL_SPACE
		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] -=
					i_array_data.array_data_spectral_space[i];

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;
#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] -= i_array_data.array_data_cartesian_space[i];

#endif

		return *this;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	DataArray<D> operator-(
			const DataArray<D> &i_array_data
	)	const
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		DataArray<D> out(resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] =
					array_data_spectral_space[i]-
					i_array_data.array_data_spectral_space[i];

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]-
					i_array_data.array_data_cartesian_space[i];

#endif

		return out;
	}



	/**
	 * Compute element-wise subtraction
	 */
	inline
	DataArray<D> operator-(
			const double i_value
	)	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] = array_data_spectral_space[i];

		double scale = resolution[0]*resolution[1];
		out.array_data_spectral_space[0] -= i_value*scale;

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]-i_value;

#endif
		return out;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	DataArray<D> valueMinusThis(
			const double i_value
	)	const
	{
		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			out.array_data_spectral_space[i] = array_data_spectral_space[i];

		double scale = resolution[0]*resolution[1];
		out.array_data_spectral_space[0] = i_value*scale - out.array_data_spectral_space[0];

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] = i_value - array_data_cartesian_space[i];

#endif
		return out;
	}

	/**
	 * Invert sign
	 */
	inline
	DataArray<D>& operator-()
	{
#if SWEET_USE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_spectral_length; i++)
			array_data_spectral_space[i] = -array_data_spectral_space[i];

		array_data_cartesian_space_valid = false;
		array_data_spectral_space_valid = true;

#else

		requestDataInCartesianSpace();

		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] = -array_data_cartesian_space[i];

#endif

		return *this;
	}



	/**
	 * Compute element-wise multiplication
	 */
	inline
	DataArray<D> operator*(
			const DataArray<D> &i_array_data
	)	const
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		DataArray<D> out(i_array_data.resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE && SWEET_USE_SPECTRAL_DEALIASING

		DataArray<D> u = aliasing_scaleUp();
		DataArray<D> v = rw_array_data.aliasing_scaleUp();

		u.requestDataInCartesianSpace();
		v.requestDataInCartesianSpace();

		DataArray<D> scaled_output(u.resolution);

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < scaled_output.array_data_cartesian_length; i++)
			scaled_output.array_data_cartesian_space[i] =
					u.array_data_cartesian_space[i]*
					v.array_data_cartesian_space[i];

		scaled_output.array_data_cartesian_space_valid = true;
		scaled_output.array_data_spectral_space_valid = false;

		out = scaled_output.aliasing_scaleDown(out.resolution);

		out.array_data_cartesian_space_valid = false;
		out.array_data_spectral_space_valid = true;

#else
		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]*
					i_array_data.array_data_cartesian_space[i];

#if SWEET_USE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
#endif

#endif
		return out;
	}



	/**
	 * Compute multiplication with a scalar
	 */
	inline
	DataArray<D> operator*(
			const double i_value
	)	const
	{
		DataArray<D> out(resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE

		if (array_data_spectral_space_valid)
		{
			#pragma omp parallel for OPENMP_SIMD
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				out.array_data_spectral_space[i] =
						array_data_spectral_space[i]*i_value;

			out.array_data_cartesian_space_valid = false;
			out.array_data_spectral_space_valid = true;
		}
		else
		{
			assert(array_data_cartesian_space_valid);

			#pragma omp parallel for OPENMP_SIMD
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				out.array_data_cartesian_space[i] =
						array_data_cartesian_space[i]*i_value;

			out.array_data_cartesian_space_valid = true;
			out.array_data_spectral_space_valid = false;
		}

#else
		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]*i_value;
#endif
		return out;
	}


	/**
	 * Compute element-wise division
	 */
	inline
	DataArray<D> operator/(
			const DataArray<D> &i_array_data
	)	const
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		DataArray<D> out(this->resolution);
		out.temporary_data = true;

#if SWEET_USE_SPECTRAL_SPACE && SWEET_USE_SPECTRAL_DEALIASING

		DataArray<D> u = aliasing_scaleUp();
		DataArray<D> v = rw_array_data.aliasing_scaleUp();

		u.requestDataInCartesianSpace();
		v.requestDataInCartesianSpace();

		DataArray<D> scaled_output(u.resolution);

		#pragma omp parallel for OPENMP_SIMD
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

#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]/
					i_array_data.array_data_cartesian_space[i];

#if SWEET_USE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
#endif

#endif
		return out;
	}


	friend
	inline
	std::ostream& operator<<(std::ostream &o_ostream, const DataArray<D> &i_dataArray)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_dataArray;

		rw_array_data.requestDataInCartesianSpace();

		assert(D == 2);
		if (D == 2)
		{
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
		}
		return o_ostream;
	}


	inline
	void printSpectrum()
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)*this;

		rw_array_data.requestDataInSpectralSpace();

		assert(D == 2);
		if (D == 2)
		{
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
DataArray<2> operator*(
		const double i_value,
		const DataArray<2> &i_array_data
)
{
	return ((DataArray<2>&)i_array_data)*i_value;
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
DataArray<2> operator-(
		const double i_value,
		const DataArray<2> &i_array_data
)
{
	return ((DataArray<2>&)i_array_data).valueMinusThis(i_value);
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
DataArray<2> operator+(
		const double i_value,
		const DataArray<2> &i_array_data
)
{
	return ((DataArray<2>&)i_array_data)+i_value;
}

#endif /* SRC_DATAARRAY_HPP_ */
