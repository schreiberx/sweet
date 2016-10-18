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
#include <sweet/plane/PlaneDataConfig.hpp>

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
	PlaneDataConfig *planeDataConfig;


	/**
	 * local data in cartesian space
	 */
	double *physical_space_data;

#if SWEET_USE_LIBFFT || SWEET_USE_PLANE_SPECTRAL_SPACE
	bool physical_space_data_valid;

	std::complex<double> *spectral_space_data;
	bool spectral_space_data_valid;
#endif


#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	int kernel_size = -1;
	double *kernel_data = nullptr;
	int kernel_id = -1;
#endif



	/**
	 * prohibit empty initialization by making this method private
	 */
private:
	PlaneData()	:
		planeDataConfig(nullptr)
	{
	}



private:
	void p_allocate_buffers()
	{
		physical_space_data = MemBlockAlloc::alloc<double>(
				planeDataConfig->physical_array_data_number_of_elements*sizeof(double)
		);

//		std::cout << planeDataConfig->physical_array_data_number_of_elements << ", " << planeDataConfig->spectral_array_data_number_of_elements << std::endl;

		spectral_space_data = MemBlockAlloc::alloc< std::complex<double> >(
				planeDataConfig->spectral_array_data_number_of_elements*sizeof(std::complex<double>)
		);
	}



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
	)
	{
		planeDataConfig = i_dataArray.planeDataConfig;

		p_allocate_buffers();


#if SWEET_USE_PLANE_SPECTRAL_SPACE
		physical_space_data_valid = i_dataArray.physical_space_data_valid;
		if (physical_space_data_valid)
#endif
		{
			// use parallel copy for first touch policy!
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
				physical_space_data[i] = i_dataArray.physical_space_data[i];
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		spectral_space_data_valid = i_dataArray.spectral_space_data_valid;

		if (spectral_space_data_valid)
		{
			// use parallel copy for first touch policy!
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
				spectral_space_data[i] = i_dataArray.spectral_space_data[i];
		}
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
			PlaneDataConfig *i_planeDataConfig
	)
	{
		planeDataConfig = i_planeDataConfig;

		p_allocate_buffers();
	}



	/**
	 * default constructor
	 */
public:
	PlaneData(
		PlaneDataConfig *i_planeDataConfig
	)
#if SWEET_DEBUG
		:
		planeDataConfig(nullptr)
#endif
	{
		setup(i_planeDataConfig);
	}



public:
	~PlaneData()
	{
		MemBlockAlloc::free(physical_space_data, planeDataConfig->physical_array_data_number_of_elements*sizeof(double));

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		MemBlockAlloc::free(spectral_space_data, planeDataConfig->spectral_array_data_number_of_elements*sizeof(std::complex<double>));
#else
		MemBlockAlloc::free(kernel_data, sizeof(double)*kernel_size*kernel_size);
#endif
	}



	inline
	void physical_set(
			std::size_t j,
			std::size_t i,
			double i_value
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		physical_space_data_valid = true;
		spectral_space_data_valid = false;
#endif

		physical_space_data[j*planeDataConfig->physical_data_size[0]+i] = i_value;
	}



#if SWEET_USE_PLANE_SPECTRAL_SPACE
	inline
	void spectral_set(
			std::size_t j,
			std::size_t i,
			std::complex<double> &i_value
	)
	{
//		requestDataInSpectralSpace();

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		std::size_t idx = (j*planeDataConfig->spectral_data_size[0])+i;
		spectral_space_data[idx] = i_value;
	}

	inline
	void spectral_set(
			std::size_t j,
			std::size_t i,
			double i_a,
			double i_b
	)
	{
//		requestDataInSpectralSpace();

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		std::size_t idx = (j*planeDataConfig->spectral_data_size[0])+i;
		spectral_space_data[idx].real(i_a);
		spectral_space_data[idx].imag(i_b);
	}
#endif



	inline
	double physical_get(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInCartesianSpace();

		return physical_space_data[j*planeDataConfig->spectral_data_size[0]+i];
	}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	inline
	const std::complex<double>& spectral_get(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInSpectralSpace();

		return spectral_space_data[j*planeDataConfig->spectral_data_size[0]+i];
	}
#endif


	inline
	void physical_set_all(
			double i_value
	)
	{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_value;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		physical_space_data_valid = true;
		spectral_space_data_valid = false;
#endif
	}



	/**
	 * Set the values in the specified row
	 */
	inline
	void physical_set_row(
			int i_row,
			double i_value
	)
	{
		if (i_row < 0)
			i_row += planeDataConfig->physical_data_size[1];

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_data_size[0]; i++)
			physical_space_data[i_row*planeDataConfig->physical_data_size[0]+i] = i_value;
	}

	/**
	 * Copy values from another row and flip sign
	 */
	inline
	void physical_copy_row_inv_sign(
			int i_src_row,
			int i_dst_row
	)
	{
		if (i_src_row < 0)
			i_src_row += planeDataConfig->physical_data_size[1];

		if (i_dst_row < 0)
			i_dst_row += planeDataConfig->physical_data_size[1];

		requestDataInCartesianSpace();

		std::size_t src_idx = i_src_row*planeDataConfig->physical_data_size[0];
		std::size_t dst_idx = i_dst_row*planeDataConfig->physical_data_size[0];

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_data_size[0]; i++)
			physical_space_data[dst_idx+i] = -physical_space_data[src_idx+i];
	}


	/**
	 * Copy values from another row
	 */
	inline
	void physical_copy_row(
			int i_src_row,
			int i_dst_row
	)
	{
		if (i_src_row < 0)
			i_src_row = planeDataConfig->physical_data_size[1]+i_src_row;

		if (i_dst_row < 0)
			i_dst_row = planeDataConfig->physical_data_size[1]+i_dst_row;

		requestDataInCartesianSpace();

		std::size_t src_idx = i_src_row*planeDataConfig->physical_data_size[0];
		std::size_t dst_idx = i_dst_row*planeDataConfig->physical_data_size[0];

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_data_size[0]; i++)
			physical_space_data[dst_idx+i] = physical_space_data[src_idx+i];
	}






#if SWEET_USE_PLANE_SPECTRAL_SPACE==1

	inline
	double spec_getRe(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInSpectralSpace();

		return spectral_space_data[j*planeDataConfig->spectral_data_size[0]+i].real();
	}



	inline
	double spec_getIm(
			std::size_t j,
			std::size_t i
	)	const
	{
		requestDataInSpectralSpace();

		return spectral_space_data[j*planeDataConfig->spectral_data_size[0]+i].imag();
	}


	inline
	void set_spec(
			std::size_t j,
			std::size_t i,
			double i_value_re,
			double i_value_im
	)
	{
		std::size_t idx =	j*planeDataConfig->spectral_data_size[0]+i;

		spectral_space_data[idx].real(i_value_re);
		spectral_space_data[idx].imag(i_value_im);

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
		{
			spectral_space_data[i].real(i_value_re);
			spectral_space_data[i].imag(i_value_im);
		}

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
	}


#if 0
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
			std::size_t idx =	j*planeDataConfig->spectral_data_size[0]+
								(i);

			spectral_space_data[idx*2+0] = i_value_re;
			spectral_space_data[idx*2+1] = i_value_im;
		}

		if (j != 0)
		{	// upper part of y
			std::size_t idx =	((resolution[1]-j))*spectral_data_size[0]+
								(i);

			spectral_space_data[idx*2+0] = i_value_re;
			// IMPORTANT! the imaginary component is mirrored!!!
			spectral_space_data[idx*2+1] = i_value_im;
		}

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

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
		double c = cos(2.0*M_PIl*i_value_phase_shift)*i_value_amplitude;
		double s = sin(2.0*M_PIl*i_value_phase_shift)*i_value_amplitude;

		set_spec_spectrum(j, i, c, s);
	}
#endif


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
		assert(i >= 0 && i < planeDataConfig->spectral_data_size[0]/2);
		assert(j >= 0 && j < planeDataConfig->spectral_data_size[1]/2);

		{	// lower part of y
			std::size_t idx = j*planeDataConfig->spectral_data_size[0]+i;

			spectral_space_data[idx*2+0] = i_value_re;
			spectral_space_data[idx*2+1] = i_value_im;
		}

		if (j != 0)
		{	// upper part of y
			std::size_t idx =	(planeDataConfig->spectral_data_size[1]-j)*planeDataConfig->spectral_data_size[0]+i;

			spectral_space_data[idx*2+0] = i_value_re;
			// IMPORTANT! the imaginary component is mirrored!!!
			spectral_space_data[idx*2+1] = -i_value_im;
		}

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
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
		assert(i >= 0 && i < planeDataConfig->spectral_data_size[0]/2);
		assert(j >= 0 && j < planeDataConfig->spectral_data_size[1]/2);

		{	// lower part of y
			std::size_t idx = j*planeDataConfig->spectral_data_size[0]+i;

			spectral_space_data[idx*2+0] = i_value_re;
			spectral_space_data[idx*2+1] = i_value_im;
		}

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
	}

#if 0
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
			std::size_t idx =	(resolution[1]-j)*spectral_data_size[0]+i;

			spectral_space_data[idx*2+0] = i_value_re;
			// IMPORTANT! the imaginary component is mirrored!!!
			spectral_space_data[idx*2+1] = -i_value_im;
		}

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		checkConsistency();
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

		if (spectral_space_data_valid)
			return *this;		// nothing to do

		if (!physical_space_data_valid)
		{
			std::cerr << "Spectral data not available! Is this maybe a non-initialized operator?" << std::endl;
			assert(false);
			exit(1);
		}

		PlaneData *rw_array_data = (PlaneData*)this;

		planeDataConfig->fft_physical_to_spectral(rw_array_data->physical_space_data, rw_array_data->spectral_space_data);

		rw_array_data->spectral_space_data_valid = true;
		rw_array_data->physical_space_data_valid = false;

#endif
		return *this;
	}


	inline
	const PlaneData& requestDataInCartesianSpace() const
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE==1

		if (physical_space_data_valid)
			return *this;		// nothing to do

		assert(spectral_space_data_valid == true);

		PlaneData *rw_array_data = (PlaneData*)this;

		planeDataConfig->fft_spectral_to_physical(rw_array_data->spectral_space_data, rw_array_data->physical_space_data);

		rw_array_data->spectral_space_data_valid = false;
		rw_array_data->physical_space_data_valid = true;
#endif
	return *this;
	}

#endif


	inline
	PlaneData return_one_if_positive()
	{
		PlaneData out(planeDataConfig);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = (physical_space_data[i] > 0 ? 1 : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
#endif
		return out;
	}



	inline
	PlaneData return_value_if_positive()	const
	{
		PlaneData out(planeDataConfig);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = (physical_space_data[i] > 0 ? physical_space_data[i] : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
#endif
		return out;
	}


	inline
	PlaneData return_one_if_negative()	const
	{
		PlaneData out(planeDataConfig);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = (physical_space_data[i] < 0 ? 1 : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
#endif
		return out;
	}


	inline
	PlaneData return_value_if_negative()	const
	{
		PlaneData out(planeDataConfig);

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = (physical_space_data[i] < 0 ? physical_space_data[i] : 0);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
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
#if SWEET_THREADING
#pragma omp parallel for reduction(&&:isallfinite)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			isallfinite = isallfinite && std::isfinite(physical_space_data[i]);
//			isallfinite = isallfinite && (physical_space[i]<1000);


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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			maxabs = std::max(maxabs, std::abs(physical_space_data[i]));

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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			sum += physical_space_data[i]*physical_space_data[i];

		sum = std::sqrt(sum/(double)(planeDataConfig->physical_array_data_number_of_elements));

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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			double value = physical_space_data[i]*physical_space_data[i];

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		sum = std::sqrt(sum/(double)(planeDataConfig->physical_array_data_number_of_elements));

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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			maxvalue = std::max(maxvalue, physical_space_data[i]);

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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			minvalue = std::min(minvalue, physical_space_data[i]);

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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			sum += physical_space_data[i];

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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			double value = physical_space_data[i];

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
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			sum += std::abs(physical_space_data[i]);


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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{

			double value = std::abs(physical_space_data[i]);
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
	 * return the sqrt of the sum of the squared values
	 */
	double reduce_norm2()	const
	{
		requestDataInCartesianSpace();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			sum += physical_space_data[i]*physical_space_data[i];


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
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			double value = physical_space_data[i]*physical_space_data[i];

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
		{
			double re = spectral_space_data[i].real();
			double im = spectral_space_data[i].imag();
			maxabs = std::max(maxabs, std::sqrt(re*re + im*im));
		}

		return maxabs;
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

		physical_set_all(0);

		// left lower corner
		//kernel_cart[    0:R[0]+1,       0:R[0]+1    ] = conv_kernel[    R[0]:,  R[0]:   ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = ky-R;
				std::size_t cx = kx-R;

				if (0 > cx || planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || planeDataConfig->physical_data_size[1] <= cy)
					continue;

				physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// right bottom corner
		//kernel_cart[    0:R[0]+1,       res-R[0]:   ] = conv_kernel[    R[0]:,  0:R[0]  ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = ky-R;
				std::size_t cx = planeDataConfig->physical_data_size[0] - R + kx;

				if (0 > cx || planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || planeDataConfig->physical_data_size[1] <= cy)
					continue;

				physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// left top corner
		//kernel_cart[    res-R[0]:,      0:R[0]+1   ] = conv_kernel[    0:R[0],  R[0]:  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = planeDataConfig->physical_data_size[1] - R + ky;
				std::size_t cx = kx-R;

				if (0 > cx || planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || planeDataConfig->physical_data_size[1] <= cy)
					continue;

				physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// right top corner
		//kernel_cart[    res-R[0]:,      res-R[0]:   ] = conv_kernel[    0:R[0], 0:R[0]  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = planeDataConfig->physical_data_size[1] - R + ky;
				std::size_t cx = planeDataConfig->physical_data_size[0] - R + kx;

				if (0 > cx || planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || planeDataConfig->physical_data_size[1] <= cy)
					continue;

				physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}

		physical_space_data_valid = true;
		spectral_space_data_valid = false;

		requestDataInSpectralSpace();	/// convert kernel_data; to spectral space
#endif
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
		if (this->spectral_space_data_valid)
		{
			this->requestDataInCartesianSpace();
			was_spectral = true;
			this->spectral_space_data_valid=false;
		}else if(!this->physical_space_data_valid){
			std::cout << "Uninitialized PlaneData in op_stencil_Re_3x3" << std::endl;
			exit(-1);
		}
		out.spectral_space_data_valid=false;
#endif


		int res_x = planeDataConfig->physical_data_size[0];
		int res_y = planeDataConfig->physical_data_size[1];


#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y, out, i_kernel_data)
#endif
		for (int y = 0; y < res_y; y++)
		{
			for (int x = 0; x < res_x; x++)
			{
				double *data_out = &out.physical_space_data[(y*res_x+x)];
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

						double dre = physical_space_data[(pos_y*res_x+pos_x)];

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
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		bool was_spectral = false;
		if (this->spectral_space_data_valid)
		{
			this->requestDataInCartesianSpace();
			was_spectral = true;
			this->spectral_space_data_valid=false;
		}else if(!this->physical_space_data_valid){
			std::cout << "Uninitialized PlaneData in op_stencil_Re_3x3" << std::endl;
			exit(-1);
		}
		out.spectral_space_data_valid=false;
#endif

		int res_x = planeDataConfig->physical_data_size[0];
		int res_y = planeDataConfig->physical_data_size[1];


#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y, out, i_kernel_data)
#endif
		for (int y = 0; y < res_y; y++)
		{
			for (int x = 0; x < res_x; x++)
			{
				double *data_out = &out.physical_space_data[(y*res_x+x)];
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

						double dre = physical_space_data[(pos_y*res_x+pos_x)];

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
		PlaneData out(planeDataConfig);

		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		// only makes sense, if this is an operator created in spectral space
		assert(i_array_data.spectral_space_data_valid == true);

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

		// determine maximum value for tolerance
		double max_value = i_array_data.reduce_spec_maxAbs();
		i_tolerance *= max_value;
		i_tolerance *= (planeDataConfig->physical_data_size[0]+planeDataConfig->physical_data_size[1]);	// the larger the matrix, the less the accuracy

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
		{
			double ar = spectral_space_data[i].real();
			double ai = spectral_space_data[i].imag();
			double br = i_array_data.spectral_space_data[i].real();
			double bi = i_array_data.spectral_space_data[i].imag();

			double den = (br*br+bi*bi);

			if (std::abs(den) <= i_tolerance)
			{
				// For inverting differential operators, this is the integration constant C
				out.spectral_space_data[i] = ar*i_denom_zeros_scalar;
				out.spectral_space_data[i+1] = ai*i_denom_zeros_scalar;
			}
			else
			{
				out.spectral_space_data[i].real((ar*br + ai*bi)/den);
				out.spectral_space_data[i].imag((ai*br - ar*bi)/den);
			}
		}

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

		return out;
	}
#endif



#if SWEET_USE_PLANE_SPECTRAL_SPACE
#if 0
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
		for (std::size_t y = planeDataConfig->spectral_data_size[1]-1; y > planeDataConfig->spectral_data_size[1]/2; y--)
		{
			i=0;
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{

				//double value_re = spec_getRe(y, x);
				//double value_im = spec_getIm(y, x);
				if( x > 2*(resolution_spec[0]-1)/3 || j > 2*(resolution_spec[1]/2)/3 )
				//if( x > 2*(resolution_spec[0]-1)/3-1 || j > 2*(resolution_spec[1]/2)/3-1 )
				{
					spectral_set( y, x, 0.0, 0.0);
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
		j=planeDataConfig->spectral_data_size[1]/2; //modenumber in y
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
					spectral_set( y, x, 0.0, 0.0);
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

		return *this;
	}
#endif

#endif


public:
	/**
	 * assignment operator
	 */
	PlaneData &operator=(double i_value)
	{
		physical_set_all(i_value);

		return *this;
	}


public:
	/**
	 * assignment operator
	 */
	PlaneData &operator=(int i_value)
	{
		physical_set_all(i_value);

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
		planeDataConfig = i_dataArray.planeDataConfig;

		planeDataConfig->physical_array_data_number_of_elements = i_dataArray.planeDataConfig->physical_array_data_number_of_elements;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (i_dataArray.physical_space_data_valid)
		{
			physical_space_data_valid = true;
#endif
			/*
			 * If this data was generated based on temporary data sets (e.g. via h = hu/u), then only swap pointers.
			 */
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
				physical_space_data[i] = i_dataArray.physical_space_data[i];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		}
		else
		{
			physical_space_data_valid = false;
		}

		planeDataConfig->spectral_array_data_number_of_elements = i_dataArray.planeDataConfig->spectral_array_data_number_of_elements;

		if (i_dataArray.spectral_space_data_valid)
		{
			spectral_space_data_valid = true;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
				spectral_space_data[i] = i_dataArray.spectral_space_data[i];
		}
		else
		{
			spectral_space_data_valid = false;
		}
#endif

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
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			out.spectral_space_data[i] = spectral_space_data[i]*i_array_data.spectral_space_data[i];

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

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
						double &data_out = out.physical_space[y*res_x+x];
						data_out = 0;

						int pos_y = y;
						assert(pos_y >= 0 && pos_y < res_y);

						if (x > 0 && x < res_x-1)
						{
							double *kernel_scalar_ptr = &kernel_data[3];
							double *data_scalar_ptr = &i_array_data.physical_space_data[pos_y*res_x+x-1];

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
								double data_scalar = i_array_data.physical_space_data[pos_y*res_x+pos_x];

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
							double &data_out = out.physical_space_data[y*res_x+x];
							data_out = 0;

							int pos_y = y+res_y;
							pos_y -= (pos_y >= res_y ? res_y : 0);
							pos_y += (pos_y < 0 ? res_y : 0);

							assert(pos_y >= 0 && pos_y < res_y);

							if (x > 0 && x < res_x-1)
							{
								double *kernel_scalar_ptr = &kernel_data[3];
								double *data_scalar_ptr = &i_array_data.physical_space_data[pos_y*res_x+x-1];

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
									double data_scalar = i_array_data.physical_space_data[pos_y*res_x+pos_x];

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
							double &data_out = out.physical_space_data[y*res_x+x];
							data_out = 0;

							if (y > 0 && y < res_y-1)
							{
								double *kernel_scalar_ptr = &kernel_data[1];
								double *data_scalar_ptr = &i_array_data.physical_space_data[(y-1)*res_x+x];

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
									double data_scalar = i_array_data.physical_space_data[pos_y*res_x+x];

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
							double &data_out = out.physical_space_data[y*res_x+x];
							data_out = 0;

							if (y > 0 && y < res_y-1)
							{
								double *kernel_scalar_ptr = &kernel_data[1];
								double *data_scalar_ptr = &i_array_data.physical_space_data[(y-1)*res_x+x];

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
									double data_scalar = i_array_data.physical_space_data[pos_y*res_x+x];

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
						double &data_out = out.physical_space_data[y*res_x+x];
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
								double data_scalar = i_array_data.physical_space_data[pos_y*res_x+pos_x];

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

		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < planeDataConfig->spectral_data_size[1]; j++)
		{
			std::size_t idx = j*planeDataConfig->spectral_data_size[0];

			std::complex<double> *src1 = &spectral_space_data[idx];
			std::complex<double> *src2 = &i_array_data.spectral_space_data[idx];

			std::complex<double> *dst = &out.spectral_space_data[idx];

			for (std::size_t i = 0; i < planeDataConfig->spectral_data_size[0]; i++)
				dst[i] = src1[i] + src2[i];
		}

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = array_data_cartesian_space[i] + i_array_data.physical_space[i];

#endif


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
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		{
			requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
				out.spectral_space_data[i] = spectral_space_data[i];

			double scale = planeDataConfig->spectral_data_size[0]*planeDataConfig->spectral_data_size[1];
			out.spectral_space_data[0] += i_value*scale;

			out.physical_space_data_valid = false;
			out.spectral_space_data_valid = true;
		}

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.array_data_cartesian_space[i] = physical_space_data[i]+i_value;

#endif

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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] +=
					i_array_data.spectral_space_data[i];

		spectral_space_data_valid = true;
		physical_space_data_valid = false;

#else

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			array_data_cartesian_space[i] += i_array_data.physical_space_data[i];
#endif

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

		double scale = planeDataConfig->spectral_array_data_number_of_elements;
		spectral_space_data[0] += i_value*scale;

		spectral_space_data_valid = true;
		physical_space_data_valid = false;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] += i_value;

#endif

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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] *= i_value;

		spectral_space_data_valid = true;
		physical_space_data_valid = false;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] *= i_value;

#endif

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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] /= i_value;


		spectral_space_data_valid = true;
		physical_space_data_valid = false;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] /= i_value;

#endif

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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] -=
					i_array_data.spectral_space_data[i];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			array_data_cartesian_space[i] -= i_array_data.physical_space[i];

#endif

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

		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			out.spectral_space_data[i] =
					spectral_space_data[i]-
					i_array_data.spectral_space_data[i];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] =
					array_data_cartesian_space[i]-
					i_array_data.physical_space_data[i];

#endif

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
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

//#pragma error "TOOD: make this depending on rexi_par_sum!"
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			out.spectral_space_data[i] = spectral_space_data[i];

		double scale = planeDataConfig->spectral_data_size[0]*planeDataConfig->spectral_data_size[1];
		out.spectral_space_data[0] -= i_value*scale;

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.array_data_cartesian_space[i] =
					physical_space_data[i]-i_value;

#endif

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
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			out.spectral_space_data[i] = -spectral_space_data[i];

		double scale = planeDataConfig->physical_data_size[0]*planeDataConfig->physical_data_size[1];
		out.spectral_space_data[0] = i_value*scale + out.spectral_space_data[0];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else

		requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_value - physical_space[i];

#endif

		return out;
	}

	/**
	 * Invert sign
	 */
	inline
	PlaneData operator-()
	{
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		requestDataInSpectralSpace();

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			out.spectral_space_data[i] = -spectral_space_data[i];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else

		requestDataInCartesianSpace();

		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.array_data_cartesian_space[i] = -physical_space_data[i];

#endif

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
		PlaneData out(planeDataConfig);

		requestDataInCartesianSpace();
		i_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] =
					physical_space_data[i]*
					i_array_data.physical_space_data[i];

	#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;
	#endif

		return out;
	}

	/**
	 * Compute element-wise multiplication
	 */
	inline
	PlaneData operator/(
			const PlaneData &i_array_data	///< this class times i_array_data
	)	const
	{
		PlaneData out(planeDataConfig);

		requestDataInCartesianSpace();
		i_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]/i_array_data.physical_space_data[i];

	#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;
	#endif

		return out;
	}


#if 0
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

		PlaneData out(planeDataConfig);

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
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;
#endif

		return out;
	}
#endif


	/**
	 * Compute multiplication with a scalar
	 */
	inline
	PlaneData operator*(
			const double i_value
	)	const
	{
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		if (spectral_space_data_valid)
		{
#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
				out.spectral_space_data[i] =
						spectral_space_data[i]*i_value;

			out.physical_space_data_valid = false;
			out.spectral_space_data_valid = true;
		}
		else
		{
			assert(physical_space_data_valid);

#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
				out.physical_space_data[i] =
						physical_space_data[i]*i_value;

			out.physical_space_data_valid = true;
			out.spectral_space_data_valid = false;
		}

#else
#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.array_data_cartesian_space[i] =
					physical_space_data[i]*i_value;
#endif

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
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (physical_space_data_valid)
		{
#endif
#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
				out.physical_space_data[i] = physical_space_data[i] / i_value;
#if SWEET_USE_PLANE_SPECTRAL_SPACE
			out.physical_space_data_valid = true;
			out.spectral_space_data_valid = false;
		}
		else
		{
			assert(spectral_space_data_valid);

#if SWEET_THREADING
			#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
				out.spectral_space_data[i] = spectral_space_data[i] / i_value;

			out.physical_space_data_valid = false;
			out.spectral_space_data_valid = true;
		}
#endif

		return out;
	}


#if 0
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

		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE && SWEET_USE_PLANE_SPECTRAL_DEALIASING

		PlaneData u = aliasing_scaleUp();
		PlaneData v = rw_array_data.aliasing_scaleUp();

		u.requestDataInCartesianSpace();
		v.requestDataInCartesianSpace();

		PlaneData scaled_output(u.resolution, true);

#if SWEET_THREADING
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < scaled_output.planeDataConfig->physical_array_data_number_of_elements; i++)
			scaled_output.physical_space_data[i] =
					u.physical_space_data[i]/
					v.physical_space_data[i];

		scaled_output.physical_space_data_valid = true;
		scaled_output.spectral_space_data_valid = false;

		out = scaled_output.aliasing_scaleDown(out.resolution);

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else
		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] =
					physical_space_data[i]/
					i_array_data.physical_space_data[i];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;
#endif

#endif

		out.checkConsistency();
		return out;
	}
#endif


	friend
	inline
	std::ostream& operator<<(
			std::ostream &o_ostream,
			const PlaneData &i_dataArray
	)
	{
		PlaneData &rw_array_data = (PlaneData&)i_dataArray;

		rw_array_data.requestDataInCartesianSpace();

		for (int y = rw_array_data.planeDataConfig->physical_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < rw_array_data.planeDataConfig->physical_data_size[0]; x++)
			{
				double value = rw_array_data.physical_get(y, x);
//					if (std::abs(value) < 1e-13)
//						value = 0;
				std::cout << value << "\t";
			}
			std::cout << std::endl;
		}

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
		PlaneData out(planeDataConfig);

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			out.spectral_space_data[i] = spectral_space_data[i] + i_value;

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}
#endif

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Return Plane Array with all spectral coefficients a+bi --> 1/(a+bi)
	 */
	inline
	PlaneData spec_invert()	const
	{
		PlaneData out(planeDataConfig);

		requestDataInSpectralSpace();

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			out.spectral_space_data[i] = 1.0/spectral_space_data[i];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}


	inline
	void printSpectrum()	const
	{
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.requestDataInSpectralSpace();

		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
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
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.requestDataInSpectralSpace();

		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
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
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.requestDataInSpectralSpace();

		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
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
		requestDataInCartesianSpace();

		std::ostream &o_ostream = std::cout;

		o_ostream << std::setprecision(i_precision);

		for (int y = planeDataConfig->physical_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				o_ostream << physical_get(y, x);

				if (x < planeDataConfig->physical_data_size[0]-1)
					o_ostream << '\t';
				else
					o_ostream << std::endl;
			}
		}

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
		requestDataInCartesianSpace();

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		for (int y = planeDataConfig->physical_res[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_res[0]; x++)
			{
				file << physical_get(y, x);

				if (x < planeDataConfig->physical_res[0]-1)
					file << i_separator;
				else
					file << std::endl;
			}
		}

		return true;
	}


#if 0
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
		requestDataInCartesianSpace();

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		file << "# vtk DataFile Version 2.0" << std::endl;
		file << "Rectangular solid example" << std::endl;
		file << "ASCII" << std::endl;
		file << "DATASET RECTILINEAR_GRID" << std::endl;
		file << "DIMENSIONS " << planeDataConfig->physical_res[0]+1 << " " << planeDataConfig->physical_res[1]+1 << " 1" << std::endl;

		file << "X_COORDINATES " << planeDataConfig->physical_res[0]+1 << " float" << std::endl;
		for (std::size_t x = 0; x < planeDataConfig->physical_res[0]+1; x++)
			file << (double)x/((double)planeDataConfig->physical_res[0]+1) << std::endl;

		file << "Y_COORDINATES " << planeDataConfig->physical_res[1]+1 << " float" << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->physical_res[1]+1; y++)
			file << (double)y/((double)planeDataConfig->physical_res[1]+1) << std::endl;

		file << "Z_COORDINATES 1 float" << std::endl;
		file << "0" << std::endl;


		std::string title = i_title;
		std::replace(title.begin(), title.end(), ' ', '_');
		file << "CELL_DATA " << planeDataConfig->physical_res[0]*planeDataConfig->physical_res[1] << std::endl;
		file << "SCALARS " << title << " float 1" << std::endl;
		file << "LOOKUP_TABLE default" << std::endl;

		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			file << physical_space_data[i] << std::endl;

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
	bool file_loadData_physical(
			const char *i_filename,		///< Name of file to load data from
			bool i_binary_data = false	///< load as binary data (disabled per default)
	)
	{
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


			std::size_t expected_size = sizeof(double)*planeDataConfig->physical_res[0]*planeDataConfig->physical_res[1];

			if (size != expected_size)
			{
				std::cerr << "Error while loading data from file " << i_filename << ":" << std::endl;
				std::cerr << "Size of file " << size << " does not match expected size of " << expected_size << std::endl;
				exit(-1);
			}

			if (!file.read((char*)physical_space_data, expected_size))
			{
				std::cerr << "Error while loading data from file " << i_filename << std::endl;
				exit(1);
			}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			physical_space_data_valid = true;
			spectral_space_data_valid = false;
#endif
			return true;
		}
		std::ifstream file(i_filename);

		for (std::size_t row = 0; row < planeDataConfig->physical_res[1]; row++)
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

				physical_set(planeDataConfig->physical_res[1]-row-1, col, i_value);

				col++;
				last_pos = pos+1;
		    }

			if (col < planeDataConfig->physical_res[0])
			{
				std::cerr << "Failed to read data from file " << i_filename << " in line " << row << ", column " << col << std::endl;
				return false;
			}
		}

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
	return ((PlaneData&)i_array_data)+i_value;
}

#endif /* SRC_DATAARRAY_HPP_ */
