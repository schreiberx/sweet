/*
 * PlaneData.hpp
 *
 *  Created on: 28 Jun 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <schreiberx@gmail.com>
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
#include <sweet/FatalError.hpp>

#include <sweet/plane/PlaneDataConfig.hpp>
#include <sweet/plane/PlaneData_Kernels.hpp>


/*
 * Precompiler helper functions to handle loops in spectral and physical space
 */
#if SWEET_USE_PLANE_SPECTRAL_SPACE

	#if SWEET_THREADING
		#define PLANE_DATA_SPECTRAL_FOR_IDX(CORE)					\
			_Pragma("omp parallel for proc_bind(spread)")			\
			for (int r = 0; r < 2; r++)								\
			{														\
				_Pragma("omp parallel for OPENMP_PAR_SIMD proc_bind(close) collapse(2)")		\
				for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
				{				\
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
					{			\
						std::size_t idx = jj*planeDataConfig->spectral_data_size[0]+ii;	\
						CORE	\
					}			\
				}				\
			}

	#else

		#define PLANE_DATA_SPECTRAL_FOR_IDX(CORE)					\
			for (int r = 0; r < 2; r++)								\
			{														\
				for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
				{				\
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
					{			\
						std::size_t idx = jj*planeDataConfig->spectral_data_size[0]+ii;	\
						CORE	\
					}			\
				}				\
			}

	#endif
#endif


#if SWEET_THREADING
	#define PLANE_DATA_PHYSICAL_FOR_IDX(CORE)				\
		_Pragma("omp parallel for OPENMP_PAR_SIMD proc_bind(close)")	\
			for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)	\
			{	CORE;	}



	#define PLANE_DATA_PHYSICAL_FOR_IDX_REDUCTION(CORE, REDUCTION)				\
		_Pragma("omp parallel for OPENMP_PAR_SIMD "##REDUCTION##" proc_bind(close)")	\
			for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)	\
			{	CORE;	}

	#define PLANE_DATA_PHYSICAL_FOR_2D_IDX(CORE)										\
			_Pragma("omp parallel for OPENMP_PAR_SIMD proc_bind(close) collapse(2)")	\
			for (std::size_t j = 0; j < planeDataConfig->physical_data_size[1]; j++)						\
			{				\
				for (std::size_t i = 0; i < planeDataConfig->physical_data_size[0]; i++)	\
				{			\
					std::size_t idx = j*planeDataConfig->physical_data_size[0]+i;	\
					CORE;	\
				}			\
			}

#else

	#define PLANE_DATA_PHYSICAL_FOR_IDX(CORE)				\
			for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)	\
			{	CORE;	}

	#define PLANE_DATA_PHYSICAL_FOR_2D_IDX(CORE)										\
			for (std::size_t j = 0; j < planeDataConfig->physical_data_size[1]; j++)						\
			{				\
				for (std::size_t i = 0; i < planeDataConfig->physical_data_size[0]; i++)	\
				{			\
					std::size_t idx = j*planeDataConfig->physical_data_size[0]+i;	\
					CORE;	\
				}			\
			}

#endif



/*
 * this option activates if data has to be allocated for the spectral space
 */
#ifndef SWEET_USE_PLANE_SPECTRAL_SPACE
	#define SWEET_USE_PLANE_SPECTRAL_SPACE	1
#endif

/*
 * This option should !!! NEVER !!! show up here in this file
 */
#ifndef SWEET_USE_PLANE_SPECTRAL_DEALIASING
	#define SWEET_USE_PLANE_SPECTRAL_DEALIASING 1
#endif


#if SWEET_USE_LIBFFT
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
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
		:
		private PlaneData_Kernels
#endif
{
public:
	PlaneDataConfig *planeDataConfig;


	/**
	 * physical space data
	 */
	double *physical_space_data;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	bool physical_space_data_valid;

	/**
	 * spectral space data
	 */
	std::complex<double> *spectral_space_data;
	bool spectral_space_data_valid;
#endif


	/**
	 * prohibit empty initialization by making this method private
	 */
private:
	PlaneData()	:
		planeDataConfig(nullptr)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		,
		physical_space_data_valid(false),
		spectral_space_data_valid(false)
#endif
	{
	}

	/**
	 * dummy initialization by handing over an unused integer
	 */
public:
	PlaneData(int i)	:
		planeDataConfig(nullptr)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		,
		physical_space_data_valid(false),
		spectral_space_data_valid(false)
#endif
	{
	}



private:
	void p_allocate_buffers()
	{
		physical_space_data = MemBlockAlloc::alloc<double>(
				planeDataConfig->physical_array_data_number_of_elements*sizeof(double)
		);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		spectral_space_data = MemBlockAlloc::alloc< std::complex<double> >(
				planeDataConfig->spectral_array_data_number_of_elements*sizeof(std::complex<double>)
		);
#endif
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
		assert(i_dataArray.planeDataConfig != nullptr);

		planeDataConfig = i_dataArray.planeDataConfig;

		p_allocate_buffers();


#if SWEET_USE_PLANE_SPECTRAL_SPACE
		physical_space_data_valid = i_dataArray.physical_space_data_valid;
		if (physical_space_data_valid)
#endif
		{
			PLANE_DATA_PHYSICAL_FOR_IDX(
					physical_space_data[idx] = i_dataArray.physical_space_data[idx];
			);
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		spectral_space_data_valid = i_dataArray.spectral_space_data_valid;

		if (spectral_space_data_valid)
		{
			PLANE_DATA_SPECTRAL_FOR_IDX(
					spectral_space_data[idx] = i_dataArray.spectral_space_data[idx];
			);
		}
#endif
	}



public:
	/**
	 * setup the PlaneData in case that the special
	 * empty constructor with int as a parameter was used.
	 *
	 * Calling this setup function should be in general avoided
	 * except in very rare circumstances
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
	)	:
		planeDataConfig(nullptr)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		,physical_space_data_valid(false),
		spectral_space_data_valid(false)
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
#endif
	}


	inline
	void p_physical_set(
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


	void physical_update_lambda_array_indices(
			std::function<void(int,int,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (spectral_space_data_valid)
			request_data_physical();
#endif

		PLANE_DATA_PHYSICAL_FOR_2D_IDX(
				i_lambda(i, j, physical_space_data[idx])
		);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		physical_space_data_valid = true;
		spectral_space_data_valid = false;
#endif
	}


	inline
	double physical_get(
			std::size_t j,
			std::size_t i
	)	const
	{
		request_data_physical();

		return physical_space_data[j*planeDataConfig->physical_data_size[0]+i];
	}


	inline
	void physical_set_all(
			double i_value
	)
	{
		PLANE_DATA_PHYSICAL_FOR_IDX(
				physical_space_data[idx] = i_value;
		);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		physical_space_data_valid = true;
		spectral_space_data_valid = false;
#endif
	}



	inline
	void physical_set_zero()
	{
		PLANE_DATA_PHYSICAL_FOR_IDX(
				physical_space_data[idx] = 0;
		);

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

		request_data_physical();

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

		request_data_physical();

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

		request_data_physical();

		std::size_t src_idx = i_src_row*planeDataConfig->physical_data_size[0];
		std::size_t dst_idx = i_dst_row*planeDataConfig->physical_data_size[0];

#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_data_size[0]; i++)
			physical_space_data[dst_idx+i] = physical_space_data[src_idx+i];
	}




#if SWEET_USE_PLANE_SPECTRAL_SPACE

	void spectral_zeroAliasingModes()	const
	{
#if SWEET_THREADING
#pragma omp parallel for proc_bind(spread)
#endif
		for (int k = 0; k < 2; k++)
		{
			if (k == 0)
			{
				/*
				 * First process part between top and bottom spectral data blocks
				 */
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD proc_bind(close) collapse(2)
#endif
				for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[0][1][1]; jj < planeDataConfig->spectral_data_iteration_ranges[1][1][0]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii++)
					{
						spectral_space_data[jj*planeDataConfig->spectral_data_size[0]+ii] = 0;
					}

			}
			else
			{
				/*
				 * Then process the aliasing block on the right side
				 */
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD proc_bind(close) collapse(2)
#endif
				for (std::size_t jj = 0; jj < planeDataConfig->spectral_data_size[1]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii < planeDataConfig->spectral_data_size[0]; ii++)
					{
						spectral_space_data[jj*planeDataConfig->spectral_data_size[0]+ii] = 0;
					}
			}
		}

	}

	inline
	const std::complex<double>& spectral_get(
			std::size_t j,
			std::size_t i
	)	const
	{
		request_data_spectral();

		return spectral_space_data[j*planeDataConfig->spectral_data_size[0]+i];
	}



#if 1
	inline
	void p_spectral_set(
			std::size_t j,
			std::size_t i,
			const std::complex<double> &i_value
	)
	{
//		request_data_spectral();

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		std::size_t idx = (j*planeDataConfig->spectral_data_size[0])+i;
		spectral_space_data[idx] = i_value;
	}
#endif



	inline
	void spectral_set_all(
			double i_value_re,
			double i_value_im
	)
	{
		PLANE_DATA_SPECTRAL_FOR_IDX(
				spectral_space_data[idx].real(i_value_re);
				spectral_space_data[idx].imag(i_value_im);
		);

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
	}


	inline
	void spectral_set_zero()
	{
		PLANE_DATA_SPECTRAL_FOR_IDX(
				spectral_space_data[idx].real(0);
				spectral_space_data[idx].imag(0);
		);

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
	}

#endif



public:
	inline
	void request_data_spectral() const
	{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE

		FatalError("request_data_spectral: spectral space is disabled");

#else

#if SWEET_DEBUG
	#if SWEET_THREADING
		if (omp_get_num_threads() > 1)
			FatalError("Threading race conditions likely");
	#endif
#endif

		if (spectral_space_data_valid)
			return;		// nothing to do

		// cast to writable PlaneData
		PlaneData *rw_array_data = (PlaneData*)this;

#if SWEET_DEBUG
		if (!physical_space_data_valid)
			FatalError("Spectral data not available! Is this maybe a non-initialized operator?");
#endif

		planeDataConfig->fft_physical_to_spectral(rw_array_data->physical_space_data, rw_array_data->spectral_space_data);

		rw_array_data->spectral_space_data_valid = true;
		rw_array_data->physical_space_data_valid = false;

#endif
	}


	inline
	void request_data_physical() const
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE

		if (physical_space_data_valid)
			return;		// nothing to do

		PlaneData *rw_array_data = (PlaneData*)this;

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		spectral_zeroAliasingModes();
#endif

#if SWEET_DEBUG
		if (!spectral_space_data_valid)
			FatalError("Physical data not available and no spectral data!");
#endif

		planeDataConfig->fft_spectral_to_physical(rw_array_data->spectral_space_data, rw_array_data->physical_space_data);

		rw_array_data->spectral_space_data_valid = false;
		rw_array_data->physical_space_data_valid = true;
#endif
	}



	inline
	PlaneData physical_query_return_one_if_positive()
	{
		request_data_physical();

		PlaneData out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] > 0 ? 1 : 0);
			);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
#endif
		return out;
	}



	inline
	PlaneData physical_query_return_value_if_positive()	const
	{
		request_data_physical();

		PlaneData out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] > 0 ? physical_space_data[idx] : 0);
			);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
#endif
		return out;
	}


	inline
	PlaneData physical_query_return_one_if_negative()	const
	{
		request_data_physical();

		PlaneData out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] < 0 ? 1 : 0);
			);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
#endif
		return out;
	}


	inline
	PlaneData physical_query_return_value_if_negative()	const
	{
		request_data_physical();

		PlaneData out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] < 0 ? physical_space_data[idx] : 0);
			);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
#endif
		return out;
	}



	/**
	 * return true, if any value is infinity
	 */
	bool reduce_boolean_all_finite() const
	{
		request_data_physical();

		bool isallfinite = true;

#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(&&:isallfinite)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			isallfinite = isallfinite && std::isfinite(physical_space_data[i]);


		return isallfinite;
	}



	/**
	 * return the maximum of all absolute values
	 */
	double reduce_maxAbs()	const
	{
		request_data_physical();

		double maxabs = -1;
#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(max:maxabs)
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
		request_data_physical();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum)
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
		request_data_physical();

		double sum = 0;
		double c = 0;

#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
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
		request_data_physical();

		double maxvalue = -std::numeric_limits<double>::max();

#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(max:maxvalue)
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
		request_data_physical();

		double minvalue = std::numeric_limits<double>::max();

#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(min:minvalue)
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
		request_data_physical();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum)
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
		request_data_physical();

		double sum = 0;
		double c = 0;
#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
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
		request_data_physical();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum)
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
		request_data_physical();

		double sum = 0;
		double c = 0;
#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
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
		request_data_physical();

		double sum = 0;
#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum)
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
		request_data_physical();

		double sum = 0.0;
		double c = 0.0;

#if SWEET_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
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


public:
	template <int S>
	void kernel_stencil_setup(
			const double i_kernel_array[S][S],
			double i_scale = 1.0
	)
	{
		((PlaneData_Kernels&)*this).kernel_stencil_setup(
				i_kernel_array,
				i_scale,

				planeDataConfig,
				physical_space_data
			);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		physical_space_data_valid = true;
		spectral_space_data_valid = false;

		request_data_spectral();
#endif
	}



#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Invert the application of a linear operator in spectral space.
	 * The operator is given in i_array_data
	 */
	inline
	PlaneData spectral_div_element_wise(
			const PlaneData &i_array_data	///< operator
	)	const
	{
		PlaneData out(planeDataConfig);

		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		// only makes sense, if this is an operator created in spectral space
		assert(i_array_data.spectral_space_data_valid == true);

		request_data_spectral();
		rw_array_data.request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(

				double ar = spectral_space_data[idx].real();
				double ai = spectral_space_data[idx].imag();
				double br = i_array_data.spectral_space_data[idx].real();
				double bi = i_array_data.spectral_space_data[idx].imag();

				double den = br*br+bi*bi;

				if (den == 0)
				{
					out.spectral_space_data[idx].real(0);
					out.spectral_space_data[idx].imag(0);
				}
				else
				{
					out.spectral_space_data[idx].real((ar*br + ai*bi)/den);
					out.spectral_space_data[idx].imag((ai*br - ar*bi)/den);
				}
		);

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

		return out;
	}
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

			PLANE_DATA_PHYSICAL_FOR_IDX(
					physical_space_data[idx] = i_dataArray.physical_space_data[idx];
			);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		}
		else
		{
			physical_space_data_valid = false;
		}

		if (i_dataArray.spectral_space_data_valid)
		{
			spectral_space_data_valid = true;

			PLANE_DATA_SPECTRAL_FOR_IDX(
					spectral_space_data[idx] = i_dataArray.spectral_space_data[idx];
				);
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

		request_data_spectral();
		rw_array_data.request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = spectral_space_data[idx]*i_array_data.spectral_space_data[idx];
		);

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

#else

		PlaneData &rw_array_data = (PlaneData&)i_array_data;

		kernel_apply(
				planeDataConfig->physical_data_size[0],
				planeDataConfig->physical_data_size[1],
				rw_array_data.physical_space_data,

				out.physical_space_data
			);
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

		request_data_spectral();
		rw_array_data.request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = spectral_space_data[idx] + i_array_data.spectral_space_data[idx];
			);

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

#else

		request_data_physical();
		rw_array_data.request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = physical_space_data[idx] + i_array_data.physical_space_data[idx];
			);
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

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		request_data_spectral();

		PlaneData out = *this;
		out.spectral_space_data[0] += i_value*(double)planeDataConfig->physical_array_data_number_of_elements;

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else
		PlaneData out = *this;

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = physical_space_data[idx]+i_value;
				);

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

		request_data_spectral();
		rw_array_data.request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				spectral_space_data[idx] += i_array_data.spectral_space_data[idx];
		);

		//spectral_space_data_valid = true;
		//physical_space_data_valid = false;

#else

		PLANE_DATA_PHYSICAL_FOR_IDX(
				physical_space_data[idx] += i_array_data.physical_space_data[idx];
			);

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

		request_data_spectral();

		double scale = planeDataConfig->spectral_array_data_number_of_elements;
		spectral_space_data[0] += i_value*scale;

		spectral_space_data_valid = true;
		physical_space_data_valid = false;

#else

		request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
			physical_space_data[idx] += i_value;
		);

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

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				spectral_space_data[idx] *= i_value;
			);

		spectral_space_data_valid = true;
		physical_space_data_valid = false;

#else

		request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
			physical_space_data[idx] *= i_value;
		);

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

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				spectral_space_data[idx] /= i_value;
			);

		spectral_space_data_valid = true;
		physical_space_data_valid = false;

#else

		request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
			physical_space_data[idx] /= i_value;
		);

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

		request_data_spectral();
		rw_array_data.request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				spectral_space_data[idx] -= i_array_data.spectral_space_data[idx];
			);

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

#else

		request_data_physical();
		rw_array_data.request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
				physical_space_data[idx] -= i_array_data.physical_space_data[idx];
			);

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

		request_data_spectral();
		rw_array_data.request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
			out.spectral_space_data[idx] = spectral_space_data[idx] - i_array_data.spectral_space_data[idx];
		);

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else

		request_data_physical();
		rw_array_data.request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = physical_space_data[idx]-i_array_data.physical_space_data[idx];
				);

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

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		request_data_spectral();
		PlaneData out = *this;

		double scale = planeDataConfig->spectral_data_size[0]*planeDataConfig->spectral_data_size[1];
		out.spectral_space_data[0] -= i_value*scale;

		//out.physical_space_data_valid = false;
		//out.spectral_space_data_valid = true;

#else
		PlaneData out(planeDataConfig);

		request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(out.physical_space_data[idx] = physical_space_data[idx]-i_value;)

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

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = -spectral_space_data[idx];
			);

		double scale = planeDataConfig->physical_data_size[0]*planeDataConfig->physical_data_size[1];
		out.spectral_space_data[0] = i_value*scale + out.spectral_space_data[0];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else

		request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = i_value - physical_space_data[idx];
			);

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

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = -spectral_space_data[idx];
			);

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

#else

		request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
			out.physical_space_data[idx] = -physical_space_data[idx];
		);

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

#if SWEET_USE_PLANE_SPECTRAL_SPACE
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (spectral_space_data_valid)
			spectral_zeroAliasingModes();
		if (i_array_data.spectral_space_data_valid)
			i_array_data.spectral_zeroAliasingModes();
#endif
#endif

		request_data_physical();
		i_array_data.request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = physical_space_data[idx]*i_array_data.physical_space_data[idx];
			);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;


		/*
		 * Request data to be in spectral space again.
		 * This automatically forces to cut off the aliasing modes when converting back to physical space
		 */
		out.request_data_spectral();
#endif

		return out;
	}



	/**
	 * Compute element-wise division
	 */
	inline
	PlaneData operator/(
			const PlaneData &i_array_data	///< this class times i_array_data
	)	const
	{
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (spectral_space_data_valid)
			spectral_zeroAliasingModes();

		if (i_array_data.spectral_space_data_valid)
			i_array_data.spectral_zeroAliasingModes();
#endif
#endif

		request_data_physical();
		i_array_data.request_data_physical();

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = physical_space_data[idx]/i_array_data.physical_space_data[idx];
		);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;


		/*
		 * Request data to be in spectral space again.
		 * This automatically forces to cut off the aliasing modes when converting back to physical space
		 */
		out.request_data_spectral();
#endif

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
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = spectral_space_data[idx]*i_value;
		);

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

#else

		PLANE_DATA_PHYSICAL_FOR_IDX(
			out.physical_space_data[idx] = physical_space_data[idx]*i_value;
		);
#endif

		return out;
	}


	/**
	 * Compute element-wise division
	 */
	inline
	const PlaneData operator/(
			const double &i_value
	)	const
	{
		PlaneData out(planeDataConfig);

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = spectral_space_data[idx] / i_value;
		);

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

#else

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = physical_space_data[idx] / i_value;
		);

#endif

		return out;
	}





#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Add scalar to all spectral modes
	 */
	inline
	PlaneData spectral_addScalarAll(
			const double &i_value
	)	const
	{
		PlaneData out(planeDataConfig);

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = spectral_space_data[idx] + i_value;
		);

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}


	/**
	 * Return Plane Array with all spectral coefficients a+bi --> 1/(a+bi)
	 */
	inline
	PlaneData spectral_invert()	const
	{
		PlaneData out(planeDataConfig);

		request_data_spectral();

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = 1.0/spectral_space_data[idx];
		);

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}


	inline
	void print_spectralData()	const
	{
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.request_data_spectral();

		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = rw_array_data.spectral_get(y, x);
				std::cout << "(" << value.real() << ", " << value.imag() << ")\t";
			}
			std::cout << std::endl;
		}
	}

	inline
	void print_spectralIndex()	const
	{
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.request_data_spectral();

		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = rw_array_data.spectral_get(y, x);
				std::cout << "(" << x << ", "<< y << ", "<< value.real() << ", " << value.imag() << ")\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

	}

	inline
	void print_spectralNonZero()	const
	{
		PlaneData &rw_array_data = (PlaneData&)*this;

		rw_array_data.request_data_spectral();

		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = rw_array_data.spectral_get(y, x);
				if (value.real()*value.real()+value.imag()*value.imag() > 1.0e-13)
					std::cout << "(" << x << ", "<< y << ", "<< value.real() << ", " << value.imag() << ")" <<std::endl;;
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
	bool print_physicalArrayData(
			int i_precision = 8		///< number of floating point digits
			)	const
	{
		request_data_physical();

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


	friend
	inline
	std::ostream& operator<<(
			std::ostream &o_ostream,
			const PlaneData &i_dataArray
	)
	{
		PlaneData &rw_array_data = (PlaneData&)i_dataArray;

		rw_array_data.request_data_physical();

		for (int j = rw_array_data.planeDataConfig->physical_data_size[1]-1; j >= 0; j--)
		{
			for (std::size_t i = 0; i < rw_array_data.planeDataConfig->physical_data_size[0]; i++)
			{
				std::cout << i_dataArray.physical_space_data[j*i_dataArray.planeDataConfig->physical_data_size[0]+i] << "\t";
			}
			std::cout << std::endl;
		}

		return o_ostream;
	}

	/**
	 * Write data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_physical_saveData_ascii(
			const char *i_filename,		///< Name of file to store data to
			char i_separator = '\t',	///< separator to use for each line
			int i_precision = 12		///< number of floating point digits
	)	const
	{
		request_data_physical();

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



	/**
	 * Write data to VTK file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_physical_saveData_vtk(
			const char *i_filename,		///< Name of file to store data to
			const char *i_title,		///< Title of scalars
			int i_precision = 12		///< number of floating point digits
	)	const
	{
		request_data_physical();

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
	bool file_physical_loadData(
			const char *i_filename,		///< Name of file to load data from
			bool i_binary_data = false	///< load as binary data (disabled per default)
	)
	{
		if (i_binary_data)
		{
			std::ifstream file(i_filename, std::ios::binary);

			if (!file)
				FatalError(std::string("Failed to open file ")+i_filename);

			file.seekg(0, std::ios::end);
			std::size_t size = file.tellg();
			file.seekg(0, std::ios::beg);


			std::size_t expected_size = sizeof(double)*planeDataConfig->physical_res[0]*planeDataConfig->physical_res[1];

			if (size != expected_size)
			{
				std::cerr << "Error while loading data from file " << i_filename << ":" << std::endl;
				std::cerr << "Size of file " << size << " does not match expected size of " << expected_size << std::endl;
				FatalError("EXIT");
			}

			if (!file.read((char*)physical_space_data, expected_size))
			{
				std::cerr << "Error while loading data from file " << i_filename << std::endl;
				FatalError("EXIT");
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

				p_physical_set(planeDataConfig->physical_res[1]-row-1, col, i_value);

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
