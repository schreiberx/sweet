/*
 * PlaneDataPhysical.hpp
 *
 *  Created on: 14 Mar 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef PLANE_DATA_PHYSICAL_HPP_
#define PLANE_DATA_PHYSICAL_HPP_


#include <complex>
#include <cfloat>
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
#include <functional>
#include <utility>
#include <cmath>


#include <sweet/MemBlockAlloc.hpp>
#include <sweet/openmp_helper.hpp>
#include <sweet/plane/PlaneDataConfig.hpp>
#include <sweet/SWEETError.hpp>
#include <sweet/plane/PlaneData_Kernels.hpp>

//#include <sweet/plane/PlaneData_Spectral.hpp>
class PlaneData_Spectral;

#if SWEET_THREADING_SPACE
#define PLANE_DATA_PHYSICAL_FOR_IDX(CORE)				\
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD			\
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)	\
		{	CORE;	}

#define PLANE_DATA_PHYSICAL_FOR_IDX_REDUCTION(CORE, REDUCTION)				\
		_Pragma("omp parallel for "##REDUCTION##" "##PROC_BIND_CLOSE##"")	\
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)	\
		{	CORE;	}

#define PLANE_DATA_PHYSICAL_FOR_2D_IDX(CORE)										\
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2	\
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

#if SWEET_THREADING_SPACE
#	include <omp.h>
#endif


class PlaneData_Physical
			:	private PlaneData_Kernels
{
////	friend class PlaneData_Spectral;
////	friend PlaneData_Physical PlaneData_Spectral::multiplication_physical_space(PlaneData_Physical& i_a, PlaneData_Physical& i_b);

public:
	const PlaneDataConfig *planeDataConfig;

public:
	double *physical_space_data;


	void swap(
			PlaneData_Physical &i_planeData
	)
	{
		assert(planeDataConfig == i_planeData.planeDataConfig);

		std::swap(physical_space_data, i_planeData.physical_space_data);
	}

public:
	PlaneData_Physical(
			const PlaneDataConfig *i_planeDataConfig
	)	:
		/// important: set this to nullptr, since a check for this will be performed by setup(...)
		planeDataConfig(i_planeDataConfig),
		physical_space_data(nullptr)
	{
		alloc_data();
	}


public:
	PlaneData_Physical(
			const PlaneDataConfig *i_planeDataConfig,
			double i_value
	)	:
		/// important: set this to nullptr, since a check for this will be performed by setup(...)
		planeDataConfig(i_planeDataConfig),
		physical_space_data(nullptr)
	{
		alloc_data();
		physical_set_all_value(i_value);
	}


public:
	PlaneData_Physical()	:
		planeDataConfig(nullptr),
		physical_space_data(nullptr)
	{
	}

	/**
	 * dummy initialization by handing over an unused integer
	 */
public:
	PlaneData_Physical(int i)	:
		planeDataConfig(nullptr),
		physical_space_data(nullptr)
{
}


public:
	PlaneData_Physical(
			const PlaneData_Physical &i_plane_data
	)	:
		planeDataConfig(i_plane_data.planeDataConfig),
		physical_space_data(nullptr)

	{
		if (i_plane_data.planeDataConfig != nullptr)
			alloc_data();

		operator=(i_plane_data);
	}


public:
	PlaneData_Physical(
			PlaneData_Physical &&i_plane_data
	)	:
		planeDataConfig(i_plane_data.planeDataConfig),
		physical_space_data(nullptr)
	{
		if (i_plane_data.planeDataConfig == nullptr)
			return;

		physical_space_data = i_plane_data.physical_space_data;
		i_plane_data.physical_space_data = nullptr;
	}



	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(
			const PlaneDataConfig *i_planeDataConfig
	)	const
	{
		assert(planeDataConfig->physical_res[0] == i_planeDataConfig->physical_res[0]);
		assert(planeDataConfig->physical_res[1] == i_planeDataConfig->physical_res[1]);
	}



public:
	PlaneData_Physical& operator=(
			const PlaneData_Physical &i_plane_data
	)
	{
		if (i_plane_data.planeDataConfig == nullptr)
			return *this;

		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		memcpy(physical_space_data, i_plane_data.physical_space_data, sizeof(double)*planeDataConfig->physical_array_data_number_of_elements);

		this->dealiasing(*this);

		return *this;
	}


public:
	PlaneData_Physical& operator=(
			PlaneData_Physical &&i_plane_data
	)
	{
		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		std::swap(physical_space_data, i_plane_data.physical_space_data);

		this->dealiasing(*this);

		return *this;
	}


public:
	/**
	 * assignment operator
	 */
	PlaneData_Physical &operator=(double i_value)
	{
		physical_set_all_value(i_value);

		return *this;
	}


public:
	/**
	 * assignment operator
	 */
	PlaneData_Physical &operator=(int i_value)
	{
		physical_set_all_value(i_value);

		return *this;
	}



	PlaneData_Physical operator+(
			const PlaneData_Physical &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] + i_plane_data.physical_space_data[idx];

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		this->dealiasing(out);
#endif

		return out;
	}



	PlaneData_Physical& operator+=(
			const PlaneData_Physical &i_plane_data
	)
	{
		check(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_plane_data.physical_space_data[idx];

		return *this;
	}


	PlaneData_Physical& operator+=(
			double i_scalar
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_scalar;

		return *this;
	}


	PlaneData_Physical& operator-=(
			const PlaneData_Physical &i_plane_data
	)
	{
		check(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_plane_data.physical_space_data[idx];

		return *this;
	}


	PlaneData_Physical& operator-=(
			double i_scalar
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_scalar;

		return *this;
	}



	PlaneData_Physical operator-(
			const PlaneData_Physical &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] - i_plane_data.physical_space_data[idx];

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		this->dealiasing(out);
#endif

		return out;
	}



	PlaneData_Physical operator-()
	{
		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = -physical_space_data[idx];

		return out;
	}



	PlaneData_Physical operator*(
			const PlaneData_Physical &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);


		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_plane_data.physical_space_data[i];


#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		this->dealiasing(out);
#endif

		return out;

	}

	void dealiasing(
		PlaneData_Physical& io_data
	) const
	{
		
		// Ugly fix to circular dependency between PlaneData_Physical and PlaneData_Spectral

		// create spectral data container
		std::complex<double> *spectral_space_data = nullptr;
		spectral_space_data = MemBlockAlloc::alloc<std::complex<double>>(planeDataConfig->spectral_array_data_number_of_elements * sizeof(std::complex<double>));

		// FFT
		planeDataConfig->fft_physical_to_spectral(io_data.physical_space_data, spectral_space_data);

		// Dealiasing
		//SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int k = 0; k < 2; k++)
		{
			if (k == 0)
			{
				/*
				 * First process part between top and bottom spectral data blocks
				 */
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[0][1][1]; jj < planeDataConfig->spectral_data_iteration_ranges[1][1][0]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii++)
					{
						//spectral_space_data[jj*planeDataConfig->spectral_data_size[0]+ii] = 0;
						spectral_space_data[planeDataConfig->getArrayIndexByModes(jj, ii)] = 0;
					}
			}
			else
			{
				/*
				 * Then process the aliasing block on the right side
				 */
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (std::size_t jj = 0; jj < planeDataConfig->spectral_data_size[1]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii < planeDataConfig->spectral_data_size[0]; ii++)
					{
						//spectral_space_data[jj*planeDataConfig->spectral_data_size[0]+ii] = 0;
						spectral_space_data[planeDataConfig->getArrayIndexByModes(jj, ii)] = 0;
					}
			}
		}

		// IFFT
		planeDataConfig->fft_spectral_to_physical(spectral_space_data, io_data.physical_space_data);

		// Free spectral data
		MemBlockAlloc::free(spectral_space_data, planeDataConfig->spectral_array_data_number_of_elements * sizeof(std::complex<double>));
		spectral_space_data = nullptr;
	}


	PlaneData_Physical multiplication_no_dealiasing(
			const PlaneData_Physical &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_plane_data.physical_space_data[i];

		return out;
	}



	PlaneData_Physical operator/(
			const PlaneData_Physical &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		check(i_plane_data.planeDataConfig);

		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]/i_plane_data.physical_space_data[i];

		this->dealiasing(out);

		return out;
	}



	PlaneData_Physical operator*(
			const double i_value
	)	const
	{
		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_value;

		return out;
	}




	const PlaneData_Physical& operator*=(
			const double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] *= i_value;

		return *this;
	}




	PlaneData_Physical operator/(
			double i_value
	)	const
	{
		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]/i_value;

		return out;
	}

	PlaneData_Physical& operator/=(
			double i_scalar
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] /= i_scalar;

		return *this;
	}




	PlaneData_Physical operator+(
			double i_value
	)	const
	{
		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]+i_value;

		return out;
	}



	PlaneData_Physical operator-(
			double i_value
	)	const
	{
		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]-i_value;


		return out;
	}


	PlaneData_Physical operator_scalar_sub_this(
			double i_value
	)	const
	{
		PlaneData_Physical out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = i_value - physical_space_data[idx];

		return out;
	}


	/**
	 * Apply a linear operator given by this class to the input data array.
	 */
	inline
	PlaneData_Physical operator()(
			const PlaneData_Physical &i_array_data
	)	const
	{
		PlaneData_Physical out(planeDataConfig);

		PlaneData_Physical &rw_array_data = (PlaneData_Physical&)i_array_data;

		kernel_apply(
				planeDataConfig->physical_data_size[0],
				planeDataConfig->physical_data_size[1],
				rw_array_data.physical_space_data,

				out.physical_space_data
		);


		return out;
	}

	friend
	inline
	std::ostream& operator<<(
			std::ostream &o_ostream,
			const PlaneData_Physical &i_dataArray
	)
	{
		PlaneData_Physical &rw_array_data = (PlaneData_Physical&)i_dataArray;

		for (int j = (int)rw_array_data.planeDataConfig->physical_data_size[1]-1; j >= 0; j--)
		{
			for (std::size_t i = 0; i < rw_array_data.planeDataConfig->physical_data_size[0]; i++)
			{
				std::cout << i_dataArray.physical_space_data[j*i_dataArray.planeDataConfig->physical_data_size[0]+i] << "\t";
			}
			std::cout << std::endl;
		}

		return o_ostream;
	}




public:
	void setup(
			const PlaneDataConfig *i_planeDataConfig
	)
	{
		if (planeDataConfig != nullptr)
			SWEETError("Setup called twice!");

		planeDataConfig = i_planeDataConfig;
		alloc_data();
	}



private:
	void alloc_data()
	{
		assert(physical_space_data == nullptr);
		physical_space_data = MemBlockAlloc::alloc<double>(planeDataConfig->physical_array_data_number_of_elements * sizeof(double));
	}




public:
	void setup_if_required(
			const PlaneDataConfig *i_planeDataConfig
	)
	{
		if (planeDataConfig != nullptr)
			return;

		setup(i_planeDataConfig);
	}



public:
	~PlaneData_Physical()
	{
		free();
	}


public:
	void free()
	{
		if (physical_space_data != nullptr)
		{
			MemBlockAlloc::free(physical_space_data, planeDataConfig->physical_array_data_number_of_elements * sizeof(double));
			physical_space_data = nullptr;
		}

		planeDataConfig = nullptr;
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

	}


public:
	void physical_update_lambda_array_indices(
			std::function<void(int,int,double&)> i_lambda,	///< lambda function to return value for lat/mu
			bool anti_aliasing = true
	)
	{
		PLANE_DATA_PHYSICAL_FOR_2D_IDX(
				i_lambda(i, j, physical_space_data[idx])
		);
	}


	void physical_update_lambda_array_idx(
			std::function<void(int,double&)> i_lambda,	///< lambda function to return value for lat/mu
			bool anti_aliasing = true
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR

		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			i_lambda(i, physical_space_data[i]);
		}
	}


	void physical_update_lambda_unit_coordinates_corner_centered(
			std::function<void(double,double,double&)> i_lambda,	///< lambda function to return value for lat/mu
			bool anti_aliasing = true
	)
	{
		PLANE_DATA_PHYSICAL_FOR_2D_IDX(
				i_lambda(
						(double)i/(double)planeDataConfig->physical_res[0],
						(double)j/(double)planeDataConfig->physical_res[1],
						physical_space_data[idx]
				)
		);
	}

	void physical_update_lambda_unit_coordinates_cell_centered(
			std::function<void(double,double,double&)> i_lambda,	///< lambda function to return value for lat/mu
			bool anti_aliasing = true
	)
	{
		PLANE_DATA_PHYSICAL_FOR_2D_IDX(
				i_lambda(
						((double)i+0.5)/(double)planeDataConfig->physical_res[0],
						((double)j+0.5)/(double)planeDataConfig->physical_res[1],
						physical_space_data[idx]
				)
		);
	}





	/*
	 * Set all values to zero
	 */
	void physical_set_zero()
	{
		PLANE_DATA_PHYSICAL_FOR_2D_IDX(
				physical_space_data[idx] = 0;
				)
		///SWEET_THREADING_SPACE_PARALLEL_FOR

		///for (int i = 0; i < planeDataConfig->physical_num_lon; i++)
		///	for (int j = 0; j < planeDataConfig->physical_num_lat; j++)
		///		physical_space_data[j*planeDataConfig->physical_num_lon + i] = 0;
	}



	/*
	 * Set all values to a specific value
	 */
	void physical_set_all_value(
			double i_value
	)
	{
		PLANE_DATA_PHYSICAL_FOR_2D_IDX(
				physical_space_data[idx] = i_value;
				)
		////SWEET_THREADING_SPACE_PARALLEL_FOR
		////for (int i = 0; i < planeDataConfig->physical_num_lon; i++)
		////	for (int j = 0; j < planeDataConfig->physical_num_lat; j++)
		////		physical_space_data[j*planeDataConfig->physical_num_lon + i] = i_value;
	}



	/*
	 * Set given point a specific value
	 */
	void physical_set_value(
			int i_y_idx,
			int i_x_idx,
			double i_value
	)
	{
		physical_space_data[i_y_idx*planeDataConfig->physical_res[0] + i_x_idx] = i_value;
	}

	const double physical_get(
			int i_y_idx,
			int i_x_idx
	)	const
	{
		return physical_space_data[i_y_idx*planeDataConfig->physical_res[0] + i_x_idx];
	}


	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_max(
			const PlaneData_Physical &i_plane_data
	)
	{
		check(i_plane_data.planeDataConfig);

		double error = -1;

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						std::abs(
								physical_space_data[j] - i_plane_data.physical_space_data[j]
							),
							error
						);
		}
		return error;
	}


	double physical_reduce_rms()
	{
		double error = 0;

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			double &d = physical_space_data[j];
			error += d*d;
		}

		return std::sqrt(error / (double)planeDataConfig->physical_array_data_number_of_elements);
	}



	double physical_reduce_sum()	const
	{
		double sum = 0;
		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
			sum += physical_space_data[j];

		return sum;
	}


	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double physical_reduce_sum_quad()	const
	{
		double sum = 0;
		double c = 0;
#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
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
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double physical_reduce_sum_quad_increasing()	const
	{
		double sum = 0;
		double c = 0;
#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			double value = physical_space_data[i]*(double)i;

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return sum;
	}


	double physical_reduce_sum_metric()
	{
		SWEETError("TODO: Implement metric-scaled summation");
		double sum = 0;
		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			sum += physical_space_data[j];
		}

		return sum;
	}

	/**
	 * reduce to root mean square
	 */
	double physical_reduce_rms_quad()
	{

		double sum = 0;
		double c = 0;

#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
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
	 * return the sqrt of the sum of the squared values
	 */
	double physical_reduce_norm2()	const
	{
		double sum = 0;
#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			sum += physical_space_data[i]*physical_space_data[i];


		return std::sqrt(sum);
	}


	/**
	 * return the sqrt of the sum of the squared values, use quad precision for reduction
	 */
	double physical_reduce_norm2_quad()	const
	{
		double sum = 0.0;
		double c = 0.0;

#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
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




	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_max_abs(
			const PlaneData_Physical &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		double error = -1;

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						std::abs(
								physical_space_data[j] - i_plane_data.physical_space_data[j]
							),
							error	// leave the error variable as the 2nd parameter. In case of NaN of the 1st parameter, std::max returns NaN
						);
		}
		return error;
	}



	/**
	 * Return the maximum absolute value
	 */
	double physical_reduce_max_abs()	const
	{
		double error = -1;

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						std::abs(physical_space_data[j]),
						error		// leave the error variable as the 2nd parameter. In case of NaN of the 1st parameter, std::max returns NaN
				);
		}
		return error;
	}


	/**
	 * Return the minimum value
	 */
	double physical_reduce_min()	const
	{
		double error = std::numeric_limits<double>::infinity();

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
			error = std::min(physical_space_data[j], error);

		return error;
	}


	/**
	 * Return the minimum value
	 */
	double physical_reduce_max()	const
	{
		double error = -std::numeric_limits<double>::infinity();

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
			error = std::max(physical_space_data[j], error);

		return error;
	}


	/**
	 * Rescale data so that max_abs returns the given value
	 */
	PlaneData_Physical physical_rescale_to_max_abs(
			double i_new_max_abs
	)	const
	{
		double max_abs = physical_reduce_max_abs();
		double scale = i_new_max_abs/max_abs;

		PlaneData_Physical out(planeDataConfig);

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
			out.physical_space_data[j] = physical_space_data[j]*scale;

		return out;
	}



	bool physical_isAnyNaNorInf()	const
	{
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			if (std::isnan(physical_space_data[i]) || std::isinf(physical_space_data[i]) != 0)
				return true;
		}

		return false;
	}

	inline
	PlaneData_Physical physical_query_return_one_if_positive()
	{
		PlaneData_Physical out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] > 0 ? 1 : 0);
		);

		return out;
	}



	inline
	PlaneData_Physical physical_query_return_value_if_positive()	const
	{
		PlaneData_Physical out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] > 0 ? physical_space_data[idx] : 0);
		);

		return out;
	}


	inline
	PlaneData_Physical physical_query_return_one_if_negative()	const
	{
		PlaneData_Physical out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] < 0 ? 1 : 0);
		);

		return out;
	}


	inline
	PlaneData_Physical physical_query_return_value_if_negative()	const
	{
		PlaneData_Physical out(planeDataConfig);

		PLANE_DATA_PHYSICAL_FOR_IDX(
				out.physical_space_data[idx] = (physical_space_data[idx] < 0 ? physical_space_data[idx] : 0);
		);

		return out;
	}



	/**
	 * return true, if any value is infinity
	 */
	bool physical_reduce_boolean_all_finite() const
	{

		bool isallfinite = true;

#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(&&:isallfinite)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			isallfinite = isallfinite && std::isfinite(physical_space_data[i]);

		return isallfinite;
	}


	void physical_print(
			int i_precision = -1
	)	const
	{
		if (i_precision >= 0)
			std::cout << std::setprecision(i_precision);

        for (int j = (int)(planeDataConfig->physical_res[1]-1); j >= 0; j--)
        {
        		for (int i = 0; i < (int)planeDataConfig->physical_res[0]; i++)
        		{
        			std::cout << physical_space_data[j*planeDataConfig->physical_res[0]+i];
        			if (i < (int)planeDataConfig->physical_res[0]-1)
        				std::cout << "\t";
        		}
        		std::cout << std::endl;
        }
	}

	/**
	 * print spectral data and zero out values which are numerically close to zero
	 */
	inline
	void print_physicalData_zeroNumZero(double i_zero_threshold = 1e-13)	const
	{
		PlaneData_Physical &rw_array_data = (PlaneData_Physical&)*this;

		for (int y = planeDataConfig->physical_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				double value = rw_array_data.physical_get(y, x);

				if (std::abs(value) < i_zero_threshold)
					value = 0.0;

				std::cout << value;

				if (x != planeDataConfig->physical_data_size[0]-1)
					std::cout << "\t";
			}
			std::cout << std::endl;
		}
	}




	void physical_file_write(
			const std::string &i_filename,
			double plane_domain_size[2],
			const char *i_title = "",
			int i_precision = 20
	)	const
	{
		std::ofstream file(i_filename, std::ios_base::trunc);

		if (i_precision >= 0)
			file << std::setprecision(i_precision);

		file << "#TI " << i_title << std::endl;
		file << "#TX" << std::endl;
		file << "#TY" << std::endl;

		//file << "lat\\lon\t";
		// Use 0 to make it processable by python
		file << "0\t";

		for (int i = 0; i < (int)planeDataConfig->physical_res[0]; i++)
		{
//			double lon_degree = ((double)i/(double)planeDataConfig->spat_num_lon)*2.0*M_PI;
			double x = ((double)i/(double)planeDataConfig->physical_res[0])*plane_domain_size[0]; // ????

			file << x;
			if (i < (int)planeDataConfig->physical_res[0]-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = (int)planeDataConfig->physical_res[1]-1; j >= 0; j--)
        {
			double y = ((double)j/(double)planeDataConfig->physical_res[1])*plane_domain_size[1]; // ????

        		file << y << "\t";

        		for (int i = 0; i < (int)planeDataConfig->physical_res[0]; i++)
        		{
        			file << physical_space_data[j*planeDataConfig->physical_res[0]+i];
        			if (i < (int)planeDataConfig->physical_res[0]-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}




	bool physical_file_load(
			const std::string &i_filename,		///< Name of file to load data from
			bool i_binary_data = false	///< load as binary data (disabled per default)
	)
	{
		if (i_binary_data)
		{
			std::ifstream file(i_filename, std::ios::binary);

			if (!file)
				SWEETError(std::string("Failed to open file ")+i_filename);

			file.seekg(0, std::ios::end);
			std::size_t size = file.tellg();
			file.seekg(0, std::ios::beg);


			std::size_t expected_size = sizeof(double)*planeDataConfig->physical_array_data_number_of_elements;

			if (size != expected_size)
			{
				std::cerr << "Error while loading data from file " << i_filename << ":" << std::endl;
				std::cerr << "Size of file " << size << " does not match expected size of " << expected_size << std::endl;
				SWEETError("EXIT");
			}

			if (!file.read((char*)physical_space_data, expected_size))
			{
				std::cerr << "Error while loading data from file " << i_filename << std::endl;
				SWEETError("EXIT");
			}

			return true;
		}


		std::ifstream file(i_filename);
		std::string line;

		bool first_data_line = true;


		/*
		 * set physical data to be valid right here!!!
		 * otherwise it might happen that physical_set_value always requests
		 * data to be converted to physical data, hence overwriting data
		 */

		int row = 0;
		while (row < (int)planeDataConfig->physical_res[1])
		{
			std::getline(file, line);
			if (!file.good())
			{
				std::cerr << "ERROR: EOF - Failed to read data from file " << i_filename << " in line " << row << std::endl;
				return false;
			}

			// skip comment lines
			if (line[0] == '#')
				continue;

			int last_pos = 0;
			int col = 0;

			if (first_data_line)
			{
				// skip first data line since these are the coordinates
				first_data_line = false;
				continue;
			}

			for (int pos = 0; pos < (int)line.size()+1; pos++)
			{
				if (pos < (int)line.size())
					if (line[pos] != '\t' && line[pos] != ' ')
						continue;

				// skip first element!
				if (last_pos > 0)
				{
					std::string strvalue = line.substr(last_pos, pos-last_pos);
					double i_value = atof(strvalue.c_str());
					int x = col;
					int y = planeDataConfig->physical_res[1]-row-1;

					if (x > (int)planeDataConfig->physical_res[0])
						return false;

					if (y > (int)planeDataConfig->physical_res[1])
						return false;

					physical_set_value(y, x, i_value);

					col++;
				}

				last_pos = pos+1;
		    }

			if (col != (int)planeDataConfig->physical_res[0])
			{
				std::cerr << "ERROR: column mismatch - Failed to read data from file " << i_filename << " in line " << row << ", column " << col << std::endl;
				return false;
			}

			row++;
		}

		if (row != (int)planeDataConfig->physical_res[1])
		{
			std::cerr << "ERROR: rows mismatch - Failed to read data from file " << i_filename << " in line " << row << std::endl;
			return false;
		}


		return true;
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
			int i_precision = 16,		///< number of floating point digits
			int dimension = 2			///< store 1D or 2D
	)	const
	{

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		file << "#SWEET" << std::endl;
		file << "#FORMAT ASCII" << std::endl;
		file << "#PRIMITIVE PLANE" << std::endl;

		std::size_t resx = planeDataConfig->physical_res[0];
		std::size_t resy = planeDataConfig->physical_res[1];

		file << "#SPACE PHYSICAL" << std::endl;
		file << "#RESX " << resx << std::endl;
		file << "#RESY " << resy << std::endl;

		std::size_t ymin = 0;
		if (dimension == 2)
			ymin = 0;
		else
			ymin = planeDataConfig->physical_res[1]-1;

		for (int y = (int) resy-1; y >= (int) ymin; y--)
		{
			for (std::size_t x = 0; x < resx; x++)
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
				SWEETError(std::string("Failed to open file ")+i_filename);

			file.seekg(0, std::ios::end);
			std::size_t size = file.tellg();
			file.seekg(0, std::ios::beg);


			std::size_t expected_size = sizeof(double)*planeDataConfig->physical_res[0]*planeDataConfig->physical_res[1];

			if (size != expected_size)
			{
				std::cerr << "Error while loading data from file " << i_filename << ":" << std::endl;
				std::cerr << "Size of file " << size << " does not match expected size of " << expected_size << std::endl;
				SWEETError("EXIT");
			}

			if (!file.read((char*)physical_space_data, expected_size))
			{
				std::cerr << "Error while loading data from file " << i_filename << std::endl;
				SWEETError("EXIT");
			}

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

				physical_set_value(planeDataConfig->physical_res[1]-row-1, col, i_value);

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



	void file_write_raw(
			const std::string &i_filename
	)	const
	{
		std::fstream file(i_filename, std::ios::out | std::ios::binary);
		file.write((const char*)physical_space_data, sizeof(double)*planeDataConfig->physical_array_data_number_of_elements);
	}



	void file_read_raw(
			const std::string &i_filename
	)	const
	{
		std::fstream file(i_filename, std::ios::in | std::ios::binary);
		file.read((char*)physical_space_data, sizeof(double)*planeDataConfig->physical_array_data_number_of_elements);
	}


	void print_debug(
			const char *name
	)	const
	{
		std::cout << name << ":" << std::endl;
		std::cout << "                min: " << this->physical_reduce_min() << std::endl;
		std::cout << "                max: " << this->physical_reduce_max() << std::endl;
		std::cout << "                sum: " << this->physical_reduce_sum() << std::endl;
		std::cout << "                suminc: " << this->physical_reduce_sum_quad_increasing() << std::endl;
		std::cout << std::endl;
	}


	void print()	const
	{
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			std::cout << physical_space_data[idx] << "\t";
		std::cout << std::endl;
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
PlaneData_Physical operator*(
		const double i_value,
		const PlaneData_Physical &i_array_data
)
{
	return i_array_data*i_value;
}



/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 */
inline
static
PlaneData_Physical operator-(
		const double i_value,
		const PlaneData_Physical &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}



/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 */
inline
static
PlaneData_Physical operator+(
		const double i_value,
		const PlaneData_Physical &i_array_data
)
{
	return i_array_data+i_value;
}





#endif
