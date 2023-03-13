/*
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PLANE_DATA_PHYSICAL_COMPLEX_HPP_
#define SRC_PLANE_DATA_PHYSICAL_COMPLEX_HPP_

#include <complex>
#include <functional>
#include <array>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <limits>
#include <utility>
#include <functional>

#include <cmath>
#include <sweet/core/MemBlockAlloc.hpp>
#include <sweet/core/openmp_helper.hpp>
#include <sweet/core/plane/PlaneData_Config.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/SWEETError.hpp>

#if SWEET_THREADING_SPACE
	#define PLANE_DATA_COMPLEX_PHYSICAL_FOR_IDX(CORE)				\
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD	\
			for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)	\
			{	CORE;	}

	#define PLANE_DATA_COMPLEX_PHYSICAL_FOR_2D_IDX(CORE)										\
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

	#define PLANE_DATA_COMPLEX_PHYSICAL_FOR_IDX(CORE)				\
			for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)	\
			{	CORE;	}

	#define PLANE_DATA_COMPLEX_PHYSICAL_FOR_2D_IDX(CORE)										\
			for (std::size_t j = 0; j < planeDataConfig->physical_data_size[1]; j++)						\
			{				\
				for (std::size_t i = 0; i < planeDataConfig->physical_data_size[0]; i++)	\
				{			\
					std::size_t idx = j*planeDataConfig->physical_data_size[0]+i;	\
					CORE;	\
				}			\
			}

#endif

namespace sweet
{

class PlaneData_PhysicalComplex
{

public:
	const PlaneData_Config *planeDataConfig;

public:
	std::complex<double> *physical_space_data;


	void swap(
			PlaneData_PhysicalComplex &i_planeData
	)
	{
		assert(planeDataConfig == i_planeData.planeDataConfig);

		std::swap(physical_space_data, i_planeData.physical_space_data);
	}

public:
	PlaneData_PhysicalComplex(
			const PlaneData_Config *i_planeDataConfig
	)	:
		planeDataConfig(i_planeDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_planeDataConfig);
	}


public:
	PlaneData_PhysicalComplex()	:
		planeDataConfig(nullptr),
		physical_space_data(nullptr)
	{
	}


public:
	PlaneData_PhysicalComplex(
			const PlaneData_PhysicalComplex &i_plane_data
	)	:
		planeDataConfig(i_plane_data.planeDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_plane_data.planeDataConfig);

		operator=(i_plane_data);
	}


public:
	PlaneData_PhysicalComplex(
			const PlaneData_Physical &i_plane_data
	)	:
		planeDataConfig(i_plane_data.planeDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_plane_data.planeDataConfig);

		operator=(i_plane_data);
	}


	/*
	 * load real and imaginary data from physical arrays
	 */
	void loadRealImag(
			const PlaneData_Physical &i_re,
			const PlaneData_Physical &i_im
	)
	{
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
#if 0
			physical_space_data[i].real(i_re.physical_space_data[i]);
			physical_space_data[i].imag(i_im.physical_space_data[i]);
#else
			physical_space_data[i] = std::complex<double>(
								i_re.physical_space_data[i],
								i_im.physical_space_data[i]
						);
#endif
		}
	}


	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(
			const PlaneData_Config *i_planeDataConfig
	)	const
	{
		assert(planeDataConfig->physical_res[0] == i_planeDataConfig->physical_res[0]);
		assert(planeDataConfig->physical_res[1] == i_planeDataConfig->physical_res[1]);
	}



public:
	PlaneData_PhysicalComplex& operator=(
			const PlaneData_PhysicalComplex &i_plane_data
	)
	{
		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		memcpy(physical_space_data, i_plane_data.physical_space_data, sizeof(std::complex<double>)*planeDataConfig->physical_array_data_number_of_elements);

		return *this;
	}


public:
	PlaneData_PhysicalComplex& operator=(
			PlaneData_PhysicalComplex &&i_plane_data
	)
	{
		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		std::swap(physical_space_data, i_plane_data.physical_space_data);

		return *this;
	}



public:
	PlaneData_PhysicalComplex& operator=(
			const PlaneData_Physical &i_plane_data
	)
	{
		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_plane_data.physical_space_data[i];

		return *this;
	}





public:
	void setup(
			const PlaneData_Config *i_planeDataConfig
	)
	{
		planeDataConfig = i_planeDataConfig;

		physical_space_data = MemBlockAlloc::alloc<std::complex<double>>(planeDataConfig->physical_array_data_number_of_elements * sizeof(std::complex<double>));
	}



public:
	void setup_if_required(
			const PlaneData_Config *i_planeDataConfig
	)
	{
		if (planeDataConfig != nullptr)
		{
			assert(physical_space_data != nullptr);
			return;
		}

		planeDataConfig = i_planeDataConfig;
		physical_space_data = MemBlockAlloc::alloc<std::complex<double>>(planeDataConfig->physical_array_data_number_of_elements * sizeof(std::complex<double>));
	}



public:
	~PlaneData_PhysicalComplex()
	{
		if (physical_space_data != nullptr)
			MemBlockAlloc::free(physical_space_data, planeDataConfig->physical_array_data_number_of_elements * sizeof(std::complex<double>));
	}


public:
	PlaneData_PhysicalComplex operator+(
			const PlaneData_PhysicalComplex &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] + i_plane_data.physical_space_data[idx];

		return out;
	}



	PlaneData_PhysicalComplex& operator+=(
			const PlaneData_PhysicalComplex &i_plane_data
	)
	{
		check(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_plane_data.physical_space_data[idx];

		return *this;
	}


	PlaneData_PhysicalComplex& operator-=(
			const PlaneData_PhysicalComplex &i_plane_data
	)
	{
		check(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_plane_data.physical_space_data[idx];

		return *this;
	}



	PlaneData_PhysicalComplex operator-(
			const PlaneData_PhysicalComplex &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] - i_plane_data.physical_space_data[idx];

		return out;
	}



	PlaneData_PhysicalComplex operator-()
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = -physical_space_data[idx];

		return out;
	}



	PlaneData_PhysicalComplex operator*(
			const PlaneData_PhysicalComplex &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_plane_data.physical_space_data[i];

		return out;
	}



	PlaneData_PhysicalComplex operator/(
			const PlaneData_PhysicalComplex &i_plane_data
	)	const
	{
		check(i_plane_data.planeDataConfig);

		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			out.physical_space_data[i] = physical_space_data[i]/i_plane_data.physical_space_data[i];
		}

		return out;
	}



	PlaneData_PhysicalComplex operator*(
			const double i_value
	)	const
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_value;

		return out;
	}



	PlaneData_PhysicalComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]*i_value;

		return out;
	}




	const PlaneData_PhysicalComplex& operator*=(
			const double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] *= i_value;

		return *this;
	}



	PlaneData_PhysicalComplex operator/(
			double i_value
	)	const
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]/i_value;

		return out;
	}



	PlaneData_PhysicalComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]+i_value;

		return out;
	}



	PlaneData_PhysicalComplex operator-(
			const std::complex<double> &i_value
	)	const
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]-i_value;

		return out;
	}


	friend
	inline
	std::ostream& operator<<(
			std::ostream &o_ostream,
			const PlaneData_PhysicalComplex &i_dataArray
	)
	{
		PlaneData_PhysicalComplex &rw_array_data = (PlaneData_PhysicalComplex&)i_dataArray;

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


	void physical_update_lambda_array_indices(
			std::function<void(int,int,std::complex<double>&)> i_lambda,	///< lambda function to return value for lat/mu
			bool i_anti_aliasing = true
	)
	{
		PLANE_DATA_COMPLEX_PHYSICAL_FOR_2D_IDX(
				i_lambda(i, j, physical_space_data[idx])
		);
	}


	void physical_update_lambda_array_idx(
			std::function<void(int,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR

		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			i_lambda(i, physical_space_data[i]);
		}
	}


	void physical_update_lambda_unit_coordinates_corner_centered(
			std::function<void(double,double,std::complex<double>&)> i_lambda,	///< lambda function to return value for lat/mu
			bool i_anti_aliasing = true
	)
	{
		PLANE_DATA_COMPLEX_PHYSICAL_FOR_2D_IDX(
				i_lambda(
						(double)i/(double)planeDataConfig->physical_res[0],
						(double)j/(double)planeDataConfig->physical_res[1],
						physical_space_data[idx]
				)
		);
	}

	void physical_update_lambda_unit_coordinates_cell_centered(
			std::function<void(double,double,std::complex<double>&)> i_lambda,	///< lambda function to return value for lat/mu
			bool i_anti_aliasing = true
	)
	{
		PLANE_DATA_COMPLEX_PHYSICAL_FOR_2D_IDX(
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
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = 0;
	}


	/*
	 * Set all values to a specific value
	 */
	void physical_set_all_value(
			std::complex<double> &i_value
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_value;
	}

	/*
	 * Set all values to a specific value
	 */
	void physical_set_all_value(
			double i_value_real,
			double i_value_imag
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			physical_space_data[i].real(i_value_real);
			physical_space_data[i].imag(i_value_imag);
		}
	}


	/*
	 * Set all values to a specific value
	 */
	void physical_set_value(
			int i_y_idx,
			int i_x_idx,
			std::complex<double> &i_value
	)
	{
		physical_space_data[i_y_idx*planeDataConfig->physical_res[0] + i_x_idx] = i_value;
	}

	/*
	 * Set all values to a specific value
	 */
	void physical_set_value(
			int i_y_idx,
			int i_x_idx,
			double i_value_real,
			double i_value_imag
	)
	{
		physical_space_data[i_y_idx*planeDataConfig->physical_res[0] + i_x_idx].real(i_value_real);
		physical_space_data[i_y_idx*planeDataConfig->physical_res[0] + i_x_idx].imag(i_value_imag);
	}


	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_max(
			const PlaneData_PhysicalComplex &i_plane_data
	)
	{
		check(i_plane_data.planeDataConfig);

		double error = -1;

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						error,
						std::abs(
								physical_space_data[j] - i_plane_data.physical_space_data[j]
							)
						);
		}
		return error;
	}

	double physical_reduce_rms()
	{
		double error = 0;

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			std::complex<double> &d = physical_space_data[j];
			error += d.real()*d.real() + d.imag()*d.imag();
		}

		return std::sqrt(error / (double)planeDataConfig->physical_array_data_number_of_elements);
	}


	/**
	 * Return the maximum error norm
	 */
	double physical_reduce_max_abs()	const
	{
		double error = -1;

		for (std::size_t j = 0; j < planeDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						error,
						std::abs(physical_space_data[j].real())
						);

			error = std::max(
						error,
						std::abs(physical_space_data[j].imag())
						);
		}

		return error;
	}

	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double physical_reduce_sum_re_quad()	const
	{
		double sum = 0;
		double c = 0;

#if SWEET_THREADING_SPACE
//#if !SWEET_THREADING_TIME_REXI
		#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			double value = physical_space_data[i].real();

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
	 * reduce to root mean square
	 */
	double physical_reduce_rms_quad()
	{
		double sum = 0;
		double c = 0;

#if SWEET_THREADING_SPACE
//#if !SWEET_THREADING_TIME_REXI
		#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			double radius2 = physical_space_data[i].real()*physical_space_data[i].real()+physical_space_data[i].imag()*physical_space_data[i].imag();
			//double value = std::sqrt(radius2);
			//value *= value;
			double value = radius2;

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		sum = std::sqrt(sum/double(planeDataConfig->physical_array_data_number_of_elements));
		return sum;
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
			double value = physical_space_data[i].real()*physical_space_data[i].real() + physical_space_data[i].imag()*physical_space_data[i].imag();

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return std::sqrt(sum);
	}


};

#include "PlaneData_PhysicalComplex_Operators.hpp"

}	// namespace sweet

//#include "PlaneData_PhysicalComplex_Operators.hpp"

#endif
