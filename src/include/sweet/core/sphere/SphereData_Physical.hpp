/*
 * SphereDataPhysical.hpp
 *
 *  Created on: 01 Jan 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SPHERE_DATA_PHYSICAL_HPP_
#define SPHERE_DATA_PHYSICAL_HPP_

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
#include <vector>
#include <iterator>

#include <cmath>
#include <sweet/core/MemBlockAlloc.hpp>
#include <sweet/core/openmp_helper.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/SWEETError.hpp>



class SphereData_Physical
{
	friend class SphereData_SpectralComplex;

public:
	const SphereData_Config *sphereDataConfig;

public:
	double *physical_space_data;


	void swap(
			SphereData_Physical &i_sphereData
	)
	{
		assert(sphereDataConfig == i_sphereData.sphereDataConfig);

		std::swap(physical_space_data, i_sphereData.physical_space_data);
	}

public:
	SphereData_Physical(
			const SphereData_Config *i_sphereDataConfig
	)	:
		/// important: set this to nullptr, since a check for this will be performed by setup(...)
		sphereDataConfig(i_sphereDataConfig),
		physical_space_data(nullptr)
	{
		alloc_data();
	}


public:
	SphereData_Physical(
			const SphereData_Config *i_sphereDataConfig,
			double i_value
	)	:
		/// important: set this to nullptr, since a check for this will be performed by setup(...)
		sphereDataConfig(i_sphereDataConfig),
		physical_space_data(nullptr)
	{
		alloc_data();
		physical_set_all_value(i_value);
	}


public:
	SphereData_Physical()	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr)
	{
	}


public:
	SphereData_Physical(
			const SphereData_Physical &i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		physical_space_data(nullptr)

	{
		if (i_sph_data.sphereDataConfig != nullptr)
			alloc_data();

		operator=(i_sph_data);
	}


public:
	SphereData_Physical(
			SphereData_Physical &&i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		physical_space_data(nullptr)
	{
		if (i_sph_data.sphereDataConfig == nullptr)
			return;

		std::swap(physical_space_data, i_sph_data.physical_space_data);
	}



	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(
			const SphereData_Config *i_sphereDataConfig
	)	const
	{
		assert(sphereDataConfig->physical_num_lat == i_sphereDataConfig->physical_num_lat);
		assert(sphereDataConfig->physical_num_lon == i_sphereDataConfig->physical_num_lon);
	}



public:
	SphereData_Physical& operator=(
			const SphereData_Physical &i_sph_data
	)
	{
		if (i_sph_data.sphereDataConfig == nullptr)
			return *this;

		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		memcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		return *this;
	}


public:
	SphereData_Physical& operator=(
			SphereData_Physical &&i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		std::swap(physical_space_data, i_sph_data.physical_space_data);

		return *this;
	}


	SphereData_Physical operator+(
			const SphereData_Physical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] + i_sph_data.physical_space_data[idx];

		return out;
	}



	SphereData_Physical& operator+=(
			const SphereData_Physical &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_sph_data.physical_space_data[idx];

		return *this;
	}


	SphereData_Physical& operator+=(
			double i_scalar
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_scalar;

		return *this;
	}


	SphereData_Physical& operator-=(
			const SphereData_Physical &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_sph_data.physical_space_data[idx];

		return *this;
	}


	SphereData_Physical& operator-=(
			double i_scalar
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_scalar;

		return *this;
	}



	SphereData_Physical operator-(
			const SphereData_Physical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] - i_sph_data.physical_space_data[idx];

		return out;
	}



	SphereData_Physical operator-()
	{
		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = -physical_space_data[idx];

		return out;
	}



	SphereData_Physical operator*(
			const SphereData_Physical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_sph_data.physical_space_data[i];

		return out;
	}



	SphereData_Physical operator/(
			const SphereData_Physical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		check(i_sph_data.sphereDataConfig);

		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]/i_sph_data.physical_space_data[i];

		return out;
	}



	SphereData_Physical operator*(
			const double i_value
	)	const
	{
		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_value;

		return out;
	}




	const SphereData_Physical& operator*=(
			const double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] *= i_value;

		return *this;
	}




	SphereData_Physical operator/(
			double i_value
	)	const
	{
		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]/i_value;

		return out;
	}



	SphereData_Physical operator+(
			double i_value
	)	const
	{
		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]+i_value;

		return out;
	}



	SphereData_Physical operator-(
			double i_value
	)	const
	{
		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]-i_value;

		return out;
	}


	SphereData_Physical operator_scalar_sub_this(
			double i_value
	)	const
	{
		SphereData_Physical out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = i_value - physical_space_data[idx];

		return out;
	}



public:
	void setup(
			const SphereData_Config *i_sphereDataConfig
	)
	{
		if (sphereDataConfig != nullptr)
			SWEETError("Setup called twice!");

		sphereDataConfig = i_sphereDataConfig;
		alloc_data();
	}



private:
	void alloc_data()
	{
		assert(physical_space_data == nullptr);
		physical_space_data = MemBlockAlloc::alloc<double>(sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));
	}




public:
	void setup_if_required(
			const SphereData_Config *i_sphereDataConfig
	)
	{
		if (sphereDataConfig != nullptr)
			if (sphereDataConfig == i_sphereDataConfig)
				return;

		setup(i_sphereDataConfig);
	}



public:
	~SphereData_Physical()
	{
		free();
	}


public:
	void free()
	{
		if (physical_space_data != nullptr)
		{
			MemBlockAlloc::free(physical_space_data, sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));
			physical_space_data = nullptr;
		}

		sphereDataConfig = nullptr;
	}





	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */
	void physical_update_lambda(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				//double colatitude = acos(shtns->ct[j]);

				/*
				 * Colatitude is 0 at the north pole and 180 at the south pole
				 *
				 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
				 */
				//double lat_degree = M_PI*0.5 - colatitude;
				double lat_degree = sphereDataConfig->lat[j];

				i_lambda(lon_degree, lat_degree, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

#else

		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			double lat_degree = sphereDataConfig->lat[jlat];

			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double lon_degree = ((double)ilon/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

				//double colatitude = acos(shtns->ct[j]);

				/*
				 * Colatitude is 0 at the north pole and 180 at the south pole
				 *
				 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
				 */
				//double lat_degree = M_PI*0.5 - colatitude;

				i_lambda(lon_degree, lat_degree, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}

#endif
	}


	void physical_update_lambda_array(
			std::function<void(int,int,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{

		SWEET_THREADING_SPACE_PARALLEL_FOR

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				i_lambda(i, j, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

#else

		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				i_lambda(ilon, jlat, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}

#endif
	}


	void physical_update_lambda_array_idx(
			std::function<void(int,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR

		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
		{
			i_lambda(i, physical_space_data[i]);
		}
	}


	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude sin(phi) \in [-1;1])
	 */
	void physical_update_lambda_gaussian_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				double sin_phi = sphereDataConfig->lat_gaussian[j];

				i_lambda(lon_degree, sin_phi, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}
#else

		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			double sin_phi = sphereDataConfig->lat_gaussian[jlat];

			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double lon_degree = ((double)ilon/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

				i_lambda(lon_degree, sin_phi, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}
#endif
	}




	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters:
	 *   (longitude \in [0;2*pi], Cogaussian latitude cos(phi) \in [0;1])
	 */
	void physical_update_lambda_cogaussian_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = (((double)i)/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				double cos_phi = sphereDataConfig->lat_cogaussian[j];

				/*
				 * IDENTITAL FORMULATION
				double mu = shtns->ct[j];
				double comu = sqrt(1.0-mu*mu);
				*/

				i_lambda(lon_degree, cos_phi, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}
#else

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			double cos_phi = sphereDataConfig->lat_cogaussian[jlat];

			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double lon_degree = (((double)ilon)/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

				/*
				 * IDENTITAL FORMULATION
				double mu = shtns->ct[j];
				double comu = sqrt(1.0-mu*mu);
				*/

				i_lambda(lon_degree, cos_phi, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}
#endif
	}


	void physical_update_lambda_sinphi_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		physical_update_lambda_gaussian_grid(i_lambda);
	}

	void physical_update_lambda_cosphi_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		physical_update_lambda_cogaussian_grid(i_lambda);
	}



	/*
	 * Set all values to zero
	 */
	void physical_set_zero()
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
				physical_space_data[j*sphereDataConfig->physical_num_lon + i] = 0;
	}



	/*
	 * Set all values to a specific value
	 */
	void physical_set_all_value(
			double i_value
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
				physical_space_data[j*sphereDataConfig->physical_num_lon + i] = i_value;
	}



	/*
	 * Set all values to a specific value
	 */
	void physical_set_value(
			int i_lon_idx,
			int i_lat_idx,
			double i_value
	)
	{
#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
		physical_space_data[i_lon_idx*sphereDataConfig->physical_num_lat + i_lat_idx] = i_value;
#else
		physical_space_data[i_lat_idx*sphereDataConfig->physical_num_lon + i_lon_idx] = i_value;
#endif
	}



	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_max(
			const SphereData_Physical &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		double error = -1;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						std::abs(
								physical_space_data[j] - i_sph_data.physical_space_data[j]
							),
							error
						);
		}
		return error;
	}


	double physical_reduce_rms()
	{
		double error = 0;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			double &d = physical_space_data[j];
			error += d*d;
		}

		return std::sqrt(error / (double)sphereDataConfig->physical_array_data_number_of_elements);
	}

	double physical_reduce_norm1()
	{
		double error = 0;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			double &d = physical_space_data[j];
			error += std::abs(d);
		}

		return error;
	}

	double physical_reduce_norm2()
	{
		double error = 0;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			double &d = physical_space_data[j];
			error += d*d;
		}

		return std::sqrt(error);
	}


	double physical_reduce_sum()	const
	{
		double sum = 0;
		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
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
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
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
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
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
		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			sum += physical_space_data[j];
		}

		return sum;
	}




	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_max_abs(
			const SphereData_Physical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		double error = -1;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						std::abs(
								physical_space_data[j] - i_sph_data.physical_space_data[j]
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

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
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

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			error = std::min(physical_space_data[j], error);

		return error;
	}


	/**
	 * Return the minimum value
	 */
	double physical_reduce_max()	const
	{
		double error = -std::numeric_limits<double>::infinity();

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			error = std::max(physical_space_data[j], error);

		return error;
	}


	/**
	 * Rescale data so that max_abs returns the given value
	 */
	SphereData_Physical physical_rescale_to_max_abs(
			double i_new_max_abs
	)	const
	{
		double max_abs = physical_reduce_max_abs();
		double scale = i_new_max_abs/max_abs;

		SphereData_Physical out(sphereDataConfig);

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			out.physical_space_data[j] = physical_space_data[j]*scale;

		return out;
	}



	bool physical_isAnyNaNorInf()	const
	{
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
		{
			if (std::isnan(physical_space_data[i]) || std::isinf(physical_space_data[i]) != 0)
				return true;
		}

		return false;
	}






	void physical_print(
			int i_precision = -1
	)	const
	{
		if (i_precision >= 0)
			std::cout << std::setprecision(i_precision);

        for (int j = (int)(sphereDataConfig->physical_num_lat-1); j >= 0; j--)
        {
        		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
        		{
#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
        			std::cout << physical_space_data[i*sphereDataConfig->physical_num_lat+j];
#else
        			std::cout << physical_space_data[j*sphereDataConfig->physical_num_lon+i];
#endif
        			if (i < sphereDataConfig->physical_num_lon-1)
        				std::cout << "\t";
        		}
        		std::cout << std::endl;
        }
	}



	void physical_file_write(
			const std::string &i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const
	{
		std::ofstream file(i_filename, std::ios_base::trunc);

		if (i_precision >= 0)
			file << std::setprecision(i_precision);

		file << "#TI " << i_title << std::endl;
		file << "#TX Longitude" << std::endl;
		file << "#TY Latitude" << std::endl;

		//file << "lat\\lon\t";
		// Use 0 to make it processable by python
		file << "0\t";

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
//			double lon_degree = ((double)i/(double)sphereDataConfig->spat_num_lon)*2.0*M_PI;
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;
			lon_degree = lon_degree/M_PI*180.0;

			file << lon_degree;
			if (i < sphereDataConfig->physical_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
        {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
        		double lat_degree = sphereDataConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		file << lat_degree << "\t";

        		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
        		{
#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
        			file << physical_space_data[i*sphereDataConfig->physical_num_lat+j];
#else
        			file << physical_space_data[j*sphereDataConfig->physical_num_lon+i];
#endif
        			if (i < sphereDataConfig->physical_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}



	void physical_file_write_lon_pi_shifted(
			const char *i_filename,
			std::string i_title = "",
			int i_precision = 20
	)	const
	{
		std::ofstream file(i_filename, std::ios_base::trunc);

		file << std::setprecision(i_precision);
		file << "#TI " << i_title << std::endl;
		file << "#TX Longitude" << std::endl;
		file << "#TY Latitude" << std::endl;

		//file << "lat\\lon\t";
		// Use 0 to make it processable by python
		file << "0\t";

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
//			double lon_degree = ((double)i/(double)sphereDataConfig->spat_num_lon)*2.0*M_PI;
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;
			lon_degree = (lon_degree-M_PI)/M_PI*180.0;

			file << lon_degree;
			if (i < sphereDataConfig->physical_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
        {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
        		double lat_degree = sphereDataConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		file << lat_degree << "\t";

        		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
        		{
        			int ia = i+sphereDataConfig->physical_num_lon/2;
        			if (ia >= sphereDataConfig->physical_num_lon)
        				ia -= sphereDataConfig->physical_num_lon;

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
        			file << physical_space_data[ia*sphereDataConfig->physical_num_lat+j];
#else
        			file << physical_space_data[j*sphereDataConfig->physical_num_lon+ia];
#endif
        			if (i < sphereDataConfig->physical_num_lon-1)
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


			std::size_t expected_size = sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements;

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
		while (row < sphereDataConfig->physical_num_lat)
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
					int lon = col;
					int lat = sphereDataConfig->physical_num_lat-row-1;

					if (lon > sphereDataConfig->physical_num_lon)
						return false;

					if (lat > sphereDataConfig->physical_num_lat)
						return false;

					physical_set_value(lon, lat, i_value);

					col++;
				}

				last_pos = pos+1;
		    }

			if (col != sphereDataConfig->physical_num_lon)
			{
				std::cerr << "ERROR: column mismatch - Failed to read data from file " << i_filename << " in line " << row << ", column " << col << std::endl;
				return false;
			}

			row++;
		}

		if (row != sphereDataConfig->physical_num_lat)
		{
			std::cerr << "ERROR: rows mismatch - Failed to read data from file " << i_filename << " in line " << row << std::endl;
			return false;
		}


		return true;
	}


#if SWEET_PARAREAL || SWEET_XBRAID

	/**
	 * Load data from ASCII file.
	 *
	 * Read csv files from reference data in parareal, in order to compute parareal erros online
	 * instead of producing csv files during the parareal simulation.
	 *
	 * \return true if data was successfully read
	 */
	bool file_physical_loadRefData_Parareal(
			const char *i_filename,		///< Name of file to load data from
			bool i_binary_data = false	///< load as binary data (disabled per default)
	)
	{

		///SphereDataConfig sphereDataConfig_ref;

		std::cout << "loading DATA from " << i_filename << std::endl;

		std::ifstream file(i_filename);

		for (int i = 0; i < 4; i++)
		{
			std::string line;
			std::getline(file, line);
			std::istringstream iss(line);
			std::vector<std::string> str_vector((std::istream_iterator<std::string>(iss)),
				std::istream_iterator<std::string>());

			if (i == 0)
			{
				assert(str_vector.size() == 1);
				assert(str_vector[0] == "#TI");
			}
			else if (i == 1)
			{
				assert(str_vector.size() == 2);
				assert(str_vector[0] == "#TX");
				assert(str_vector[1] == "Longitude");
			}
			else if (i == 2)
			{
				assert(str_vector.size() == 2);
				assert(str_vector[0] == "#TY");
				assert(str_vector[1] == "Latitude");
			}
			// Line 4: longitude values
			// First character is "0"
			else if (i == 3)
			{
				assert((int)str_vector.size() == (int)sphereDataConfig->physical_num_lon + 1);
				assert(str_vector[0] == "0");
				for (int l = 0; l < sphereDataConfig->physical_num_lon; l++)
					assert(std::abs(atof(str_vector[l + 1].c_str()) - ((double)l/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI/M_PI*180.0  ) < 1e-13);
			}
		}


		///sphereDataConfig_ref.setup(resx_ref, resy_ref, (int)((resx_ref * 2) / 3), (int)((resx_ref * 2) / 3), planeDataConfig->reuse_spectral_transformation_plans);
		///*this = SphereData_Physical(&sphereDataConfig_ref);


		for (int j = sphereDataConfig->physical_num_lat - 1; j >= 0; j--)
		{
			std::string line;
			std::getline(file, line);
			std::istringstream iss(line);
			std::vector<std::string> str_vector((std::istream_iterator<std::string>(iss)),
				std::istream_iterator<std::string>());

			if (!file.good())
			{
				std::cerr << "Failed to read data from file " << i_filename << " in line " << sphereDataConfig->physical_num_lat - 1 - j + 3 << std::endl;
				return false;
			}


			assert((int)str_vector.size() == (int)sphereDataConfig->physical_num_lon + 1);
			assert(std::abs(atof(str_vector[0].c_str()) - sphereDataConfig->lat[j]/M_PI*180.0 ) < 1e-13);

			for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			{
				double val = atof(str_vector[i + 1].c_str());
#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
        			physical_space_data[i*sphereDataConfig->physical_num_lat+j] = val;
#else
        			physical_space_data[j*sphereDataConfig->physical_num_lon+i] = val;
#endif
			}
		}

		file.close();
		std::cout << "DATA loaded OK" << std::endl;

		return true;
	}
#endif












	void file_write_raw(
			const std::string &i_filename
	)	const
	{
		std::fstream file(i_filename, std::ios::out | std::ios::binary);
		file.write((const char*)physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);
	}



	void file_read_raw(
			const std::string &i_filename
	)	const
	{
		std::fstream file(i_filename, std::ios::in | std::ios::binary);
		file.read((char*)physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);
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
		for (std::size_t idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
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
SphereData_Physical operator*(
		const double i_value,
		const SphereData_Physical &i_array_data
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
SphereData_Physical operator-(
		const double i_value,
		const SphereData_Physical &i_array_data
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
SphereData_Physical operator+(
		const double i_value,
		const SphereData_Physical &i_array_data
)
{
	return i_array_data+i_value;
}




#endif
