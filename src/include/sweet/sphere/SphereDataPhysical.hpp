/*
 * SphereDataPhysical.hpp
 *
 *  Created on: 01 Jan 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com> Schreiber <M.Schreiber@exeter.ac.uk>
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

#include <sweet/sweetmath.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/FatalError.hpp>
#include <sweet/openmp_helper.hpp>



class SphereDataPhysical
{
	friend class SphereDataComplex;

public:
	const SphereDataConfig *sphereDataConfig;

public:
	double *physical_space_data;


	void swap(
			SphereDataPhysical &i_sphereData
	)
	{
		assert(sphereDataConfig == i_sphereData.sphereDataConfig);

		std::swap(physical_space_data, i_sphereData.physical_space_data);
	}

public:
	SphereDataPhysical(
			const SphereDataConfig *i_sphereDataConfig
	)	:
		sphereDataConfig(i_sphereDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_sphereDataConfig);
	}


public:
	SphereDataPhysical()	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr)
	{
	}


public:
	SphereDataPhysical(
			const SphereDataPhysical &i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);

		operator=(i_sph_data);
	}



	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(
			const SphereDataConfig *i_sphereDataConfig
	)	const
	{
		assert(sphereDataConfig->physical_num_lat == i_sphereDataConfig->physical_num_lat);
		assert(sphereDataConfig->physical_num_lon == i_sphereDataConfig->physical_num_lon);
	}



public:
	SphereDataPhysical& operator=(
			const SphereDataPhysical &i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		memcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		return *this;
	}


public:
	SphereDataPhysical& operator=(
			SphereDataPhysical &&i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		std::swap(physical_space_data, i_sph_data.physical_space_data);

		return *this;
	}


	SphereDataPhysical robert_convertToRobert()
	{
		SphereDataPhysical out = *this;

		out.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					io_data *= std::cos(i_lat);
				}
		);

		return out;
	}



	SphereDataPhysical robert_convertToNonRobert()
	{
		SphereDataPhysical out = *this;

		out.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					io_data /= std::cos(i_lat);
				}
		);

		return out;
	}


	SphereDataPhysical operator+(
			const SphereDataPhysical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] + i_sph_data.physical_space_data[idx];

		return out;
	}



	SphereDataPhysical& operator+=(
			const SphereDataPhysical &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_sph_data.physical_space_data[idx];

		return *this;
	}


	SphereDataPhysical& operator-=(
			const SphereDataPhysical &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_sph_data.physical_space_data[idx];

		return *this;
	}



	SphereDataPhysical operator-(
			const SphereDataPhysical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] - i_sph_data.physical_space_data[idx];

		return out;
	}



	SphereDataPhysical operator-()
	{
		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = -physical_space_data[idx];

		return out;
	}



	SphereDataPhysical operator*(
			const SphereDataPhysical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_sph_data.physical_space_data[i];

		return out;
	}



	SphereDataPhysical operator/(
			const SphereDataPhysical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		check(i_sph_data.sphereDataConfig);

		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]/i_sph_data.physical_space_data[i];

		return out;
	}



	SphereDataPhysical operator*(
			const double i_value
	)	const
	{
		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_value;

		return out;
	}




	const SphereDataPhysical& operator*=(
			const double i_value
	)	const
	{
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] *= i_value;

		return *this;
	}





	const SphereDataPhysical& operator*=(
			const double &i_value
	)	const
	{
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] *= i_value;

		return *this;
	}



	SphereDataPhysical operator/(
			double i_value
	)	const
	{
		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]/i_value;

		return out;
	}



	SphereDataPhysical operator+(
			double i_value
	)	const
	{
		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]+i_value;

		return out;
	}



	SphereDataPhysical operator-(
			double i_value
	)	const
	{
		SphereDataPhysical out(sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]-i_value;

		return out;
	}




public:
	void setup(
			const SphereDataConfig *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		physical_space_data = MemBlockAlloc::alloc<double>(sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));
	}



public:
	~SphereDataPhysical()
	{
		if (physical_space_data != nullptr)
			MemBlockAlloc::free(physical_space_data, sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));
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
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif

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

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif

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


	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude sin(phi) \in [-1;1])
	 */
	void physical_update_lambda_gaussian_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif

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

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
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
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif

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
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif

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
			const SphereDataPhysical &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		double error = -1;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						error,
						std::abs(
								physical_space_data[j] - i_sph_data.physical_space_data[j]
							)
						);
		}
		return error;
	}


	double physical_reduce_rms()
	{
		double error = 0;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			double &d = physical_space_data[j];
			error += d*d;
		}

		return std::sqrt(error / (double)sphereDataConfig->physical_array_data_number_of_elements);
	}



	double physical_reduce_sum()
	{
		double sum = 0;
		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
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
#if SWEET_SPACE_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
#endif
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
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
#if SWEET_SPACE_THREADING
#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
#endif
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
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
		FatalError("TODO: Implement metric-scaled summation");
		double sum = 0;
		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			sum += physical_space_data[j];
		}

		return sum;
	}




	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_max_abs(
			const SphereDataPhysical &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		double error = -1;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						error,
						std::abs(
								physical_space_data[j] - i_sph_data.physical_space_data[j]
							)
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

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						error,
						std::abs(physical_space_data[j])
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

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			error = std::min(error, physical_space_data[j]);

		return error;
	}


	/**
	 * Return the minimum value
	 */
	double physical_reduce_max()	const
	{
		double error = -std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			error = std::max(error, physical_space_data[j]);

		return error;
	}


	/**
	 * Rescale data so that max_abs returns the given value
	 */
	SphereDataPhysical physical_rescale_to_max_abs(
			double i_new_max_abs
	)	const
	{
		double max_abs = physical_reduce_max_abs();
		double scale = i_new_max_abs/max_abs;

		SphereDataPhysical out(sphereDataConfig);

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			out.physical_space_data[j] = physical_space_data[j]*scale;

		return out;
	}


	void physical_print(
			int i_precision = -1
	)	const
	{
		if (i_precision >= 0)
			std::cout << std::setprecision(i_precision);

        for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
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
				FatalError(std::string("Failed to open file ")+i_filename);

			file.seekg(0, std::ios::end);
			std::size_t size = file.tellg();
			file.seekg(0, std::ios::beg);


			std::size_t expected_size = sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements;

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
SphereDataPhysical operator*(
		const double i_value,
		const SphereDataPhysical &i_array_data
)
{
	return ((SphereDataPhysical&)i_array_data)*i_value;
}


/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 *
 * Otherwise, we'd have to write it as arrayData-1.5
 *
 */
#if 0
inline
static
SphereDataPhysical operator-(
		const double i_value,
		const SphereDataPhysical &i_array_data
)
{
//	return ((SphereDataPhysical&)i_array_data).operator-(i_value);
//	return -(((SPHData&)i_array_data).operator-(i_value));
}
#endif
/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */

#if 0
inline
static
SphereDataPhysical operator+(
		const double i_value,
		const SphereDataPhysical &i_array_data
)
{
	i_array_data.checkConsistency();
	return ((SphereDataPhysical&)i_array_data)+i_value;
}
#endif

#endif /* SPHDATA_HPP_ */
