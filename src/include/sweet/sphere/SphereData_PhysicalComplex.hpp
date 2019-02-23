/*
 * SphereDataPhysicalComplex.hpp
 *
 *  Created on: 3 Mar 2017
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SPHERE_DATA_PHYSICAL_COMPLEX_HPP_
#define SPHERE_DATA_PHYSICAL_COMPLEX_HPP_

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
#include <sweet/FatalError.hpp>
#include <sweet/openmp_helper.hpp>
#include <sweet/sphere/SphereData_Config.hpp>



class SphereData_PhysicalComplex
{
	friend class SphereData_SpectralComplex;

public:
	const SphereData_Config *sphereDataConfig;

public:
	std::complex<double> *physical_space_data;


	void swap(
			SphereData_PhysicalComplex &i_sphereData
	)
	{
		assert(sphereDataConfig == i_sphereData.sphereDataConfig);

		std::swap(physical_space_data, i_sphereData.physical_space_data);
	}

public:
	SphereData_PhysicalComplex(
			const SphereData_Config *i_sphereDataConfig
	)	:
		sphereDataConfig(i_sphereDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_sphereDataConfig);
	}


public:
	SphereData_PhysicalComplex()	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr)
	{
	}


public:
	SphereData_PhysicalComplex(
			const SphereData_PhysicalComplex &i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);

		operator=(i_sph_data);
	}


public:
	SphereData_PhysicalComplex(
			const SphereData_Physical &i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		physical_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);

		operator=(i_sph_data);
	}


	/*
	 * load real and imaginary data from physical arrays
	 */
	void loadRealImag(
			const SphereData_Physical &i_re,
			const SphereData_Physical &i_im
	)
	{
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
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
			const SphereData_Config *i_sphereDataConfig
	)	const
	{
		assert(sphereDataConfig->physical_num_lat == i_sphereDataConfig->physical_num_lat);
		assert(sphereDataConfig->physical_num_lon == i_sphereDataConfig->physical_num_lon);
	}



public:
	SphereData_PhysicalComplex& operator=(
			const SphereData_PhysicalComplex &i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		memcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(std::complex<double>)*sphereDataConfig->physical_array_data_number_of_elements);

		return *this;
	}


public:
	SphereData_PhysicalComplex& operator=(
			SphereData_PhysicalComplex &&i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		std::swap(physical_space_data, i_sph_data.physical_space_data);

		return *this;
	}



public:
	SphereData_PhysicalComplex& operator=(
			const SphereData_Physical &i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_sph_data.physical_space_data[i];

		return *this;
	}

#if 0
	SphereData_PhysicalComplex robert_convertToRobert()
	{
		SphereData_PhysicalComplex out = *this;

		out.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					io_data *= std::cos(i_lat);
				}
		);

		return out;
	}



	SphereData_PhysicalComplex robert_convertToNonRobert()
	{
		SphereData_PhysicalComplex out = *this;

		out.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					io_data /= std::cos(i_lat);
				}
		);

		return out;
	}
#endif

	SphereData_PhysicalComplex operator+(
			const SphereData_PhysicalComplex &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] + i_sph_data.physical_space_data[idx];

		return out;
	}



	SphereData_PhysicalComplex& operator+=(
			const SphereData_PhysicalComplex &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_sph_data.physical_space_data[idx];

		return *this;
	}


	SphereData_PhysicalComplex& operator-=(
			const SphereData_PhysicalComplex &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_sph_data.physical_space_data[idx];

		return *this;
	}



	SphereData_PhysicalComplex operator-(
			const SphereData_PhysicalComplex &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx] - i_sph_data.physical_space_data[idx];

		return out;
	}



	SphereData_PhysicalComplex operator-()
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = -physical_space_data[idx];

		return out;
	}



	SphereData_PhysicalComplex operator*(
			const SphereData_PhysicalComplex &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_sph_data.physical_space_data[i];

		return out;
	}



	SphereData_PhysicalComplex operator/(
			const SphereData_PhysicalComplex &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]/i_sph_data.physical_space_data[i];

		return out;
	}



	SphereData_PhysicalComplex operator*(
			const double i_value
	)	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = physical_space_data[i]*i_value;

		return out;
	}



	SphereData_PhysicalComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]*i_value;

		return out;
	}




	const SphereData_PhysicalComplex& operator*=(
			const double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] *= i_value;

		return *this;
	}



	SphereData_PhysicalComplex operator/(
			double i_value
	)	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]/i_value;

		return out;
	}



	SphereData_PhysicalComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]+i_value;

		return out;
	}



	SphereData_PhysicalComplex operator-(
			const std::complex<double> &i_value
	)	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			out.physical_space_data[idx] = physical_space_data[idx]-i_value;

		return out;
	}




public:
	void setup(
			const SphereData_Config *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		physical_space_data = MemBlockAlloc::alloc<std::complex<double>>(sphereDataConfig->physical_array_data_number_of_elements * sizeof(std::complex<double>));
	}



public:
	~SphereData_PhysicalComplex()
	{
		if (physical_space_data != nullptr)
			MemBlockAlloc::free(physical_space_data, sphereDataConfig->physical_array_data_number_of_elements * sizeof(std::complex<double>));
	}





	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */
	void physical_update_lambda(
			std::function<void(double,double,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
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
			std::function<void(int,int,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{


#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				i_lambda(i, j, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

#else

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
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
			std::function<void(double,double,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR
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

		SWEET_THREADING_SPACE_PARALLEL_FOR
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
			std::function<void(double,double,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
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
			std::function<void(double,double,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		physical_update_lambda_gaussian_grid(i_lambda);
	}

	void physical_update_lambda_cosphi_grid(
			std::function<void(double,double,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		physical_update_lambda_cogaussian_grid(i_lambda);
	}



	/*
	 * Set all values to zero
	 */
	void physical_set_zero()
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
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
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_value;
	}



	/*
	 * Set all values to a specific value
	 */
	void physical_set_value(
			int i_lon_idx,
			int i_lat_idx,
			std::complex<double> &i_value
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
			const SphereData_PhysicalComplex &i_sph_data
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

#if 0
	double physical_reduce_rms()
	{
		double error = 0;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			std::complex<double> &d = physical_space_data[j];
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
#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
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
#if SWEET_THREADING_SPACE
#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
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
#endif




	/**
	 * Return the maximum error norm
	 */
	double physical_reduce_max_abs()	const
	{
		double error = -1;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
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


#if 0

	/**
	 * Return the minimum value
	 */
	double physical_reduce_min()	const
	{
		double error = std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::min(
						error,
						std::abs(physical_space_data[j].real())
						);

			error = std::min(
						error,
						std::abs(physical_space_data[j].imag())
						);
		}


		return error;
	}


	/**
	 * Rescale data so that max_abs returns the given value
	 */
	SphereData_PhysicalComplex physical_rescale_to_max_abs(
			double i_new_max_abs
	)	const
	{
		double max_abs = physical_reduce_max_abs();
		double scale = i_new_max_abs/max_abs;

		SphereData_PhysicalComplex out(sphereDataConfig);

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			out.physical_space_data[j] = physical_space_data[j]*scale;

		return out;
	}


	void physical_print(
			int i_precision = 20
	)	const
	{
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
#endif
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
SphereData_PhysicalComplex operator*(
		const double i_value,
		const SphereData_PhysicalComplex &i_array_data
)
{
	return ((SphereData_PhysicalComplex&)i_array_data)*i_value;
}


inline
static
SphereData_PhysicalComplex operator*(
		const std::complex<double> &i_value,
		const SphereData_PhysicalComplex &i_array_data
)
{
	return ((SphereData_PhysicalComplex&)i_array_data)*i_value;
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
SphereData_PhysicalComplex operator-(
		const double i_value,
		const SphereData_PhysicalComplex &i_array_data
)
{
	return ((SphereData_PhysicalComplex&)i_array_data).valueMinusThis(i_value);
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
SphereData_PhysicalComplex operator+(
		const double i_value,
		const SphereData_PhysicalComplex &i_array_data
)
{
	i_array_data.checkConsistency();
	return ((SphereData_PhysicalComplex&)i_array_data)+i_value;
}
#endif

#endif /* SPHDATA_HPP_ */
