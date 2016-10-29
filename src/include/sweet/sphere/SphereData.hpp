/*
 * SphereHelper.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHERE_HELPER_HPP_
#define SPHERE_HELPER_HPP_

#include <complex>
#include <functional>
#include <array>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <limits>

#include <sweet/sweetmath.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>

#include <sweet/openmp_helper.hpp>


class SphereData
{
	friend class SphereDataComplex;

public:
	SphereDataConfig *sphereDataConfig;

public:
	double *physical_space_data;
	std::complex<double> *spectral_space_data;

	bool spectral_space_data_valid;
	bool physical_space_data_valid;

public:
	SphereData(
			SphereDataConfig *i_sphConfig
	)	:
		sphereDataConfig(i_sphConfig),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		assert(i_sphConfig != 0);

		setup(i_sphConfig);
	}


public:
	SphereData()	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
	}


public:
	SphereData(
			const SphereData &i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);

		operator=(i_sph_data);
	}



	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(SphereDataConfig *i_sphConfig)	const
	{
		assert(sphereDataConfig->physical_num_lat == i_sphConfig->physical_num_lat);
		assert(sphereDataConfig->physical_num_lon == i_sphConfig->physical_num_lon);

		assert(sphereDataConfig->spectral_modes_m_max == i_sphConfig->spectral_modes_m_max);
		assert(sphereDataConfig->spectral_modes_n_max == i_sphConfig->spectral_modes_n_max);
	}



public:
	SphereData& operator=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		if (i_sph_data.physical_space_data_valid)
			memcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		if (i_sph_data.spectral_space_data_valid)
			memcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(cplx)*sphereDataConfig->spectral_array_data_number_of_elements);

		physical_space_data_valid = i_sph_data.physical_space_data_valid;
		spectral_space_data_valid = i_sph_data.spectral_space_data_valid;

		return *this;
	}



public:
	const SphereData& spectral_copyToDifferentModes(
			SphereData &o_sph_data
	)	const
	{
		request_data_spectral();

		assert (o_sph_data.sphereDataConfig != sphereDataConfig);

		/*
		 *  0 = invalid
		 * -1 = scale down
		 *  1 = scale up
		 */
		int scaling_mode = 0;

		if (sphereDataConfig->spectral_modes_m_max < o_sph_data.sphereDataConfig->spectral_modes_m_max)
		{
			scaling_mode = 1;
		}
		else if (sphereDataConfig->spectral_modes_m_max > o_sph_data.sphereDataConfig->spectral_modes_m_max)
		{
			scaling_mode = -1;
		}
//		assert(scaling_mode != 0);


		if (sphereDataConfig->spectral_modes_n_max < o_sph_data.sphereDataConfig->spectral_modes_n_max)
		{
			assert(scaling_mode != -1);
			scaling_mode = 1;
		}
		else if (sphereDataConfig->spectral_modes_n_max > o_sph_data.sphereDataConfig->spectral_modes_n_max)
		{
			assert(scaling_mode != 1);
			scaling_mode = -1;
		}
		assert(scaling_mode != 0);


		if (scaling_mode == -1)
		{
			/*
			 * more modes -> less modes
			 */

#if SWEET_THREADING
#pragma omp parallel for
#endif
			for (int m = 0; m <= o_sph_data.sphereDataConfig->spectral_modes_m_max; m++)
			{
				cplx *dst = &o_sph_data.spectral_space_data[o_sph_data.sphereDataConfig->getArrayIndexByModes(m, m)];
				cplx *src = &spectral_space_data[sphereDataConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(cplx)*(o_sph_data.sphereDataConfig->spectral_modes_n_max-m+1);
				memcpy(dst, src, size);
			}
		}
		else
		{
			/*
			 * less modes -> more modes
			 */

			// zero all values
			o_sph_data.spectral_set_zero();

#if SWEET_THREADING
#pragma omp parallel for
#endif
			for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
			{
				cplx *dst = &o_sph_data.spectral_space_data[o_sph_data.sphereDataConfig->getArrayIndexByModes(m, m)];
				cplx *src = &spectral_space_data[sphereDataConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(cplx)*(sphereDataConfig->spectral_modes_n_max-m+1);
				memcpy(dst, src, size);
			}
		}

		o_sph_data.physical_space_data_valid = false;
		o_sph_data.spectral_space_data_valid = true;

		return *this;
	}


	void request_data_spectral()	const
	{
		if (spectral_space_data_valid)
			return;

		assert(physical_space_data_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		spat_to_SH(sphereDataConfig->shtns, physical_space_data, spectral_space_data);

		SphereData *this_var = (SphereData*)this;

		this_var->physical_space_data_valid = false;
		this_var->spectral_space_data_valid = true;
	}


	void request_data_physical()	const
	{
		if (physical_space_data_valid)
			return;

		assert(spectral_space_data_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		SH_to_spat(sphereDataConfig->shtns, spectral_space_data, physical_space_data);

		SphereData *this_var = (SphereData*)this;

		this_var->physical_space_data_valid = true;
		this_var->spectral_space_data_valid = false;
	}



	SphereData operator+(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData out_sph_data(sphereDataConfig);


#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	SphereData& operator+=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] += i_sph_data.spectral_space_data[idx];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}


	SphereData& operator-=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}



	SphereData operator-(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereData operator-()
	{
		request_data_spectral();

		SphereData out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = -spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereData operator*(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SphereData out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out_sph_data.physical_space_data[i] = i_sph_data.physical_space_data[i]*physical_space_data[i];

		out_sph_data.spectral_space_data_valid = false;
		out_sph_data.physical_space_data_valid = true;

		return out_sph_data;
	}



	SphereData operator*(
			const double i_value
	)	const
	{
		request_data_spectral();

		SphereData out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}




	const SphereData& operator*=(
			const double i_value
	)	const
	{
		request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}





	const SphereData& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}



	SphereData operator/(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereData operator+(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out_sph_data(*this);

		out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}




public:
	void setup(SphereDataConfig *i_sphConfig)
	{
		sphereDataConfig = i_sphConfig;

		spectral_space_data_valid = false;
		physical_space_data_valid = false;

		physical_space_data = MemBlockAlloc::alloc<double>(sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));
		spectral_space_data = MemBlockAlloc::alloc<cplx>(sphereDataConfig->spectral_array_data_number_of_elements * sizeof(cplx));
	}


public:
	~SphereData()
	{
		MemBlockAlloc::free(physical_space_data, sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));
		MemBlockAlloc::free(spectral_space_data, sphereDataConfig->spectral_array_data_number_of_elements * sizeof(cplx));
	}


	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	inline
	SphereData spectral_solve_helmholtz(
			const double &i_a,
			const double &i_b,
			double r
	)
	{
		SphereData out(*this);

		const double a = i_a;
		const double b = i_b/(r*r);

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				io_data /= (a + (-b*(double)n*((double)n+1.0)));
			}
		);

		return out;
	}


public:
	/**
	 * Truncate modes which are not representable in spectral space
	 */
	SphereData physical_truncate()
	{
		request_data_physical();

		SphereData out_sph_data(sphereDataConfig);
		spat_to_SH(sphereDataConfig->shtns, physical_space_data, out_sph_data.spectral_space_data);
		SH_to_spat(sphereDataConfig->shtns, out_sph_data.spectral_space_data, out_sph_data.physical_space_data);

		out_sph_data.physical_space_data_valid = true;
		out_sph_data.spectral_space_data_valid = false;

		return out_sph_data;
	}

	/**
	 * Truncate modes which are not representable in spectral space
	 */
	SphereData spectral_truncate()
	{
		request_data_spectral();

		SphereData out_sph_data(sphereDataConfig);
		SH_to_spat(sphereDataConfig->shtns, spectral_space_data, out_sph_data.physical_space_data);
		spat_to_SH(sphereDataConfig->shtns, out_sph_data.physical_space_data, out_sph_data.spectral_space_data);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	void spectral_update_lambda(
			std::function<void(int,int,cplx&)> i_lambda
	)
	{
		if (physical_space_data_valid)
			request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				i_lambda(n, m, spectral_space_data[idx]);
				idx++;
			}
		}

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
	}



	bool physical_isAnyNaNorInf()
	{
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
		{
			if (std::isnan(physical_space_data[i]) || std::isinf(physical_space_data[i]) != 0)
				return true;
		}

		return false;
	}


	const std::complex<double>& spectral_get(
			int in,
			int im
	)	const
	{
		static const std::complex<double> zero = {0,0};

		assert(spectral_space_data_valid);

		if (in < 0 ||  im < 0)
			return zero;

		if (in > sphereDataConfig->spectral_modes_n_max)
			return zero;

		if (im > sphereDataConfig->spectral_modes_m_max)
			return zero;

		if (im > in)
			return zero;

		assert (im <= sphereDataConfig->spectral_modes_m_max);
		return spectral_space_data[sphereDataConfig->getArrayIndexByModes(in, im)];
	}


	/*
	 * Set all values to zero
	 */
	void spectral_set_zero()
	{
#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int i = 0; i < sphereDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] = {0,0};

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
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
		if (spectral_space_data_valid)
			request_data_physical();

#if SWEET_THREADING
#pragma omp parallel for
#endif

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

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
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
		if (spectral_space_data_valid)
			request_data_physical();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				double sin_phi = sphereDataConfig->lat_gaussian[j];

				i_lambda(lon_degree, sin_phi, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

		spectral_space_data_valid = false;
		physical_space_data_valid = true;
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
		if (spectral_space_data_valid)
			request_data_physical();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = (((double)i)/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				double cos_phi = sphereDataConfig->lat_cogaussian[j];

//				std::cout << cos_phi << std::endl;
				/*
				 * IDENTITAL FORMULATION
				double mu = shtns->ct[j];
				double comu = sqrt(1.0-mu*mu);
				*/

				i_lambda(lon_degree, cos_phi, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

		spectral_space_data_valid = false;
		physical_space_data_valid = true;
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
#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
				physical_space_data[i*sphereDataConfig->physical_num_lat + j] = 0;

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}


	/*
	 * Set all values to a specific value
	 */
	void physical_set_value(
			double i_value
	)
	{
#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
				physical_space_data[i*sphereDataConfig->physical_num_lat + j] = i_value;

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}


	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_error_max(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

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
	double reduce_abs_max()
	{
		request_data_physical();

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
	 * Return the maximum value
	 */
	double reduce_min()
	{
		request_data_physical();

		double error = std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::min(error, physical_space_data[j]);
		}
		return error;
	}


	/**
	 * Return the minimum value
	 */
	double reduce_max()
	{
		request_data_physical();

		double error = -std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(error, physical_space_data[j]);
		}
		return error;
	}


	void spectral_print(
			int i_precision = 20
	)	const
	{
		request_data_spectral();

		std::cout << std::setprecision(i_precision);

		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				std::cout << spectral_space_data[idx] << "\t";
				idx++;
			}
			std::cout << std::endl;
		}
	}


	void physical_print(
			int i_precision = 20
	)	const
	{
		request_data_physical();

		std::cout << std::setprecision(i_precision);

        for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
        {
        		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
        		{
        			std::cout << physical_space_data[i*sphereDataConfig->physical_num_lat+j];
        			if (i < sphereDataConfig->physical_num_lon-1)
        				std::cout << "\t";
        		}
        		std::cout << std::endl;
        }
	}


	void physical_write_file(
			const char *i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const
	{
		request_data_physical();

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
//			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
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
        			file << physical_space_data[i*sphereDataConfig->physical_num_lat+j];
        			if (i < sphereDataConfig->physical_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}

	void file_physical_writeFile_lon_pi_shifted(
			const char *i_filename,
			std::string i_title = "",
			int i_precision = 20
	)
	{
		request_data_physical();

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
//			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
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

        			file << physical_space_data[ia*sphereDataConfig->physical_num_lat+j];
        			if (i < sphereDataConfig->physical_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
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
SphereData operator*(
		const double i_value,
		const SphereData &i_array_data
)
{
	return ((SphereData&)i_array_data)*i_value;
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
SphereData operator-(
		const double i_value,
		const SphereData &i_array_data
)
{
	return ((SphereData&)i_array_data).valueMinusThis(i_value);
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
SphereData operator+(
		const double i_value,
		const SphereData &i_array_data
)
{
	i_array_data.checkConsistency();
	return ((SphereData&)i_array_data)+i_value;
}
#endif

#endif /* SPHDATA_HPP_ */
