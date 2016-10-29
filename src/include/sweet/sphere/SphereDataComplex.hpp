/*
 * SPHDataComplex.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHDATACOMPLEX_HPP_
#define SPHDATACOMPLEX_HPP_

#include <complex>
#include <functional>
#include <array>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/MemBlockAlloc.hpp>



class SphereDataComplex
{
public:
	SphereDataConfig *sphereDataConfig;

public:
	std::complex<double> *physical_space_data;
	std::complex<double> *spectral_space_data;

	bool physical_space_data_valid;
	bool spectral_space_data_valid;

public:
	SphereDataComplex()	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
	}


public:
	SphereDataComplex(
			SphereDataConfig *i_sphConfig
	)	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		assert(i_sphConfig != 0);

		setup(i_sphConfig);
	}


public:
	SphereDataComplex(
			const SphereDataComplex &i_sph_data
	)	:
		sphereDataConfig(nullptr),
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
	inline void check_sphereDataConfig_identical_res(SphereDataConfig *i_sphConfig)	const
	{
		assert(sphereDataConfig->physical_num_lat == i_sphConfig->physical_num_lat);
		assert(sphereDataConfig->physical_num_lon == i_sphConfig->physical_num_lon);

		assert(sphereDataConfig->spectral_modes_m_max == i_sphConfig->spectral_modes_m_max);
		assert(sphereDataConfig->spectral_modes_n_max == i_sphConfig->spectral_modes_n_max);
	}



#if 0
public:
	SphereDataComplex& operator=(
			const SphereData &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		assert(sphereDataConfig == i_sph_data.sphereDataConfig);

#warning	"TODO: maybe replace this with assignment in spectral space!"
		physical_fromSphereData(i_sph_data);
		return *this;
	}
#endif



public:
	SphereDataComplex& operator=(
			const SphereDataComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		if (i_sph_data.physical_space_data_valid)
			memcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(cplx)*sphereDataConfig->physical_array_data_number_of_elements);

		if (i_sph_data.spectral_space_data_valid)
			memcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(cplx)*sphereDataConfig->spectral_complex_array_data_number_of_elements);

		physical_space_data_valid = i_sph_data.physical_space_data_valid;
		spectral_space_data_valid = i_sph_data.spectral_space_data_valid;

		return *this;
	}


#if 0
public:
	void physical_RealToSphereData(
			SphereData &o_sph_data
	)
	{
		check_sphereDataConfig_identical_res(o_sph_data.sphereDataConfig);

		request_data_physical();

		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			o_sph_data.physical_space_data[i] = physical_space_data[i].real();

		o_sph_data.physical_space_data_valid = true;
		o_sph_data.spectral_space_data_valid = false;
	}

public:
	void physical_ImagToSphereData(
			SphereData &o_sph_data
	)
	{
		check_sphereDataConfig_identical_res(o_sph_data.sphereDataConfig);

		request_data_physical();

		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			o_sph_data.physical_space_data[i] = physical_space_data[i].imag();

		o_sph_data.physical_space_data_valid = true;
		o_sph_data.spectral_space_data_valid = false;
	}


public:
	void physical_fromSphereData(
			const SphereData &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		i_sph_data.request_data_physical();

		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_sph_data.physical_space_data[i];

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}
#endif


public:
	void request_data_spectral()	const
	{
		if (spectral_space_data_valid)
			return;

		assert(physical_space_data_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		spat_cplx_to_SH(sphereDataConfig->shtns, physical_space_data, spectral_space_data);

		SphereDataComplex *this_var = (SphereDataComplex*)this;

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
		SH_to_spat_cplx(sphereDataConfig->shtns, spectral_space_data, physical_space_data);

		SphereDataComplex *this_var = (SphereDataComplex*)this;

		this_var->physical_space_data_valid = true;
		this_var->spectral_space_data_valid = false;
	}


	SphereDataComplex operator+(
			const SphereDataComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereDataComplex out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereDataComplex& operator+=(
			const SphereDataComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] += i_sph_data.spectral_space_data[idx];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}


	SphereDataComplex& operator-=(
			const SphereDataComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}



	SphereDataComplex operator-(
			const SphereDataComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereDataComplex out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereDataComplex operator-()
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = -spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereDataComplex operator*(
			const SphereDataComplex &i_sph_data
	)	const
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SphereDataComplex out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int i = 0; i < sphereDataConfig->spectral_complex_array_data_number_of_elements; i++)
			out_sph_data.physical_space_data[i] = i_sph_data.physical_space_data[i]*physical_space_data[i];

		out_sph_data.spectral_space_data_valid = false;
		out_sph_data.physical_space_data_valid = true;

		return out_sph_data;
	}


	SphereDataComplex operator*(
			const double i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	const SphereDataComplex& operator*=(
			const double i_value
	)	const
	{
		request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}


	const SphereDataComplex& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}


	SphereDataComplex operator/(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereDataComplex operator+(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(*this);

		out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	SphereDataComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		SphereDataComplex out_sph_data(*this);

		out_sph_data.request_data_spectral();
		out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	SphereDataComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphereDataConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



public:
	void setup(SphereDataConfig *i_sphConfig)
	{
		// assure that the initialization is not done twice!
		assert(sphereDataConfig == nullptr);

		sphereDataConfig = i_sphConfig;

		physical_space_data_valid = false;
		spectral_space_data_valid = false;

		physical_space_data = MemBlockAlloc::alloc<cplx>(sphereDataConfig->physical_array_data_number_of_elements * sizeof(cplx));
		spectral_space_data = MemBlockAlloc::alloc<cplx>(sphereDataConfig->spectral_complex_array_data_number_of_elements * sizeof(cplx));
	}

public:
	~SphereDataComplex()
	{
		MemBlockAlloc::free(physical_space_data, sphereDataConfig->physical_array_data_number_of_elements * sizeof(cplx));
		MemBlockAlloc::free(spectral_space_data, sphereDataConfig->spectral_complex_array_data_number_of_elements * sizeof(cplx));
	}



	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	inline
	SphereDataComplex spectral_solve_helmholtz(
			const std::complex<double> &i_a,
			const std::complex<double> &i_b,
			double r
	)	const
	{
		SphereDataComplex out(*this);

		out.request_data_spectral();

		const std::complex<double> a = i_a;
		const std::complex<double> b = i_b/(r*r);

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
	SphereDataComplex physical_truncate()
	{
		request_data_physical();

		SphereDataComplex out_sph_data(sphereDataConfig);
		spat_cplx_to_SH(sphereDataConfig->shtns, physical_space_data, out_sph_data.spectral_space_data);
		SH_to_spat_cplx(sphereDataConfig->shtns, out_sph_data.spectral_space_data, out_sph_data.physical_space_data);

		out_sph_data.physical_space_data_valid = true;
		out_sph_data.spectral_space_data_valid = false;

		return out_sph_data;
	}



	/**
	 * Truncate modes which are not representable in spectral space
	 */
	SphereDataComplex spectral_truncate()
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphereDataConfig);
		SH_to_spat_cplx(sphereDataConfig->shtns, spectral_space_data, out_sph_data.physical_space_data);
		spat_cplx_to_SH(sphereDataConfig->shtns, out_sph_data.physical_space_data, out_sph_data.spectral_space_data);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	inline
	void spectral_update_lambda(
			std::function<void(int,int,cplx&)> i_lambda
	)
	{
		if (physical_space_data_valid)
			request_data_spectral();

		assert(spectral_space_data_valid);

#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				i_lambda(n, m, spectral_space_data[idx]);
				idx++;
			}
		}
	}



	void physical_update_lambda_sinphi_grid(
			std::function< void(double,double,std::complex<double>&) > i_lambda	///< lambda function to return value for lat/mu
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



	inline
	const std::complex<double>& spectral_get(
			int in,
			int im
	)	const
	{
		static const std::complex<double> zero = {0,0};

		assert(spectral_space_data_valid);

		if (in < 0)
			return zero;

		if (in > sphereDataConfig->spectral_modes_n_max)
			return zero;

		if (std::abs(im) > sphereDataConfig->spectral_modes_m_max)
			return zero;

		if (std::abs(im) > in)
			return zero;

		assert (im <= sphereDataConfig->spectral_modes_m_max);
		return spectral_space_data[sphereDataConfig->getArrayIndexByModes_Complex(in, im)];
	}



	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */
public:
	void physical_update_lambda(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
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
				 * WARNING: The latitude degrees are not aequidistantly spaced in the angles!!!! We have to use the shtns->ct lookup table
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
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-1;1])
	 */
	void physical_update_lambda_gaussian_grid(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
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
				double mu = sphereDataConfig->lat_gaussian[j];

				i_lambda(lon_degree, mu, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}




	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-1;1])
	 */

	void physical_update_lambda_cogaussian_grid(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
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
				double comu = sphereDataConfig->shtns->st[j];
				/*
				 * IDENTITAL FORMULATION
				double mu = shtns->ct[j];
				double comu = sqrt(1.0-mu*mu);
				*/

				i_lambda(lon_degree, comu, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
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


	/**
	 * Return the maximum error norm
	 */
	double physical_reduce_error_max(
			const SphereDataComplex &i_data
	)
	{
		check_sphereDataConfig_identical_res(i_data.sphereDataConfig);

		request_data_physical();
		i_data.request_data_physical();

		double error = -1;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						error,
						std::abs(
								physical_space_data[j].real() - i_data.physical_space_data[j].real()
							)
						);

			error = std::max(
						error,
						std::abs(
								physical_space_data[j].imag() - i_data.physical_space_data[j].imag()
							)
						);
		}
		
		return error;
	}


	/**
	 * Return the maximum error norm
	 */
	double physical_reduce_error_max()
	{
		request_data_physical();

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




	/*
	 * Set all values to zero in spectral space
	 */
	void spectral_set_zero()
	{
#if SWEET_THREADING
#pragma omp parallel for
#endif

		for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				spectral_space_data[idx] = 0;
				idx++;
			}
		}

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
	}


	void spectral_print(
			int i_precision = 8
	)	const
	{
		request_data_spectral();

		std::cout << std::setprecision(i_precision);

		/**
		 * WARNING: This follows a different order contrast to how it is stored
		 */
		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
				std::cout << spectral_space_data[idx] << "\t";
			}
			std::cout << std::endl;
		}
	}


	void physical_print(
			int i_precision = 8
	)	const
	{
		request_data_physical();

		std::cout << std::setprecision(i_precision);

#if 0
		for (std::size_t i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;
			lon_degree = lon_degree/M_PI*180.0;

			std::cout << lon_degree;
			if (i < sphereDataConfig->physical_num_lon-1)
				std::cout << "\t";
		}
		std::cout << std::endl;
#endif

        for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
        {
#if 0
        		double lat_degree = sphereDataConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		std::cout << lat_degree << "\t";
#endif
        		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
        		{
        			std::cout << physical_space_data[i*sphereDataConfig->physical_num_lat+j];
        			if (i < sphereDataConfig->physical_num_lon-1)
        				std::cout << "\t";
        		}
        		std::cout << std::endl;
        }
	}



public:

	void physical_write_file(
			const char *i_filename,
			const char *i_title = "",
			int i_precision = 8
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

	void physical_write_file_lon_pi_shifted(
			const char *i_filename,
			std::string i_title = "",
			int i_precision = 8
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
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SphereDataComplex operator*(
		const double i_value,
		const SphereDataComplex &i_array_data
)
{
	return ((SphereDataComplex&)i_array_data)*i_value;
}


/**
 * operator to support operations such as:
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SphereDataComplex operator*(
		const std::complex<double> &i_value,
		const SphereDataComplex &i_array_data
)
{
	return i_array_data*i_value;
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
SphereDataComplex operator-(
		const double i_value,
		const SphereDataComplex &i_array_data
)
{
	return ((SphereDataComplex&)i_array_data).valueMinusThis(i_value);
//	return -(((SPHDataComplex&)i_array_data).operator-(i_value));
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

inline
static
SphereDataComplex operator+(
		const double i_value,
		const SphereDataComplex &i_array_data
)
{
	return ((SphereDataComplex&)i_array_data)+i_value;
}

inline
static
SphereDataComplex operator+(
		const std::complex<double> &i_value,
		const SphereDataComplex &i_array_data
)
{
	return i_array_data+i_value;
}

#endif /* SPHDATA_HPP_ */
