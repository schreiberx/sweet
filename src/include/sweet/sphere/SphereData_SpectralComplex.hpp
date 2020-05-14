/*
 * SphereDataSpectralComplex.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SPHERE_DATA_SPECTRAL_COMPLEX_HPP_
#define SPHERE_DATA_SPECTRAL_COMPLEX_HPP_

#include <complex>
#include <functional>
#include <array>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <utility>
#include <functional>

#include <sweet/MemBlockAlloc.hpp>
#include <sweet/parmemcpy.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_PhysicalComplex.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>



class SphereData_SpectralComplex
{
public:
	const SphereData_Config *sphereDataConfig;

public:
	std::complex<double> *physical_space_data;
	std::complex<double> *spectral_space_data;

	bool physical_space_data_valid;
	bool spectral_space_data_valid;

public:
	SphereData_SpectralComplex()	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
	}


public:
	SphereData_SpectralComplex(
			const SphereData_Config *i_sphereDataConfig
	)	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		assert(i_sphereDataConfig != 0);

		setup(i_sphereDataConfig);
	}


public:
	SphereData_SpectralComplex(
			const SphereData_SpectralComplex &i_sph_data
	)	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);

		operator=(i_sph_data);
	}




public:
	SphereData_SpectralComplex(
			SphereData_SpectralComplex &&i_sph_data
	)	:
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		if (i_sph_data.physical_space_data_valid)
			std::swap(physical_space_data, i_sph_data.physical_space_data);

		if (i_sph_data.spectral_space_data_valid)
			std::swap(spectral_space_data, i_sph_data.spectral_space_data);

		physical_space_data_valid = i_sph_data.physical_space_data_valid;
		spectral_space_data_valid = i_sph_data.spectral_space_data_valid;
	}





	SphereData_SpectralComplex(
			const SphereData_PhysicalComplex &i_sph_data
	):
		sphereDataConfig(nullptr),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_sph_data.physical_space_data[i];

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}




	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check_sphereDataConfig_identical_res(const SphereData_Config *i_sphereDataConfig)	const
	{
		assert(sphereDataConfig->physical_num_lat == i_sphereDataConfig->physical_num_lat);
		assert(sphereDataConfig->physical_num_lon == i_sphereDataConfig->physical_num_lon);

		assert(sphereDataConfig->spectral_modes_m_max == i_sphereDataConfig->spectral_modes_m_max);
		assert(sphereDataConfig->spectral_modes_n_max == i_sphereDataConfig->spectral_modes_n_max);
	}





public:
	SphereData_SpectralComplex& operator=(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		if (i_sph_data.physical_space_data_valid)
			parmemcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(cplx)*sphereDataConfig->physical_array_data_number_of_elements);

		if (i_sph_data.spectral_space_data_valid)
			parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(cplx)*sphereDataConfig->spectral_complex_array_data_number_of_elements);

		physical_space_data_valid = i_sph_data.physical_space_data_valid;
		spectral_space_data_valid = i_sph_data.spectral_space_data_valid;

		return *this;
	}




public:
	SphereData_SpectralComplex& operator=(
			SphereData_SpectralComplex &&i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		if (i_sph_data.physical_space_data_valid)
			std::swap(physical_space_data, i_sph_data.physical_space_data);

		if (i_sph_data.spectral_space_data_valid)
			std::swap(spectral_space_data, i_sph_data.spectral_space_data);

		physical_space_data_valid = i_sph_data.physical_space_data_valid;
		spectral_space_data_valid = i_sph_data.spectral_space_data_valid;

		return *this;
	}



	SphereData_SpectralComplex spectral_returnWithTruncatedModes(
			const SphereData_Config *i_sphereDataConfigTargetTruncation
	)	const
	{
		return spectral_returnWithDifferentModes(i_sphereDataConfigTargetTruncation).spectral_returnWithDifferentModes(sphereDataConfig);
	}


public:
	SphereData_SpectralComplex spectral_returnWithDifferentModes(
			const SphereData_Config *i_sphereDataConfigNew
	)	const
	{
		SphereData_SpectralComplex out(i_sphereDataConfigNew);

		/*
		 *  0 = invalid
		 * -1 = scale down
		 *  1 = scale up
		 */
		int scaling_mode = 0;

		if (sphereDataConfig->spectral_modes_m_max < out.sphereDataConfig->spectral_modes_m_max)
		{
			scaling_mode = 1;
		}
		else if (sphereDataConfig->spectral_modes_m_max > out.sphereDataConfig->spectral_modes_m_max)
		{
			scaling_mode = -1;
		}


		if (sphereDataConfig->spectral_modes_n_max < out.sphereDataConfig->spectral_modes_n_max)
		{
			assert(scaling_mode != -1);
			scaling_mode = 1;
		}
		else if (sphereDataConfig->spectral_modes_n_max > out.sphereDataConfig->spectral_modes_n_max)
		{
			assert(scaling_mode != 1);
			scaling_mode = -1;
		}

		if (scaling_mode == 0)
		{
			// Just copy the data
			out = *this;
			return out;
		}

		request_data_spectral();

		if (scaling_mode == -1)
		{
			/*
			 * more modes -> less modes
			 */


	SWEET_THREADING_SPACE_PARALLEL_FOR
			for (int n = 0; n <= out.sphereDataConfig->spectral_modes_n_max; n++)
			{
				int src_idx = sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
				int dst_idx = out.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);

				for (int m = -n; m <= n; m++)
				{
					out.spectral_space_data[dst_idx] = spectral_space_data[src_idx];
					src_idx++;
					dst_idx++;
				}
			}
		}
		else
		{
			/*
			 * less modes -> more modes
			 */

			// zero all values
			out.spectral_set_zero();


	SWEET_THREADING_SPACE_PARALLEL_FOR
			for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				int src_idx = sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
				int dst_idx = out.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);

				for (int m = -n; m <= n; m++)
				{
					out.spectral_space_data[dst_idx] = spectral_space_data[src_idx];
					src_idx++;
					dst_idx++;
				}
			}
		}

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}


	SphereData_PhysicalComplex getSphereDataPhysical()	const
	{
		return getSphereDataPhysicalComplex();
	}

	SphereData_PhysicalComplex getSphereDataPhysicalComplex()	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		if (physical_space_data_valid)
		{

	SWEET_THREADING_SPACE_PARALLEL_FOR
			for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
				out.physical_space_data[i] = physical_space_data[i];

			return out;
		}

		/*
		 * WARNING:
		 * We have to use a temporary array here because of destructive SH transformations
		 */
		SphereData_SpectralComplex tmp = *this;
		tmp.request_data_spectral();
		SH_to_spat_cplx(sphereDataConfig->shtns, tmp.spectral_space_data, tmp.physical_space_data);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = tmp.physical_space_data[i];

		return out;
	}



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

		SphereData_SpectralComplex *this_var = (SphereData_SpectralComplex*)this;

		this_var->physical_space_data_valid = false;
		this_var->spectral_space_data_valid = true;
	}


	const SphereData_SpectralComplex& request_data_physical()	const
	{
		if (physical_space_data_valid)
			return *this;

		assert(spectral_space_data_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		SH_to_spat_cplx(sphereDataConfig->shtns, spectral_space_data, physical_space_data);

		SphereData_SpectralComplex *this_var = (SphereData_SpectralComplex*)this;

		this_var->physical_space_data_valid = true;
		this_var->spectral_space_data_valid = false;

		return *this;
	}


	SphereData_SpectralComplex operator+(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	SphereData_SpectralComplex& operator+=(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] += i_sph_data.spectral_space_data[idx];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}


	SphereData_SpectralComplex& operator-=(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}



	SphereData_SpectralComplex operator-(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereData_SpectralComplex operator-()	const
	{
		request_data_spectral();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = -spectral_space_data[idx];

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	SphereData_SpectralComplex operator*(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out_sph_data.physical_space_data[i] = i_sph_data.physical_space_data[i]*physical_space_data[i];

		out_sph_data.spectral_space_data_valid = false;
		out_sph_data.physical_space_data_valid = true;

		return out_sph_data;
	}


	SphereData_SpectralComplex operator/(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out_sph_data.physical_space_data[i] = physical_space_data[i]/i_sph_data.physical_space_data[i];

		out_sph_data.spectral_space_data_valid = false;
		out_sph_data.physical_space_data_valid = true;

		return out_sph_data;
	}


	SphereData_SpectralComplex operator*(
			const double i_value
	)	const
	{
		request_data_spectral();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR

		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}

#if 0
	const SphereData_SpectralComplex& operator*=(
			const double i_value
	)	const
	{
		request_data_spectral();


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}
#endif


	const SphereData_SpectralComplex& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}


	SphereData_SpectralComplex operator/(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	const SphereData_SpectralComplex& operator/=(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] /= i_value;

		return *this;
	}



	SphereData_SpectralComplex operator+(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData_SpectralComplex out_sph_data(*this);

		out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	SphereData_SpectralComplex& operator+=(
			double i_value
	)
	{
		request_data_spectral();

		spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		return *this;
	}



	SphereData_SpectralComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		SphereData_SpectralComplex out_sph_data(*this);

		out_sph_data.request_data_spectral();
		out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	SphereData_SpectralComplex operator-(
			const std::complex<double> &i_value
	)	const
	{
		std::cout << "MINUS OPERATOR" << std::endl;

		SphereData_SpectralComplex out_sph_data(*this);
		out_sph_data.request_data_spectral();

		out_sph_data.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	SphereData_SpectralComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



public:
	void setup(
			const SphereData_Config *i_sphereConfig
	)
	{
		// assure that the initialization is not done twice!
		assert(sphereDataConfig == nullptr);

		sphereDataConfig = i_sphereConfig;

		physical_space_data_valid = false;
		spectral_space_data_valid = false;

		physical_space_data = MemBlockAlloc::alloc<cplx>(sphereDataConfig->physical_array_data_number_of_elements * sizeof(cplx));
		spectral_space_data = MemBlockAlloc::alloc<cplx>(sphereDataConfig->spectral_complex_array_data_number_of_elements * sizeof(cplx));
	}

public:
	~SphereData_SpectralComplex()
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
	SphereData_SpectralComplex spectral_solve_helmholtz(
			const std::complex<double> &i_a,
			const std::complex<double> &i_b,
			double r
	)	const
	{
		SphereData_SpectralComplex out(*this);

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
	SphereData_SpectralComplex& physical_truncate()
	{
		request_data_physical();

		spat_cplx_to_SH(sphereDataConfig->shtns, physical_space_data, spectral_space_data);
		SH_to_spat_cplx(sphereDataConfig->shtns, spectral_space_data, physical_space_data);

		physical_space_data_valid = true;
		spectral_space_data_valid = false;

		return *this;
	}



	/**
	 * Truncate modes which are not representable in spectral space
	 */
	SphereData_SpectralComplex& spectral_truncate()
	{
		request_data_spectral();

		SH_to_spat_cplx(sphereDataConfig->shtns, spectral_space_data, physical_space_data);
		spat_cplx_to_SH(sphereDataConfig->shtns, physical_space_data, spectral_space_data);

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}


	inline
	void spectral_update_lambda(
			std::function<void(int,int,cplx&)> i_lambda
	)
	{
		if (physical_space_data_valid)
			request_data_spectral();

		assert(spectral_space_data_valid);


		SWEET_THREADING_SPACE_PARALLEL_FOR
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
				 * WARNING: The latitude degrees are not aequidistantly spaced in the angles!!!! We have to use the shtns->ct lookup table
				 */
				//double lat_degree = M_PI*0.5 - colatitude;
				double lat_degree = sphereDataConfig->lat[j];

				i_lambda(lon_degree, lat_degree, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}
#else


		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			//double colatitude = acos(shtns->ct[j]);

			/*
			 * Colatitude is 0 at the north pole and 180 at the south pole
			 *
			 * WARNING: The latitude degrees are not aequidistantly spaced in the angles!!!! We have to use the shtns->ct lookup table
			 */
			//double lat_degree = M_PI*0.5 - colatitude;
			double lat_degree = sphereDataConfig->lat[jlat];

			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double lon_degree = ((double)ilon/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

				i_lambda(lon_degree, lat_degree, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}
#endif

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


	SWEET_THREADING_SPACE_PARALLEL_FOR

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				double mu = sphereDataConfig->lat_gaussian[j];

				i_lambda(lon_degree, mu, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}
#else
		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			double mu = sphereDataConfig->lat_gaussian[jlat];

			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double lon_degree = ((double)ilon/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

				i_lambda(lon_degree, mu, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}
#endif

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


	SWEET_THREADING_SPACE_PARALLEL_FOR

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = (((double)i)/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				double comu = sphereDataConfig->shtns->st[j];

				i_lambda(lon_degree, comu, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}
#else
		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			double comu = sphereDataConfig->shtns->st[jlat];

			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double lon_degree = (((double)ilon)/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

				i_lambda(lon_degree, comu, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}
#endif

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}




	/*
	 * Set all values to a specific value
	 */
	void physical_set_all_value(
			std::complex<double> i_value
	)
	{

	SWEET_THREADING_SPACE_PARALLEL_FOR

		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
				physical_space_data[i*sphereDataConfig->physical_num_lat + j] = i_value;

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}


	/*
	 * Set all values to zero
	 */
	void physical_set_zero()
	{
		DebugAssertSWEETError(sphereDataConfig != nullptr, "sphereDataConfig not initialized!");


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
				physical_space_data[i*sphereDataConfig->physical_num_lat + j] = 0;

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}



	SphereData_SpectralComplex physical_diff_realconst(
			const SphereData_SpectralComplex &i_sphereData
	)	const
	{
		// make a copy to avoid modifications
		SphereData_SpectralComplex a = SphereData_SpectralComplex(i_sphereData);
		SphereData_SpectralComplex b = SphereData_SpectralComplex(*this);

		a.request_data_physical();
		b.request_data_physical();

		SphereData_SpectralComplex out(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			out.physical_space_data[j] = a.physical_space_data[j] - b.physical_space_data[j];

		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;

		return out;
	}



	/**
	 * Return the maximum error norm
	 */
	double physical_reduce_max(
			const SphereData_SpectralComplex &i_data
	)
	{
		check_sphereDataConfig_identical_res(i_data.sphereDataConfig);

		request_data_physical();
		i_data.request_data_physical();

		double error = -1;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
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
	double physical_reduce_max_abs()
	{
		request_data_physical();

		double error = -1;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
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


	double physical_reduce_rms()
	{
		request_data_physical();

		double error = 0;

		for (std::size_t j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			std::complex<double> &d = physical_space_data[j];

			// use amplitude
			error += (std::conj(d)*d).real();
		}

		return std::sqrt(error / (double)sphereDataConfig->physical_array_data_number_of_elements);
	}




	/*
	 * Set all values to zero in spectral space
	 */
	void spectral_set_zero()
	{

	SWEET_THREADING_SPACE_PARALLEL_FOR

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



public:

	void physical_file_write(
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
        			file << physical_space_data[i*sphereDataConfig->physical_num_lat+j];
        			if (i < sphereDataConfig->physical_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}



	void physical_file_write_real(
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
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;
			lon_degree = lon_degree/M_PI*180.0;

			file << lon_degree;
			if (i < sphereDataConfig->physical_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
        {
        		double lat_degree = sphereDataConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		file << lat_degree << "\t";

        		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
        		{
        			file << physical_space_data[i*sphereDataConfig->physical_num_lat+j].real();
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
			double lon_degree = ((double)i/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;
			lon_degree = (lon_degree-M_PI)/M_PI*180.0;

			file << lon_degree;
			if (i < sphereDataConfig->physical_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
        {
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
SphereData_SpectralComplex operator*(
		const double i_value,
		const SphereData_SpectralComplex &i_array_data
)
{
	return ((SphereData_SpectralComplex&)i_array_data)*i_value;
}


/**
 * operator to support operations such as:
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SphereData_SpectralComplex operator*(
		const std::complex<double> &i_value,
		const SphereData_SpectralComplex &i_array_data
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
SphereData_SpectralComplex operator-(
		const double i_value,
		const SphereData_SpectralComplex &i_array_data
)
{
	return ((SphereData_SpectralComplex&)i_array_data).valueMinusThis(i_value);
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
SphereData_SpectralComplex operator+(
		const double i_value,
		const SphereData_SpectralComplex &i_array_data
)
{
	return ((SphereData_SpectralComplex&)i_array_data)+i_value;
}

inline
static
SphereData_SpectralComplex operator+(
		const std::complex<double> &i_value,
		const SphereData_SpectralComplex &i_array_data
)
{
	return i_array_data+i_value;
}

inline
static
SphereData_SpectralComplex operator-(
		const std::complex<double> &i_value,
		const SphereData_SpectralComplex &i_array_data
)
{
	i_array_data.request_data_spectral();

	SphereData_SpectralComplex out_sph_data(i_array_data.sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < i_array_data.sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = -i_array_data.spectral_space_data[idx];

	out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	out_sph_data.physical_space_data_valid = false;
	out_sph_data.spectral_space_data_valid = true;

	return out_sph_data;

}

#endif /* SPHDATA_HPP_ */
