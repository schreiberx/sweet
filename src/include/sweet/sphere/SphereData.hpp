/*
 * SphereData.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SPHERE_DATA_HPP_
#define SPHERE_DATA_HPP_

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

#include <sweet/sweetmath.hpp>
#include <sweet/parmemcpy.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereDataPhysical.hpp>
#include <sweet/sphere/SphereDataPhysicalComplex.hpp>
#include <sweet/FatalError.hpp>
#include <sweet/openmp_helper.hpp>



class SphereData
{
	friend class SphereDataComplex;

public:
	const SphereDataConfig *sphereDataConfig = nullptr;

public:
	double *physical_space_data = nullptr;
	std::complex<double> *spectral_space_data = nullptr;

	bool physical_space_data_valid;
	bool spectral_space_data_valid;


	void swap(
			SphereData &i_sphereData
	)
	{
		assert(sphereDataConfig == i_sphereData.sphereDataConfig);

		std::swap(physical_space_data, i_sphereData.physical_space_data);
		std::swap(spectral_space_data, i_sphereData.spectral_space_data);
		std::swap(spectral_space_data_valid, i_sphereData.spectral_space_data_valid);
		std::swap(physical_space_data_valid, i_sphereData.physical_space_data_valid);
	}

public:
	SphereData(
			const SphereDataConfig *i_sphereDataConfig
	)	:
		sphereDataConfig(i_sphereDataConfig),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		assert(i_sphereDataConfig != 0);

		setup(i_sphereDataConfig);
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
		// Dummy initialization
		if (i_sph_data.sphereDataConfig == nullptr)
			return;

		setup(i_sph_data.sphereDataConfig);

		operator=(i_sph_data);
	}



public:
	SphereData(
			SphereData &&i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		physical_space_data(nullptr),
		spectral_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);

		if (i_sph_data.physical_space_data_valid)
			std::swap(physical_space_data, i_sph_data.physical_space_data);

		if (i_sph_data.spectral_space_data_valid)
			std::swap(spectral_space_data, i_sph_data.spectral_space_data);

		physical_space_data_valid = i_sph_data.physical_space_data_valid;
		spectral_space_data_valid = i_sph_data.spectral_space_data_valid;
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

		assert(sphereDataConfig->spectral_modes_m_max == i_sphereDataConfig->spectral_modes_m_max);
		assert(sphereDataConfig->spectral_modes_n_max == i_sphereDataConfig->spectral_modes_n_max);
	}


public:
	SphereData& operator=(
			const SphereData &i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		if (i_sph_data.physical_space_data_valid)
			parmemcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		if (i_sph_data.spectral_space_data_valid)
			parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(cplx)*sphereDataConfig->spectral_array_data_number_of_elements);

		physical_space_data_valid = i_sph_data.physical_space_data_valid;
		spectral_space_data_valid = i_sph_data.spectral_space_data_valid;

		return *this;
	}




	/**
	 * This function implements copying the spectral data only.
	 *
	 * This becomes handy if coping with data which should be only transformed without dealiasing.
	 */
public:
	SphereData& load_nodealiasing(
			const SphereData &i_sph_data		///< data to be converted to sphereDataConfig_nodealiasing
	)
	{
		if (sphereDataConfig == nullptr)
			FatalError("sphereDataConfig not initialized");

		// only copy spectral data
		i_sph_data.request_data_spectral();

		parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(cplx)*sphereDataConfig->spectral_array_data_number_of_elements);

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}


public:
	SphereData& operator=(
			SphereData &&i_sph_data
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


public:
	SphereData spectral_returnWithDifferentModes(
			const SphereDataConfig *i_sphereDataConfig
	)	const
	{
		SphereData out(i_sphereDataConfig);

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
			for (int m = 0; m <= out.sphereDataConfig->spectral_modes_m_max; m++)
			{
				cplx *dst = &out.spectral_space_data[out.sphereDataConfig->getArrayIndexByModes(m, m)];
				cplx *src = &spectral_space_data[sphereDataConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(cplx)*(out.sphereDataConfig->spectral_modes_n_max-m+1);
				memcpy(dst, src, size);
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
			for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
			{
				cplx *dst = &out.spectral_space_data[out.sphereDataConfig->getArrayIndexByModes(m, m)];
				cplx *src = &spectral_space_data[sphereDataConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(cplx)*(sphereDataConfig->spectral_modes_n_max-m+1);
				memcpy(dst, src, size);
			}
		}

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
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


	SphereDataPhysical getSphereDataPhysical()	const
	{
		SphereDataPhysical out(sphereDataConfig);

		if (physical_space_data_valid)
		{
			parmemcpy(out.physical_space_data, physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);
			return out;
		}

		/*
		 * WARNING:
		 * We have to use a temporary array here because of destructive SH transformations
		 */
		SphereData tmp = *this;
		tmp.request_data_spectral();
		SH_to_spat(sphereDataConfig->shtns, tmp.spectral_space_data, out.physical_space_data);

		return out;
	}


	SphereDataPhysicalComplex getSphereDataPhysicalComplex()	const
	{
		SphereDataPhysicalComplex out(sphereDataConfig);

		if (physical_space_data_valid)
		{
			parmemcpy(out.physical_space_data, physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);
			return out;
		}

		/*
		 * WARNING:
		 * We have to use a temporary array here because of destructive SH transformations
		 */
		SphereData tmp = *this;
		tmp.request_data_spectral();
		SH_to_spat(sphereDataConfig->shtns, tmp.spectral_space_data, tmp.physical_space_data);

		parmemcpy(out.physical_space_data, tmp.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		return out;
	}



	SphereData(
			const SphereDataPhysical &i_sph_data
	)
	{
		setup(i_sph_data.sphereDataConfig);

		parmemcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}



	SphereData(
			const SphereDataPhysicalComplex &i_sph_data
	)
	{
		setup(i_sph_data.sphereDataConfig);

		parmemcpy(physical_space_data, i_sph_data.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}



	SphereData robert_convertToRobert()
	{
		SphereData out = *this;

		out.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					io_data *= std::cos(i_lat);
				}
		);

		return out;
	}



	SphereData robert_convertToNonRobert()
	{
		SphereData out = *this;

		out.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					io_data /= std::cos(i_lat);
				}
		);

		return out;
	}

	SphereData robert_convertToNonRobertSquared()
	{
		SphereData out = *this;

		out.physical_update_lambda_cosphi_grid(
				[](double i_lon, double i_cosphi, double &io_data)
				{
					io_data /= i_cosphi*i_cosphi;
				}
		);

		return out;
	}


	void physical_add(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] += i_sph_data.physical_space_data[idx];
	}



	void physical_sub(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->physical_array_data_number_of_elements; idx++)
			physical_space_data[idx] -= i_sph_data.physical_space_data[idx];
	}


	SphereData operator+(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData out(sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}



	SphereData& operator+=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
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

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

		physical_space_data_valid = false;
		spectral_space_data_valid = true;

		return *this;
	}



	SphereData operator-(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}



	SphereData operator-()
	{
		request_data_spectral();

		SphereData out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = -spectral_space_data[idx];

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}



	SphereData operator*(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData a = *this;
		SphereData b = i_sph_data;

		a.request_data_physical();
		b.request_data_physical();

		SphereData out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = a.physical_space_data[i]*b.physical_space_data[i];

		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;

		// directly convert back to spectral space to truncate modes
		//out.request_data_spectral();

		return out;
	}



	SphereData operator/(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData a = *this;
		SphereData b = i_sph_data;

		a.request_data_physical();
		b.request_data_physical();

		SphereData out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = a.physical_space_data[i]/b.physical_space_data[i];

		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;

		// directly convert back to spectral space to truncate modes
		out.request_data_spectral();

		return out;
	}



	SphereData operator*(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}




	const SphereData& operator*=(
			double i_value
	)	const
	{
		request_data_spectral();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}





	const SphereData& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}



	SphereData operator/(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}


	SphereData operator+(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out(*this);

		out.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}




	SphereData operator-(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out(*this);

		out.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}




public:
	void setup(
			const SphereDataConfig *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		spectral_space_data_valid = false;
		physical_space_data_valid = false;

		physical_space_data = MemBlockAlloc::alloc<double>(sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));
		spectral_space_data = MemBlockAlloc::alloc<cplx>(sphereDataConfig->spectral_array_data_number_of_elements * sizeof(cplx));
	}



public:
	~SphereData()
	{
		if (physical_space_data != nullptr)
			MemBlockAlloc::free(physical_space_data, sphereDataConfig->physical_array_data_number_of_elements * sizeof(double));

		if (spectral_space_data != nullptr)
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



	/**
	 * Solve a Laplace problem given by
	 *
	 * (D^2) x = rhs
	 */
	inline
	SphereData spectral_solve_laplace(
			double r
	)
	{
		SphereData out(*this);

		const double b = 1.0/(r*r);

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				if (n == 0)
					io_data = 0;
				else
					io_data /= (-b*(double)n*((double)n+1.0));
			}
		);

		return out;
	}



public:
	/**
	 * Truncate modes which are not representable in spectral space
	 */
	const SphereData& physical_truncate()
	{
		request_data_physical();

		spat_to_SH(sphereDataConfig->shtns, physical_space_data, spectral_space_data);
		SH_to_spat(sphereDataConfig->shtns, spectral_space_data, physical_space_data);

		return *this;
	}



	/**
	 * Truncate modes which are not representable in spectral space
	 */
	const SphereData& spectral_truncate()	const
	{
		request_data_spectral();

		SH_to_spat(sphereDataConfig->shtns, spectral_space_data, physical_space_data);
		spat_to_SH(sphereDataConfig->shtns, physical_space_data, spectral_space_data);

		return *this;
	}




	void spectral_update_lambda(
			std::function<void(int,int,cplx&)> i_lambda
	)
	{
		if (physical_space_data_valid)
			request_data_spectral();

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
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
		SphereData a = *this;
		a.request_data_physical();

		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
		{
#if __GNUC__ == 5
			if (isnan(a.physical_space_data[i]) || isinf(a.physical_space_data[i]) != 0)
#else
			if (std::isnan(a.physical_space_data[i]) || std::isinf(a.physical_space_data[i]) != 0)
#endif
				return true;
		}

		return false;
	}


	const std::complex<double>& spectral_get(
			int i_n,
			int i_m
	)	const
	{
		static const std::complex<double> zero = {0,0};

		assert(spectral_space_data_valid);

//		if (i_n < 0)
//			i_n = i_n-1;

		if (i_n < 0 ||  i_m < 0)
			return zero;

		if (i_n > sphereDataConfig->spectral_modes_n_max)
			return zero;

		if (i_m > sphereDataConfig->spectral_modes_m_max)
			return zero;

		if (i_m > i_n)
			return zero;

		assert (i_m <= sphereDataConfig->spectral_modes_m_max);
		return spectral_space_data[sphereDataConfig->getArrayIndexByModes(i_n, i_m)];
	}


	double p_physical_get(
			int i_lon,
			int i_lat
	)	const
	{
		request_data_physical();

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
		return physical_space_data[i_lon*sphereDataConfig->physical_num_lat + i_lat];
#else
		return physical_space_data[i_lat*sphereDataConfig->physical_num_lon + i_lon];
#endif
	}


	/*
	 * Set all values to zero
	 */
	void spectral_set_zero()
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] = {0,0};

		physical_space_data_valid = false;
		spectral_space_data_valid = true;
	}




	/*
	 * Set all values to value
	 */
	void spectral_set_value(
			const std::complex<double> &i_value
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] = i_value;

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


#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR
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

		SWEET_THREADING_SPACE_PARALLEL_FOR
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
		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}


	void physical_update_lambda_array(
			std::function<void(int,int,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		if (spectral_space_data_valid)
			request_data_physical();

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				i_lambda(i, j, physical_space_data[i*sphereDataConfig->physical_num_lat + j]);
			}
		}

#else

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				i_lambda(ilon, jlat, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}

#endif

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
		physical_space_data_valid = true;
		spectral_space_data_valid = false;
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

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
		{
			double lon_degree = (((double)i)/(double)sphereDataConfig->physical_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			{
				double cos_phi = sphereDataConfig->lat_cogaussian[j];
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
				i_lambda(lon_degree, cos_phi, physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon]);
			}
		}
#endif

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
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
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = 0;

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}



	/*
	 * Set all values to a specific value
	 */
	void physical_set_all_value(
			double i_value
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			physical_space_data[i] = i_value;

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
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
		if (spectral_space_data_valid)
			request_data_physical();


#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
		physical_space_data[i_lon_idx*sphereDataConfig->physical_num_lat + i_lat_idx] = i_value;
#else
		physical_space_data[i_lat_idx*sphereDataConfig->physical_num_lon + i_lon_idx] = i_value;
#endif

		physical_space_data_valid = true;
		spectral_space_data_valid = false;
	}



	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double physical_reduce_max_abs(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		double error = -1;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
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
		request_data_physical();

		double error = 0;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			double &d = physical_space_data[j];
			error += d*d;
		}

		return std::sqrt(error / (double)sphereDataConfig->physical_array_data_number_of_elements);
	}



#if 0
	SphereData physical_diff_realconst(
			const SphereData &i_sphereData
	)	const
	{
		// make a copy to avoid modifications
		SphereData a = i_sphereData;
		SphereData b = *this;

		a.request_data_physical();
		b.request_data_physical();

		SphereData out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			out.physical_space_data[j] = a.physical_space_data[j] - b.physical_space_data[j];

		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;

		return out;
	}
#endif


	/**
	 * return the sum of all values=
	 */
	double physical_reduce_sum()	const
	{
		request_data_physical();

		double sum = 0;
		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			sum += physical_space_data[j];

		return sum;
	}


	/**
	 * return the sum of all values, use quad precision for reduction
	 */
	double physical_reduce_sum_quad()	const
	{
		request_data_physical();

		double sum = 0;
		double c = 0;
#if SWEET_THREADING_SPACE
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
	 * return the sum of all values multiplied with an increasing index, use quad precision for reduction
	 *
	 * This is only helpful for debugging purpose!
	 */
	double physical_reduce_debug_sum_quad_mul_increasing()	const
	{
		request_data_physical();

		double sum = 0;
		double c = 0;
#if SWEET_THREADING_SPACE
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


	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	std::complex<double> spectral_reduce_sum_quad_increasing()	const
	{
		request_data_spectral();

		std::complex<double> sum = 0;
		std::complex<double> c = 0;
#if SWEET_THREADING_SPACE
//#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
#endif
		for (int i = 0; i < sphereDataConfig->spectral_array_data_number_of_elements; i++)
		{
			std::complex<double> value = spectral_space_data[i]*(double)i;

			// Use Kahan summation
			std::complex<double> y = value - c;
			std::complex<double> t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return sum;
	}


	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	std::complex<double> spectral_reduce_sum_quad()	const
	{
		request_data_spectral();

		std::complex<double> sum = 0;
		std::complex<double> c = 0;
#if SWEET_THREADING_SPACE
//#pragma omp parallel for proc_bind(close) reduction(+:sum,c)
#endif
		for (int i = 0; i < sphereDataConfig->spectral_array_data_number_of_elements; i++)
		{
			std::complex<double> value = spectral_space_data[i];

			// Use Kahan summation
			std::complex<double> y = value - c;
			std::complex<double> t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return sum;
	}


	/**
	 * Return the minimum value
	 */
	std::complex<double> spectral_reduce_min()	const
	{
		request_data_spectral();

		std::complex<double> error = std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->spectral_array_data_number_of_elements; j++)
		{
			error.real(std::min(spectral_space_data[j].real(), error.real()));
			error.imag(std::min(spectral_space_data[j].imag(), error.imag()));
		}

		return error;
	}


	/**
	 * Return the minimum value
	 */
	std::complex<double> spectral_reduce_max()	const
	{
		request_data_spectral();

		std::complex<double> error = -std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->spectral_array_data_number_of_elements; j++)
		{
			error.real(std::max(spectral_space_data[j].real(), error.real()));
			error.imag(std::max(spectral_space_data[j].imag(), error.imag()));
		}

		return error;
	}


	double physical_reduce_sum_metric()
	{
		request_data_physical();

		FatalError("TODO: Implement metric-scaled summation");
		double sum = 0;
		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			sum += physical_space_data[j];
		}

		return sum;
	}



	/**
	 * Return the maximum absolute value
	 */
	double physical_reduce_max_abs()	const
	{
		request_data_physical();

		double error = -1;

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
		{
			error = std::max(
						std::abs(physical_space_data[j]),
						error
					);
		}
		return error;
	}


	/**
	 * Return the minimum value
	 */
	double physical_reduce_min()	const
	{
		request_data_physical();

		double error = std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			error = std::min(physical_space_data[j], error);

		return error;
	}


	/**
	 * Return the minimum value
	 */
	double physical_reduce_max()	const
	{
		request_data_physical();

		double error = -std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			error = std::max(physical_space_data[j], error);

		return error;
	}


	/**
	 * Rescale data so that max_abs returns the given value
	 */
	SphereData physical_rescale_to_max_abs(
			double i_new_max_abs
	)	const
	{
		double max_abs = physical_reduce_max_abs();
		double scale = i_new_max_abs/max_abs;

		SphereData out(sphereDataConfig);
		request_data_physical();

		for (int j = 0; j < sphereDataConfig->physical_array_data_number_of_elements; j++)
			out.physical_space_data[j] = physical_space_data[j]*scale;

		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;

		return out;
	}



	void spectral_print(
			int i_precision = 20
	)	const
	{
		request_data_spectral();

		std::cout << std::setprecision(i_precision);

		std::cout << "m \\ n ----->" << std::endl;
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
		request_data_physical();

		std::ofstream file(i_filename, std::ios_base::trunc);

		file << std::setprecision(i_precision);
		file << "#TI " << i_title << std::endl;
		file << "#TX Longitude" << std::endl;
		file << "#TY Latitude" << std::endl;

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

  
  	void spectrum_file_write(
			const std::string &i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const
	{
  		request_data_spectral();

  		std::ofstream file(i_filename, std::ios_base::trunc);

  		file << std::setprecision(i_precision);
  		file << "#TI " << i_title << std::endl;

  		// Use 0 to make it processable by python
  		file << "0\t";

  		std::complex<double> w = {0,0};
  		std::vector<double> sum(sphereDataConfig->spectral_modes_n_max+1,0);
  		std::vector<double> sum_squared(sphereDataConfig->spectral_modes_n_max+1,0);
  		std::vector<double> max(sphereDataConfig->spectral_modes_n_max+1,0);

  		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
  		{
  			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
  			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
  			{
  				w = spectral_space_data[idx];
  				idx++;

  				sum[n]         += std::abs(w);
  				sum_squared[n] += std::abs(w * w);
  				if (std::abs(w) >= max[n]) max[n] = std::abs(w);
  			}
  		}

  		file << sphereDataConfig->spectral_modes_m_max << " "
  				<< sphereDataConfig->spectral_modes_n_max << std::endl;
  		for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
  			file << n << " " << sum[n] << " " << max[n] << " " << std::sqrt(sum_squared[n]) << std::endl;

  		file.close();
	}




  	/**
  	 * Write the spectral data to a file in binary format
  	 */
  	void file_write_binary_spectral(
			const std::string &i_filename
	)	const
	{
  		request_data_spectral();

  		std::ofstream file(i_filename, std::ios_base::trunc | std::ios_base::binary);

		if (!file.is_open())
			FatalError("Error while opening file");

  		file << "SWEET" << std::endl;
  		file << "DATA_TYPE SH_DATA" << std::endl;
  		file << "NUM_LON " << sphereDataConfig->spectral_modes_m_max << std::endl;
  		file << "NUM_LAT " << sphereDataConfig->spectral_modes_n_max << std::endl;
  		file << "SIZE " << sphereDataConfig->spectral_array_data_number_of_elements << std::endl;
  		file << "FIN" << std::endl;
  		std::cout << file.tellp() << std::endl;

  		file.write((const char*)spectral_space_data, sizeof(std::complex<double>)*sphereDataConfig->spectral_array_data_number_of_elements);

  		file.close();
	}


  	void file_read_binary_spectral(
			const std::string &i_filename
	)
	{
  		request_data_spectral();

  		std::ifstream file(i_filename, std::ios_base::binary);

		if (!file.is_open())
			FatalError("Error while opening file");

  		std::string magic;
  		std::getline(file, magic);

  		if (magic != "SWEET")
  			FatalError("Magic code 'SWEET' not found");

  		std::string data_type;
  		int num_lon = -1;
  		int num_lat = -1;
  		int size = -1;

  		while (true)
  		{
  			std::string buf;
  			file >> buf;

  			if (buf == "FIN")
  				break;

  			if (buf == "DATA_TYPE")
  			{
  				// load data type
  	  			file >> data_type;
  				std::cout << data_type << std::endl;
  	  			continue;
  			}

  			if (buf == "NUM_LON")
  			{
  				file >> buf;
  				num_lon = std::stoi(buf);
  				std::cout << num_lon << std::endl;
  				continue;
  			}

  			if (buf == "NUM_LAT")
  			{
  				file >> buf;
  				num_lat = std::stoi(buf);
  				std::cout << num_lat << std::endl;
  				continue;
  			}

  			if (buf == "SIZE")
  			{
  				file >> buf;
  				size = std::stoi(buf);
  				std::cout << size << std::endl;
  				continue;
  			}

  			FatalError("Unknown Tag '"+buf+"'");
  		}

  		// read last newline
  		char nl;
  		file.read(&nl, 1);
  		std::cout << file.tellg() << std::endl;

  		if (data_type != "SH_DATA")
  			FatalError("Unknown data type '"+data_type+"'");

  		if (num_lon != sphereDataConfig->spectral_modes_m_max)
  			FatalError("NUM_LON "+std::to_string(num_lon)+" doesn't match SphereDataConfig");

  		if (num_lat != sphereDataConfig->spectral_modes_n_max)
  			FatalError("NUM_LAT "+std::to_string(num_lat)+" doesn't match SphereDataConfig");

  		file.read((char*)spectral_space_data, sizeof(std::complex<double>)*sphereDataConfig->spectral_array_data_number_of_elements);

  		file.close();

  		spectral_space_data_valid = true;
  		physical_space_data_valid = false;
	}

	void physical_file_write_lon_pi_shifted(
			const char *i_filename,
			std::string i_title = "",
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

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
			physical_space_data_valid = true;
			spectral_space_data_valid = false;
#endif
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
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		physical_space_data_valid = true;
		spectral_space_data_valid = false;
#endif

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


	/*
	 * Debug helper to compare output with the output of another time integrator
	 */
	void print_debug(
			const char *name
	)
	{
		SphereData tmp = *this;
		std::cout << name << ":" << std::endl;
		std::cout << "                min: " << tmp.physical_reduce_min() << std::endl;
		std::cout << "                max: " << tmp.physical_reduce_max() << std::endl;
		std::cout << "                sum: " << tmp.physical_reduce_sum() << std::endl;
		std::cout << "                suminc: " << tmp.physical_reduce_debug_sum_quad_mul_increasing() << std::endl;
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
SphereData operator*(
		double i_value,
		const SphereData &i_array_data
)
{
	return ((SphereData&)i_array_data)*i_value;
}




#endif /* SPHDATA_HPP_ */
