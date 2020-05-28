/*
 * SphereData.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SWEET_SPHERE_DATA_SPECTRAL_HPP_
#define SWEET_SPHERE_DATA_SPECTRAL_HPP_

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
#include <sweet/openmp_helper.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/SphereData_PhysicalComplex.hpp>
#include <sweet/SWEETError.hpp>



class SphereData_Spectral
{
	friend class SphereData_SpectralComplex;

	typedef std::complex<double> Tcomplex;


public:
	const SphereData_Config *sphereDataConfig = nullptr;

public:
	std::complex<double> *spectral_space_data = nullptr;

	std::complex<double>& operator[](std::size_t i)
	{
		return spectral_space_data[i];
	}

	const std::complex<double>& operator[](std::size_t i)	const
	{
		return spectral_space_data[i];
	}

	void swap(
			SphereData_Spectral &i_sphereData
	)
	{
		assert(sphereDataConfig == i_sphereData.sphereDataConfig);

		std::swap(spectral_space_data, i_sphereData.spectral_space_data);
	}


public:
	SphereData_Spectral(
			const SphereData_Config *i_sphereDataConfig
	)	:
		sphereDataConfig(i_sphereDataConfig),
		spectral_space_data(nullptr)
	{
		assert(i_sphereDataConfig != 0);

		setup(i_sphereDataConfig);
	}


public:
	SphereData_Spectral(
			const SphereData_Config *i_sphereDataConfig,
			const std::complex<double> &i_value
	)	:
		sphereDataConfig(i_sphereDataConfig),
		spectral_space_data(nullptr)
	{
		assert(i_sphereDataConfig != 0);

		setup(i_sphereDataConfig);
		spectral_set_value(i_value);
	}



public:
	SphereData_Spectral(
			const SphereData_Config *i_sphereDataConfig,
			double &i_value
	)	:
		sphereDataConfig(i_sphereDataConfig),
		spectral_space_data(nullptr)
	{
		assert(i_sphereDataConfig != 0);

		setup(i_sphereDataConfig);
		spectral_set_value(i_value);
	}



public:
	SphereData_Spectral()	:
		sphereDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
	}



public:
	SphereData_Spectral(
			const SphereData_Spectral &i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		spectral_space_data(nullptr)
	{
		if (i_sph_data.sphereDataConfig == nullptr)
			return;

		alloc_data();

		operator=(i_sph_data);
	}



public:
	SphereData_Spectral(
			SphereData_Spectral &&i_sph_data
	)	:
		sphereDataConfig(i_sph_data.sphereDataConfig),
		spectral_space_data(nullptr)
	{
		if (i_sph_data.sphereDataConfig == nullptr)
			return;

		alloc_data();

		spectral_space_data = i_sph_data.spectral_space_data;
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

		assert(sphereDataConfig->spectral_modes_m_max == i_sphereDataConfig->spectral_modes_m_max);
		assert(sphereDataConfig->spectral_modes_n_max == i_sphereDataConfig->spectral_modes_n_max);
	}


public:
	SphereData_Spectral& operator=(
			const SphereData_Spectral &i_sph_data
	)
	{
		if (i_sph_data.sphereDataConfig == nullptr)
			return *this;

		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(Tcomplex)*sphereDataConfig->spectral_array_data_number_of_elements);

		return *this;
	}




	/**
	 * This function implements copying the spectral data only.
	 *
	 * This becomes handy if coping with data which should be only transformed without dealiasing.
	 */
public:
	SphereData_Spectral& load_nodealiasing(
			const SphereData_Spectral &i_sph_data		///< data to be converted to sphereDataConfig_nodealiasing
	)
	{
		if (sphereDataConfig == nullptr)
			SWEETError("sphereDataConfig not initialized");

		parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(Tcomplex)*sphereDataConfig->spectral_array_data_number_of_elements);

		return *this;
	}


public:
	SphereData_Spectral& operator=(
			SphereData_Spectral &&i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		std::swap(spectral_space_data, i_sph_data.spectral_space_data);

		return *this;
	}


public:
	SphereData_Spectral spectral_returnWithDifferentModes(
			const SphereData_Config *i_sphereDataConfig
	)	const
	{
		SphereData_Spectral out(i_sphereDataConfig);

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

		if (scaling_mode == -1)
		{
			/*
			 * more modes -> less modes
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR
			for (int m = 0; m <= out.sphereDataConfig->spectral_modes_m_max; m++)
			{
				Tcomplex *dst = &out.spectral_space_data[out.sphereDataConfig->getArrayIndexByModes(m, m)];
				Tcomplex *src = &spectral_space_data[sphereDataConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(Tcomplex)*(out.sphereDataConfig->spectral_modes_n_max-m+1);
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
				Tcomplex *dst = &out.spectral_space_data[out.sphereDataConfig->getArrayIndexByModes(m, m)];
				Tcomplex *src = &spectral_space_data[sphereDataConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(Tcomplex)*(sphereDataConfig->spectral_modes_n_max-m+1);
				memcpy(dst, src, size);
			}
		}

		return out;
	}



	/*
	 * Setup spectral sphere data based on data in physical space
	 */
	void loadSphereDataPhysical(
			const SphereData_Physical &i_sphereDataPhysical
	)
	{
		/**
		 * Warning: The sphat_to_SH function is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 * Hence, we create a copy
		 */
		SphereData_Physical tmp(i_sphereDataPhysical);
		spat_to_SH(sphereDataConfig->shtns, tmp.physical_space_data, spectral_space_data);
	}



	/*
	 * Return the data converted to physical space
	 */
	SphereData_Physical getSphereDataPhysical()	const
	{
		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		SphereData_Spectral tmp(*this);
		SphereData_Physical retval(sphereDataConfig);
		SH_to_spat(sphereDataConfig->shtns, tmp.spectral_space_data, retval.physical_space_data);

		return retval;
	}


	/*
	 * Return the data converted to physical space
	 *
	 * alias for "getSphereDataPhysical"
	 */
	SphereData_Physical toPhys()	const
	{
		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		SphereData_Spectral tmp(*this);
		SphereData_Physical retval(sphereDataConfig);
		SH_to_spat(sphereDataConfig->shtns, tmp.spectral_space_data, retval.physical_space_data);

		return retval;
	}


	SphereData_PhysicalComplex getSphereDataPhysicalComplex()	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		/*
		 * WARNING:
		 * We have to use a temporary array here because of destructive SH transformations
		 */
		SphereData_Spectral tmp_spectral(*this);
		SphereData_Physical tmp_physical(sphereDataConfig);
		SH_to_spat(sphereDataConfig->shtns, tmp_spectral.spectral_space_data, tmp_physical.physical_space_data);

		parmemcpy(out.physical_space_data, tmp_physical.physical_space_data, sizeof(double)*sphereDataConfig->physical_array_data_number_of_elements);

		return out;
	}



	SphereData_Spectral(
			const SphereData_Physical &i_sphere_data_physical
	)
	{
		setup(i_sphere_data_physical.sphereDataConfig);

		loadSphereDataPhysical(i_sphere_data_physical);
	}



	SphereData_Spectral operator+(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_Spectral out(sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

		return out;
	}



	SphereData_Spectral& operator+=(
			const SphereData_Spectral &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] += i_sph_data.spectral_space_data[idx];

		return *this;
	}


	SphereData_Spectral& operator-=(
			const SphereData_Spectral &i_sph_data
	)
	{
		check(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

		return *this;
	}



	SphereData_Spectral operator-(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_Spectral out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

		return out;
	}



	SphereData_Spectral operator-()	const
	{
		SphereData_Spectral out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = -spectral_space_data[idx];

		return out;
	}



	SphereData_Spectral operator*(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_Physical a = getSphereDataPhysical();
		SphereData_Physical b = i_sph_data.toPhys();

		SphereData_Physical mul(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			mul.physical_space_data[i] = a.physical_space_data[i]*b.physical_space_data[i];

		SphereData_Spectral out(mul);

		return out;
	}



	SphereData_Spectral operator/(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		check(i_sph_data.sphereDataConfig);

		SphereData_Physical a = getSphereDataPhysical();
		SphereData_Physical b = i_sph_data.toPhys();

		SphereData_Physical div(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < sphereDataConfig->physical_array_data_number_of_elements; i++)
			div.physical_space_data[i] = a.physical_space_data[i]/b.physical_space_data[i];

		SphereData_Spectral out(div);

		return out;
	}



	SphereData_Spectral operator*(
			double i_value
	)	const
	{
		SphereData_Spectral out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		return out;
	}




	const SphereData_Spectral& operator*=(
			double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}




	const SphereData_Spectral& operator/=(
			double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] /= i_value;

		return *this;
	}





	const SphereData_Spectral& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}



	SphereData_Spectral operator/(
			double i_value
	)	const
	{
		SphereData_Spectral out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		return out;
	}


	SphereData_Spectral operator+(
			double i_value
	)	const
	{
		SphereData_Spectral out(*this);
		out.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
		return out;
	}


	SphereData_Spectral operator_scalar_sub_this(
			double i_value
	)	const
	{
		SphereData_Spectral out(sphereDataConfig);
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int idx = 0; idx < sphereDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = -spectral_space_data[idx];

		out.spectral_space_data[0] = i_value*std::sqrt(4.0*M_PI) + out.spectral_space_data[0];
		return out;
	}



	const SphereData_Spectral& operator+=(
			double i_value
	)	const
	{
		spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
		return *this;
	}




	SphereData_Spectral operator-(
			double i_value
	)	const
	{
		SphereData_Spectral out(*this);
		out.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);
		return out;
	}



	SphereData_Spectral& operator-=(
			double i_value
	)
	{
		spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);
		return *this;
	}


public:
	void setup(
		const SphereData_Config *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		alloc_data();
	}


public:
	void setup(
		const SphereData_Config *i_sphereDataConfig,
		double i_value
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		alloc_data();
		spectral_set_value(i_value);
	}


private:
	void alloc_data()
	{
		assert(spectral_space_data == nullptr);
		spectral_space_data = MemBlockAlloc::alloc<Tcomplex>(sphereDataConfig->spectral_array_data_number_of_elements * sizeof(Tcomplex));
	}


public:
	void setup_if_required(
		const SphereData_Config *i_sphereDataConfig
	)
	{
		if (sphereDataConfig != nullptr)
			return;

		setup(i_sphereDataConfig);
	}



public:
	void free()
	{
		if (spectral_space_data != nullptr)
		{
			MemBlockAlloc::free(spectral_space_data, sphereDataConfig->spectral_array_data_number_of_elements * sizeof(Tcomplex));
			spectral_space_data = nullptr;

			sphereDataConfig = nullptr;
		}
	}



public:
	~SphereData_Spectral()
	{
		free();
	}



	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	inline
	SphereData_Spectral spectral_solve_helmholtz(
			const double &i_a,
			const double &i_b,
			double r
	)
	{
		SphereData_Spectral out(*this);

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
	SphereData_Spectral spectral_solve_laplace(
			double r
	)
	{
		SphereData_Spectral out(*this);

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



	/**
	 * Truncate modes which are not representable in spectral space
	 */
	const SphereData_Spectral& spectral_truncate()	const
	{
		SphereData_Physical tmp(sphereDataConfig);

		SH_to_spat(sphereDataConfig->shtns, spectral_space_data, tmp.physical_space_data);
		spat_to_SH(sphereDataConfig->shtns, tmp.physical_space_data, spectral_space_data);

		return *this;
	}


	void spectral_update_lambda(
			std::function<void(int,int,Tcomplex&)> i_lambda
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				i_lambda(n, m, spectral_space_data[idx]);
				idx++;
			}
		}
	}

	const std::complex<double>& spectral_get_DEPRECATED(
			int i_n,
			int i_m
	)	const
	{
		static const std::complex<double> zero = {0,0};

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


	const std::complex<double>& spectral_get_(
			int i_n,
			int i_m
	)	const
	{
		assert(i_n >= 0 && i_m >= 0);
		assert(i_n <= sphereDataConfig->spectral_modes_n_max);
		assert(i_m <= sphereDataConfig->spectral_modes_m_max);
		assert(i_m <= i_n);

		return spectral_space_data[sphereDataConfig->getArrayIndexByModes(i_n, i_m)];
	}



	void spectral_set(
			int i_n,
			int i_m,
			const std::complex<double> &i_data
	)	const
	{
#if SWEET_DEBUG
		if (i_n < 0 ||  i_m < 0)
			SWEETError("Out of boundary a");

		if (i_n > sphereDataConfig->spectral_modes_n_max)
			SWEETError("Out of boundary b");

		if (i_m > sphereDataConfig->spectral_modes_m_max)
			SWEETError("Out of boundary c");

		if (i_m > i_n)
			SWEETError("Out of boundary d");

		assert (i_m <= sphereDataConfig->spectral_modes_m_max);
#endif

		spectral_space_data[sphereDataConfig->getArrayIndexByModes(i_n, i_m)] = i_data;
	}


	/*
	 * Set all values to zero
	 */
	void spectral_set_zero()
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < sphereDataConfig->spectral_array_data_number_of_elements; i++)
			spectral_space_data[i] = {0,0};
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
	}



	/*
	 * Add a constant in physical space by adding the corresponding value in spectral space
	 */
	void spectral_add_physical_constant(
			double i_value
	)	const
	{
		this->spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
	}



	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	std::complex<double> spectral_reduce_sum_quad_increasing()	const
	{
		std::complex<double> sum = 0;
		std::complex<double> c = 0;
#if SWEET_THREADING_SPACE
//#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
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
		std::complex<double> sum = 0;
		std::complex<double> c = 0;
#if SWEET_THREADING_SPACE
//#pragma omp parallel for PROC_BIND_CLOSE reduction(+:sum,c)
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
		std::complex<double> error = -std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphereDataConfig->spectral_array_data_number_of_elements; j++)
		{
			error.real(std::max(spectral_space_data[j].real(), error.real()));
			error.imag(std::max(spectral_space_data[j].imag(), error.imag()));
		}

		return error;
	}


	void spectral_print(
			int i_precision = 16,
			double i_abs_threshold = -1
	)	const
	{
		std::cout << std::setprecision(i_precision);

		std::cout << "m \\ n ----->" << std::endl;
		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				if (std::abs(spectral_space_data[idx]) < i_abs_threshold)
					std::cout << 0 << "\t";
				else
					std::cout << spectral_space_data[idx] << "\t";
				idx++;
			}
			std::cout << std::endl;
		}
	}


  
  	void spectrum_file_write(
			const std::string &i_filename,
			const char *i_title = "",
			int i_precision = 20
	)	const
	{
  		std::ofstream file(i_filename, std::ios_base::trunc);

  		file << std::setprecision(i_precision);
		file << "#SWEET_SPHERE_SPECTRAL_DATA_ASCII" << std::endl;

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
  		std::ofstream file(i_filename, std::ios_base::trunc | std::ios_base::binary);

		if (!file.is_open())
			SWEETError("Error while opening file");

  		file << "SWEET" << std::endl;
  		file << "DATA_TYPE SH_DATA" << std::endl;
  		file << "MODES_M_MAX " << sphereDataConfig->spectral_modes_m_max << std::endl;
  		file << "MODES_N_MAX " << sphereDataConfig->spectral_modes_n_max << std::endl;
  		file << "GRID_TYPE GAUSSIAN" << std::endl;
  		file << "NUM_ELEMENTS " << sphereDataConfig->spectral_array_data_number_of_elements << std::endl;
  		file << "FIN" << std::endl;
  		//std::cout << file.tellp() << std::endl;

  		file.write((const char*)spectral_space_data, sizeof(std::complex<double>)*sphereDataConfig->spectral_array_data_number_of_elements);

  		file.close();
	}


  	void file_read_binary_spectral(
			const std::string &i_filename
	)
	{
  		std::ifstream file(i_filename, std::ios_base::binary);

		if (!file.is_open())
			SWEETError("Error while opening file");

  		std::string magic;
  		std::getline(file, magic);

  		if (magic != "SWEET")
  			SWEETError("Magic code 'SWEET' not found");

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

  			SWEETError("Unknown Tag '"+buf+"'");
  		}

  		// read last newline
  		char nl;
  		file.read(&nl, 1);
  		std::cout << file.tellg() << std::endl;

  		if (data_type != "SH_DATA")
  			SWEETError("Unknown data type '"+data_type+"'");

  		if (num_lon != sphereDataConfig->spectral_modes_m_max)
  			SWEETError("NUM_LON "+std::to_string(num_lon)+" doesn't match SphereDataConfig");

  		if (num_lat != sphereDataConfig->spectral_modes_n_max)
  			SWEETError("NUM_LAT "+std::to_string(num_lat)+" doesn't match SphereDataConfig");

  		file.read((char*)spectral_space_data, sizeof(std::complex<double>)*sphereDataConfig->spectral_array_data_number_of_elements);

  		file.close();
	}




	void normalize(
			const std::string &normalization = ""
	)
	{
		if (normalization == "avg_zero")
		{
			// move average value to 0
			double phi_min = getSphereDataPhysical().physical_reduce_min();
			double phi_max = getSphereDataPhysical().physical_reduce_max();

			double avg = 0.5*(phi_max+phi_min);

			operator-=(avg);
		}
		else if (normalization == "min_zero")
		{
			// move minimum value to zero
			double phi_min = getSphereDataPhysical().physical_reduce_min();
			operator-=(phi_min);
		}
		else if (normalization == "max_zero")
		{
			// move maximum value to zero
			double phi_max = getSphereDataPhysical().physical_reduce_max();
			operator-=(phi_max);
		}
		else if (normalization == "")
		{
		}
		else
		{
			SWEETError("Normalization not supported!");
		}
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
SphereData_Spectral operator*(
		double i_value,
		const SphereData_Spectral &i_array_data
)
{
	return ((SphereData_Spectral&)i_array_data)*i_value;
}






/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData
 *
 */
inline
static
SphereData_Spectral operator+(
		double i_value,
		const SphereData_Spectral &i_array_data
)
{
	return ((SphereData_Spectral&)i_array_data)+i_value;
}




/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData
 *
 */
inline
static
SphereData_Spectral operator-(
		double i_value,
		const SphereData_Spectral &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}




#endif /* SWEET_SPHERE_DATA_HPP_ */
