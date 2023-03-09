/*
 * SphereData.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
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
#include <sweet/core/parmemcpy.hpp>
#include <sweet/core/MemBlockAlloc.hpp>
#include <sweet/core/openmp_helper.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>
#include <sweet/core/sphere/SphereData_PhysicalComplex.hpp>
#include <sweet/core/SWEETError.hpp>



namespace sweet
{

class SphereData_Spectral
{
	friend class SphereData_SpectralComplex;

	typedef std::complex<double> Tcomplex;


public:
	const SphereData_Config *sphereDataConfig = nullptr;

public:
	std::complex<double> *spectral_space_data = nullptr;


public:
	std::complex<double>& operator[](std::size_t i)
	{
		return spectral_space_data[i];
	}

	const std::complex<double>& operator[](std::size_t i)	const
	{
		return spectral_space_data[i];
	}

public:
	void swap(
			SphereData_Spectral &i_sphereData
	)
	{
		assert(sphereDataConfig == i_sphereData.sphereDataConfig);

		std::swap(spectral_space_data, i_sphereData.spectral_space_data);
	}

	void swapWithConfig(
			SphereData_Spectral &i_sphereData
	)
	{
		std::swap(sphereDataConfig, i_sphereData.sphereDataConfig);

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


#if 0
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
#endif

#if 1
public:
	SphereData_Spectral(
			const SphereData_Config *i_sphereDataConfig,
			const double &i_value
	)	:
		sphereDataConfig(i_sphereDataConfig),
		spectral_space_data(nullptr)
	{
		assert(i_sphereDataConfig != 0);

		setup(i_sphereDataConfig);
		spectral_set_value(i_value);
	}
#endif


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

		std::swap(spectral_space_data, i_sph_data.spectral_space_data);
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
	bool setup(
		const SphereData_Config *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		return alloc_data();
	}


public:
	bool setup(
		const SphereData_Config &i_sphereDataConfig
	)
	{
		return setup(&i_sphereDataConfig);
	}


public:
	bool setup(
		const SphereData_Config *i_sphereDataConfig,
		double i_value
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		bool retval = alloc_data();
		spectral_set_value(i_value);

		return retval;
	}


private:
	bool alloc_data()
	{
		assert(spectral_space_data == nullptr);
		spectral_space_data = MemBlockAlloc::alloc<Tcomplex>(sphereDataConfig->spectral_array_data_number_of_elements * sizeof(Tcomplex));
		return true;
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
	void clear()
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
		clear();
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
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b0 D^2 + b1 D^4 + b2 D^6 + b3 D^8) x = rhs
	 */
	inline
	SphereData_Spectral spectral_solve_helmholtz_higher_order(
			const double & a,
			const std::array<double, 4> & b,
			double r
	)
	{
		SphereData_Spectral out(*this);

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				double laplace_op_2 = - (double)n*((double)n+1.0)/(r*r);
				double laplace_op_4 = laplace_op_2 * laplace_op_2;
				double laplace_op_6 = laplace_op_4 * laplace_op_2;
				double laplace_op_8 = laplace_op_4 * laplace_op_4;
				io_data /= (a + b[0] * laplace_op_2 + b[1] * laplace_op_4 + b[2] * laplace_op_6 + b[3] * laplace_op_8);
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
	 * return the sum of all values, use quad precision for reduction
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
	 * return the sum of squares of all values, use quad precision for reduction
	 * Important: Since m=0  modes appear only once and m>0 appear twice in the full spectrum
	 */
	double spectral_reduce_sum_sqr_quad()	const
	{
		std::complex<double> sum = 0;
		//std::complex<double> c = 0;

		//m=0 case - weight 1
		std::size_t idx = sphereDataConfig->getArrayIndexByModes(0, 0);
		for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
		{
			sum += spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
			idx++;
		}
		
		//m>0 case - weight 2, as they appear twice
		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				sum += 2.0*spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
				idx++;
			}
		}

		return sum.real();
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
	 * Return the max value
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


	/**
	 * Return the max abs value
	 */
	double spectral_reduce_max_abs()	const
	{
		double error = -std::numeric_limits<double>::infinity();
		std::complex<double> w = {0,0};
		for (int j = 0; j < sphereDataConfig->spectral_array_data_number_of_elements; j++)
		{
			w = spectral_space_data[j]*std::conj(spectral_space_data[j]);
			error = std::max(std::abs(w), error);
		}

		return error;
	}


	/**
	 * Return the max abs value for the first rnorm spectral coefficients
	 */
	double spectral_reduce_max_abs(std::size_t rnorm)	const
	{

		assert ((int)rnorm <= sphereDataConfig->spectral_modes_m_max);
		assert ((int)rnorm <= sphereDataConfig->spectral_modes_n_max);

		double error = -std::numeric_limits<double>::infinity();
		std::complex<double> w = {0,0};

		for (std::size_t m = 0; m <= rnorm; m++)
		{
			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (std::size_t n = m; n <= rnorm; n++)
			{
				w = spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
				error = std::max(std::abs(w), error);
				idx++;
			}
		}

		return error;
	}



	/**
	 * Return the minimum abs value
	 */
	double spectral_reduce_min_abs()	const
	{
		double error = std::numeric_limits<double>::infinity();
		std::complex<double> w = {0,0};
		for (int j = 0; j < sphereDataConfig->spectral_array_data_number_of_elements; j++)
		{
			w = spectral_space_data[j]*std::conj(spectral_space_data[j]);
			error = std::min(std::abs(w), error);
		}

		return error;
	}

	bool spectral_reduce_is_any_nan_or_inf()	const
	{
		bool retval = false;

#if SWEET_THREADING_SPACE
		#pragma omp parallel for simd reduction(|:retval)
#endif
		for (int j = 0; j < sphereDataConfig->spectral_array_data_number_of_elements; j++)
		{
			retval |= std::isnan(spectral_space_data[j].real());
			retval |= std::isinf(spectral_space_data[j].real());
			retval |= std::isnan(spectral_space_data[j].imag());
			retval |= std::isinf(spectral_space_data[j].imag());
		}

		return retval;
	}

	bool spectral_is_first_nan_or_inf()	const
	{
		bool retval = false;

		retval |= std::isnan(spectral_space_data[0].real());
		retval |= std::isinf(spectral_space_data[0].real());
		retval |= std::isnan(spectral_space_data[0].imag());
		retval |= std::isinf(spectral_space_data[0].imag());

		return retval;
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

	void spectral_structure_print(
			int i_precision = 16,
			double i_abs_threshold = -1
	)	const
	{
		std::cout << std::setprecision(i_precision);

		std::cout << "m \\ n ----->" << std::endl;
		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			//std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				std::cout << "(" << m << "," << n << ")" << "\t";
			}
			std::cout << std::endl;
		}
	}


	/**
	 * Write spectral data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_spectral_saveData_ascii(
			const char *i_filename,		///< Name of file to store data to
			char i_separator = '\t',	///< separator to use for each line
			int i_precision = 12,		///< number of floating point digits
			int dimension = 2			///< store 1D or 2D
	)	const
	{

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		file << "#SWEET_SPHERE_SPECTRAL_CPLX_DATA_ASCII" << std::endl;


		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				const std::complex<double> &value = spectral_space_data[idx];
				file << "(" << value.real() << ", " << value.imag() << ")";
				idx++;

				if (n < sphereDataConfig->spectral_modes_n_max)
					file << i_separator;
				else
					file << std::endl;

			}
		}

		return true;
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


	void spectrum_abs_file_write_line(
			const std::string &i_filename,
			const char *i_title = "",
			const double i_time = 0.0,
			int i_precision = 20,
			double i_abs_threshold = -1,
			int i_reduce_mode_factor = 1
	)	const
	{

		std::ofstream file;

		if (i_time == 0.0){
			file.open(i_filename, std::ios_base::trunc);
			file << std::setprecision(i_precision);
			file << "#SWEET_SPHERE_SPECTRAL_ABS_EVOL_ASCII" << std::endl;
  			file << "#TI " << i_title << std::endl;
			file << "0\t" << std::endl; // Use 0 to make it processable by python
			file << "(n_max=" <<sphereDataConfig->spectral_modes_n_max << " m_max="
					<< sphereDataConfig->spectral_modes_n_max << ")" << std::endl;
			file << "timestamp\t" ; 
			for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
			{
				//std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
				for (int n = m; n <= sphereDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
				{
					file << "(" << n << ";" << m << ")\t" ;
				}
			}
			file<< "SpectralSum" <<std::endl;
		}
		else{
			file.open(i_filename, std::ios_base::app);
			file << std::setprecision(i_precision);
		}  		

  		std::complex<double> w = {0,0};
		double wabs = 0.0;
		double sum = 0.0;
		//std::cout << "n" << " " << "m" << " " << "norm" <<std::endl;
		file << i_time << "\t";
  		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
  		{
  			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
  			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
  			{
  				w = spectral_space_data[idx];
				wabs = std::abs(w * std::conj(w));
				if ( m > 0 ) sum += 2*wabs; //term appears twice in the spectrum
				else sum += wabs;      // term appears only once
  				
				if ( wabs < i_abs_threshold){
					//file << "(" << n << "," << m << ")\t" <<std::endl;
					file <<  0 << "\t"; //<<std::endl;
				}
				else{
					//file << "(" << n << "," << m << ")\t" <<std::endl;
					file <<  wabs << "\t"; //<<std::endl;;
					//std::cout << n << " " << m << " " << wabs <<std::endl;
				}
				idx++;
  			}
  		}
		file<< sum << std::endl;
  		file.close();
	}


	void spectrum_phase_file_write_line(
			const std::string &i_filename,
			const char *i_title = "",
			const double i_time = 0.0,
			int i_precision = 20,
			double i_abs_threshold = -1,
			int i_reduce_mode_factor = 1
	)	const
	{

		std::ofstream file;

		if (i_time == 0.0){
			file.open(i_filename, std::ios_base::trunc);
			file << std::setprecision(i_precision);
			file << "#SWEET_SPHERE_SPECTRAL_PHASE_EVOL_ASCII" << std::endl;
  			file << "#TI " << i_title << std::endl;
			file << "0\t" << std::endl; // Use 0 to make it processable by python
			file << "(n_max=" <<sphereDataConfig->spectral_modes_n_max << " m_max="
					<< sphereDataConfig->spectral_modes_n_max << ")" << std::endl;
			file << "timestamp\t" ; 
			for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
			{
				//std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
				for (int n = m; n <= sphereDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
				{
					file << "(" << n << ";" << m << ")\t" ;
				}
			}
			file << std::endl;
		}
		else{
			file.open(i_filename, std::ios_base::app);
			file << std::setprecision(i_precision);
		}  		

  		std::complex<double> w = {0,0};
		double wphase = 0.0;
		
		//std::cout << "n" << " " << "m" << " " << "norm" <<std::endl;
		file << i_time << "\t";
  		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max/i_reduce_mode_factor; m++)
  		{
  			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
  			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max/i_reduce_mode_factor; n++)
  			{
  				w = spectral_space_data[idx];
				wphase = std::arg(w); // std::abs(w * std::conj(w));
				
				//file << "(" << n << "," << m << ")\t" <<std::endl;
				file <<  wphase << "\t"; //<<std::endl;;
				//std::cout << n << " " << m << " " << wabs <<std::endl;
				
				idx++;
  			}
  		}
		file<< std::endl;
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

  		file.write((const char*)spectral_space_data, sizeof(std::complex<double>)*sphereDataConfig->spectral_array_data_number_of_elements);

  		file.close();
	}


  	void file_read_binary_spectral(
			const std::string &i_filename
	)
	{
  		std::ifstream file(i_filename, std::ios_base::binary);

		if (!file.is_open())
			SWEETError("Error while opening file " + i_filename);

  		std::string magic;
  		std::getline(file, magic);

  		if (magic != "SWEET")
  			SWEETError("Magic code 'SWEET' not found");

  		std::string data_type;
  		int num_lon = -1;
  		int num_lat = -1;
  		//int size = -1;
                std::string grid_type = "";

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
  				///std::cout << data_type << std::endl;
  	  			continue;
  			}

  			//if (buf == "NUM_LON")
  			if (buf == "MODES_M_MAX")
  			{
  				file >> buf;
  				num_lon = std::stoi(buf);
  				///std::cout << num_lon << std::endl;
  				continue;
  			}

  			if (buf == "MODES_N_MAX")
  			////if (buf == "NUM_LAT")
  			{
  				file >> buf;
  				num_lat = std::stoi(buf);
  				///std::cout << num_lat << std::endl;
  				continue;
  			}

  			if (buf == "GRID_TYPE")
  			//if (buf == "SIZE")
  			{
  				file >> buf;
  				grid_type = buf;
  				///std::cout << grid_type << std::endl;
  				continue;
  			}

  			if (buf == "NUM_ELEMENTS")
  			//if (buf == "SIZE")
  			{
  				file >> buf;
  				//size = std::stoi(buf);
  				//std::cout << size << std::endl;
  				continue;
  			}

  			SWEETError("Unknown Tag '"+buf+"'");
  		}

  		// read last newline
  		char nl;
  		file.read(&nl, 1);
  		///std::cout << file.tellg() << std::endl;

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


	/**
	 * Interpolate from a finer mesh. Remove highest frequency modes.
	 */
	SphereData_Spectral restrict(
			const SphereData_Spectral &i_array_data
	)
	{
		SphereData_Spectral out = *this;
		out = i_array_data.spectral_returnWithDifferentModes(out.sphereDataConfig);
		return out;

	}


	/**
	 * Interpolate from a coarser mesh. Pad zeros corresponding to highest frequency modes.
	 */
	SphereData_Spectral pad_zeros(
			const SphereData_Spectral &i_array_data
	)
	{

		SphereData_Spectral out = *this;
		out = i_array_data.spectral_returnWithDifferentModes(out.sphereDataConfig);
		return out;

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

}


#endif
