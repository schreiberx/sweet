/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_SPHERE_DATA_SPECTRAL_COMPLEX_HPP_
#define SRC_SPHERE_DATA_SPECTRAL_COMPLEX_HPP_

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

	typedef std::complex<double> Tcomplex;

public:
	std::complex<double> *spectral_space_data;

	std::complex<double>& operator[](std::size_t i)
	{
		return spectral_space_data[i];
	}

	const std::complex<double>& operator[](std::size_t i)	const
	{
		return spectral_space_data[i];
	}


public:
	SphereData_SpectralComplex()	:
		sphereDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
	}


public:
	SphereData_SpectralComplex(
			const SphereData_Config *i_sphereDataConfig
	)	:
		sphereDataConfig(nullptr),
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
		spectral_space_data(nullptr)
	{
		setup(i_sph_data.sphereDataConfig);

		operator=(i_sph_data);
	}




public:
	SphereData_SpectralComplex(
			SphereData_SpectralComplex &&i_data
	)	:
		sphereDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
		if (sphereDataConfig == nullptr)
			setup(i_data.sphereDataConfig);

		assert(i_data.spectral_space_data != nullptr);

		std::swap(spectral_space_data, i_data.spectral_space_data);
	}





	SphereData_SpectralComplex(
			const SphereData_PhysicalComplex &i_sph_data
	):
		sphereDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
		SphereData_PhysicalComplex tmp = i_sph_data;

		setup(i_sph_data.sphereDataConfig);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		spat_cplx_to_SH(sphereDataConfig->shtns, tmp.physical_space_data, spectral_space_data);
	}




	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check_sphereDataConfig_identical_res(const SphereData_Config *i_sphereDataConfig)	const
	{
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

		assert(i_sph_data.spectral_space_data);
		parmemcpy(spectral_space_data, i_sph_data.spectral_space_data, sizeof(Tcomplex)*sphereDataConfig->spectral_complex_array_data_number_of_elements);

		return *this;
	}




public:
	SphereData_SpectralComplex& operator=(
			SphereData_SpectralComplex &&i_sph_data
	)
	{
		if (sphereDataConfig == nullptr)
			setup(i_sph_data.sphereDataConfig);

		assert(i_sph_data.spectral_space_data);
		std::swap(spectral_space_data, i_sph_data.spectral_space_data);

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

		return out;
	}


	/*
	 * Setup spectral sphere data based on data in physical space
	 */
	void loadSphereDataPhysical(
			const SphereData_PhysicalComplex &i_sphereDataPhysical
	)
	{
		/**
		 * Warning: The sphat_to_SH function is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 * Hence, we create a copy
		 */
		SphereData_PhysicalComplex tmp(i_sphereDataPhysical);
		spat_cplx_to_SH(sphereDataConfig->shtns, tmp.physical_space_data, spectral_space_data);
	}



	SphereData_PhysicalComplex toPhys()	const
	{
		return getSphereDataPhysicalComplex();
	}

	SphereData_PhysicalComplex getSphereDataPhysicalComplex()	const
	{
		SphereData_PhysicalComplex out(sphereDataConfig);

		/*
		 * WARNING:
		 * We have to use a temporary array here because of destructive SH transformations
		 */
		SphereData_SpectralComplex tmp = *this;
		SH_to_spat_cplx(sphereDataConfig->shtns, tmp.spectral_space_data, out.physical_space_data);

		return out;
	}


#if 0
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
#endif

	SphereData_SpectralComplex operator+(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] + i_sph_data.spectral_space_data[idx];

		return out_sph_data;
	}



	SphereData_SpectralComplex& operator+=(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] += i_sph_data.spectral_space_data[idx];

		return *this;
	}


	SphereData_SpectralComplex& operator-=(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_sph_data.spectral_space_data[idx];

		return *this;
	}



	SphereData_SpectralComplex operator-(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		check_sphereDataConfig_identical_res(i_sph_data.sphereDataConfig);

		SphereData_SpectralComplex out_sph_data(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx] - i_sph_data.spectral_space_data[idx];

		return out_sph_data;
	}


	SphereData_SpectralComplex operator-()	const
	{
		SphereData_SpectralComplex out_sph_data(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = -spectral_space_data[idx];

		return out_sph_data;
	}




	SphereData_SpectralComplex operator*(
			const double i_value
	)	const
	{
		SphereData_SpectralComplex out_sph_data(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR

		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		return out_sph_data;
	}



	const SphereData_SpectralComplex& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}


	SphereData_SpectralComplex operator/(
			double i_value
	)	const
	{
		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		return out_sph_data;
	}


	const SphereData_SpectralComplex& operator/=(
			const std::complex<double> &i_value
	)	const
	{

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] /= i_value;

		return *this;
	}



	SphereData_SpectralComplex operator+(
			double i_value
	)	const
	{
		SphereData_SpectralComplex out_sph_data(*this);

		out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		return out_sph_data;
	}



	SphereData_SpectralComplex& operator+=(
			double i_value
	)
	{
		spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
		return *this;
	}



	SphereData_SpectralComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		SphereData_SpectralComplex out_sph_data(*this);

		out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		return out_sph_data;
	}



	SphereData_SpectralComplex operator-(
			const std::complex<double> &i_value
	)	const
	{
		SphereData_SpectralComplex out_sph_data(*this);

		out_sph_data.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);

		return out_sph_data;
	}



	SphereData_SpectralComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		SphereData_SpectralComplex out_sph_data(sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_sph_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		return out_sph_data;
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
	void setup(
			const SphereData_Config *i_sphereConfig
	)
	{
		// assure that the initialization is not done twice!
		assert(sphereDataConfig == nullptr);

		sphereDataConfig = i_sphereConfig;

		spectral_space_data = MemBlockAlloc::alloc<Tcomplex>(sphereDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
	}


public:
	~SphereData_SpectralComplex()
	{
		MemBlockAlloc::free(spectral_space_data, sphereDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
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



	inline
	void spectral_update_lambda(
			std::function<void(int,int,Tcomplex&)> i_lambda
	)
	{
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



#if 0
	inline
	const std::complex<double>& spectral_get_DEPRECATED(
			int in,
			int im
	)	const
	{
		static const std::complex<double> zero = {0,0};

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
#endif



	inline
	const std::complex<double>& spectral_get_(
			int in,
			int im
	)	const
	{
		assert(in >= 0);
		assert(in <= sphereDataConfig->spectral_modes_n_max);
		assert(std::abs(im) <= sphereDataConfig->spectral_modes_m_max);
		assert(std::abs(im) <= in);

		return spectral_space_data[sphereDataConfig->getArrayIndexByModes_Complex(in, im)];
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
	}


	void spectral_print(
			int i_precision = 8
	)	const
	{
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

};




inline
static
SphereData_SpectralComplex operator*(
		const double i_value,
		const SphereData_SpectralComplex &i_array_data
)
{
	return i_array_data*i_value;
}


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
	SphereData_SpectralComplex out_sph_data(i_array_data.sphereDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < i_array_data.sphereDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_sph_data.spectral_space_data[idx] = -i_array_data.spectral_space_data[idx];

	out_sph_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_sph_data;

}

#endif
