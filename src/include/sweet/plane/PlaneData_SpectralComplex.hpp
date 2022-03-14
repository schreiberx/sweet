/*
 * PlaneData_SpectralComplex.hpp
 *
 *  Created on: 14 Mar 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PLANE_DATA_SPECTRAL_COMPLEX_HPP_
#define SRC_PLANE_DATA_SPECTRAL_COMPLEX_HPP_

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
#include <sweet/plane/PlaneDataConfig.hpp>
#include <sweet/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>

	#if SWEET_THREADING_SPACE
		#define PLANE_DATA_COMPLEX_SPECTRAL_FOR_IDX(CORE)					\
			SWEET_THREADING_SPACE_PARALLEL_FOR			\
			for (int r = 0; r < 4; r++)								\
			{														\
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2		\
				for (std::size_t jj = planeDataConfig->spectral_complex_ranges[r][1][0]; jj < planeDataConfig->spectral_complex_ranges[r][1][1]; jj++)		\
				{				\
					for (std::size_t ii = planeDataConfig->spectral_complex_ranges[r][0][0]; ii < planeDataConfig->spectral_complex_ranges[r][0][1]; ii++)	\
					{			\
						std::size_t idx = jj*planeDataConfig->spectral_complex_data_size[0]+ii;	\
						CORE	\
					}			\
				}				\
			}



		#define PLANE_DATA_COMPLEX_SPECTRAL_FOR_IDX_ALL(CORE)				\
			SWEET_THREADING_SPACE_PARALLEL_FOR			\
			for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)		\
			{			\
				CORE	\
			}
	#else

		#define PLANE_DATA_COMPLEX_SPECTRAL_FOR_IDX(CORE)					\
			for (int r = 0; r < 4; r++)								\
			{														\
				for (std::size_t jj = planeDataConfig->spectral_complex_ranges[r][1][0]; jj < planeDataConfig->spectral_complex_ranges[r][1][1]; jj++)		\
				{				\
					for (std::size_t ii = planeDataConfig->spectral_complex_ranges[r][0][0]; ii < planeDataConfig->spectral_complex_ranges[r][0][1]; ii++)	\
					{			\
						std::size_t idx = jj*planeDataConfig->spectral_complex_data_size[0]+ii;	\
						CORE	\
					}			\
				}				\
			}


		#define PLANE_DATA_COMPLEX_SPECTRAL_FOR_IDX_ALL(CORE)				\
			for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)		\
			{			\
				CORE	\
			}

	#endif


class PlaneData_SpectralComplex
{
public:
	const PlaneDataConfig *planeDataConfig;

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
	PlaneData_SpectralComplex()	:
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
	}


public:
	PlaneData_SpectralComplex(
			const PlaneDataConfig *i_planeDataConfig
	)	:
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
		assert(i_planeDataConfig != 0);

		setup(i_planeDataConfig);
	}


public:
	PlaneData_SpectralComplex(
			const PlaneData_SpectralComplex &i_plane_data
	)	:
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
		setup(i_plane_data.planeDataConfig);

		operator=(i_plane_data);
	}




public:
	PlaneData_SpectralComplex(
			PlaneData_SpectralComplex &&i_data
	)	:
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
		if (planeDataConfig == nullptr)
			setup(i_data.planeDataConfig);

		assert(i_data.spectral_space_data != nullptr);

		std::swap(spectral_space_data, i_data.spectral_space_data);
	}





	PlaneData_SpectralComplex(
			const PlaneData_PhysicalComplex &i_plane_data
	):
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
		PlaneData_PhysicalComplex tmp = i_plane_data;

		setup(i_plane_data.planeDataConfig);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		planeDataConfig->fft_complex_physical_to_spectral(tmp.physical_space_data, this->spectral_space_data);
	}




	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check_planeDataConfig_identical_res(const PlaneDataConfig *i_planeDataConfig)	const
	{
		assert(planeDataConfig->spectral_complex_data_size[0] == i_planeDataConfig->spectral_complex_data_size[0]);
		assert(planeDataConfig->spectral_complex_data_size[1] == i_planeDataConfig->spectral_complex_data_size[1]);
	}





public:
	PlaneData_SpectralComplex& operator=(
			const PlaneData_SpectralComplex &i_plane_data
	)
	{
		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		assert(i_plane_data.spectral_space_data);
		parmemcpy(spectral_space_data, i_plane_data.spectral_space_data, sizeof(Tcomplex)*planeDataConfig->spectral_complex_array_data_number_of_elements);

		return *this;
	}




public:
	PlaneData_SpectralComplex& operator=(
			PlaneData_SpectralComplex &&i_plane_data
	)
	{
		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		assert(i_plane_data.spectral_space_data);
		std::swap(spectral_space_data, i_plane_data.spectral_space_data);

		return *this;
	}



	PlaneData_SpectralComplex spectral_returnWithTruncatedModes(
			const PlaneDataConfig *i_planeDataConfigTargetTruncation
	)	const
	{
		return spectral_returnWithDifferentModes(i_planeDataConfigTargetTruncation).spectral_returnWithDifferentModes(planeDataConfig);
	}


public:
	PlaneData_SpectralComplex spectral_returnWithDifferentModes(
			const PlaneDataConfig *i_planeDataConfigNew
	)	const
	{
		PlaneData_SpectralComplex out(i_planeDataConfigNew);

		/*
		 *  0 = invalid
		 * -1 = scale down
		 *  1 = scale up
		 */
		int scaling_mode = 0;

		if (planeDataConfig->spectral_complex_data_size[0] < out.planeDataConfig->spectral_complex_data_size[0])
		{
			scaling_mode = 1;
		}
		else if (planeDataConfig->spectral_complex_data_size[0] > out.planeDataConfig->spectral_complex_data_size[0])
		{
			scaling_mode = -1;
		}


		if (planeDataConfig->spectral_complex_data_size[1] < out.planeDataConfig->spectral_complex_data_size[1])
		{
			assert(scaling_mode != -1);
			scaling_mode = 1;
		}
		else if (planeDataConfig->spectral_complex_data_size[1] > out.planeDataConfig->spectral_complex_data_size[1])
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
			for (int n = 0; n <= out.planeDataConfig->spectral_complex_data_size[1]; n++)
			{
				int src_idx = planeDataConfig->getArrayIndexByModes_Complex(n, -n);
				int dst_idx = out.planeDataConfig->getArrayIndexByModes_Complex(n, -n);

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
			for (int n = 0; n <= planeDataConfig->spectral_complex_data_size[1]; n++)
			{
				int src_idx = planeDataConfig->getArrayIndexByModes_Complex(n, -n);
				int dst_idx = out.planeDataConfig->getArrayIndexByModes_Complex(n, -n);

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
	 * Setup spectral plane data based on data in physical space
	 */
	void loadSphereDataPhysical(
			const PlaneData_PhysicalComplex &i_planeDataPhysical
	)
	{
		/**
		 * Warning: The sphat_to_SH function is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 * Hence, we create a copy
		 */
		PlaneData_PhysicalComplex tmp(i_planeDataPhysical);
		planeDataConfig->fft_complex_physical_to_spectral(tmp.physical_space_data, this->spectral_space_data);
	}



	PlaneData_PhysicalComplex toPhys()	const
	{
		return getSphereDataPhysicalComplex();
	}

	PlaneData_PhysicalComplex getSphereDataPhysicalComplex()	const
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		/*
		 * WARNING:
		 * We have to use a temporary array here because of destructive SH transformations
		 */
		PlaneData_SpectralComplex tmp = *this;
		planeDataConfig->fft_complex_spectral_to_physical(tmp->spectral_space_data, out.spectral_space_data);

		return out;
	}


	PlaneData_SpectralComplex operator+(
			const PlaneData_SpectralComplex &i_plane_data
	)	const
	{
		check_planeDataConfig_identical_res(i_plane_data.planeDataConfig);

		PlaneData_SpectralComplex out_plane_data(planeDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_plane_data.spectral_space_data[idx] = spectral_space_data[idx] + i_plane_data.spectral_space_data[idx];

		return out_plane_data;
	}



	PlaneData_SpectralComplex& operator+=(
			const PlaneData_SpectralComplex &i_plane_data
	)
	{
		check_planeDataConfig_identical_res(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] += i_plane_data.spectral_space_data[idx];

		return *this;
	}


	PlaneData_SpectralComplex& operator-=(
			const PlaneData_SpectralComplex &i_plane_data
	)
	{
		check_planeDataConfig_identical_res(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_plane_data.spectral_space_data[idx];

		return *this;
	}



	PlaneData_SpectralComplex operator-(
			const PlaneData_SpectralComplex &i_plane_data
	)	const
	{
		check_planeDataConfig_identical_res(i_plane_data.planeDataConfig);

		PlaneData_SpectralComplex out_plane_data(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_plane_data.spectral_space_data[idx] = spectral_space_data[idx] - i_plane_data.spectral_space_data[idx];

		return out_plane_data;
	}


	PlaneData_SpectralComplex operator-()	const
	{
		PlaneData_SpectralComplex out_plane_data(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_plane_data.spectral_space_data[idx] = -spectral_space_data[idx];

		return out_plane_data;
	}




	PlaneData_SpectralComplex operator*(
			const double i_value
	)	const
	{
		PlaneData_SpectralComplex out_plane_data(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR

		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_plane_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		return out_plane_data;
	}



	const PlaneData_SpectralComplex& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}


	PlaneData_SpectralComplex operator/(
			double i_value
	)	const
	{
		PlaneData_SpectralComplex out_plane_data(planeDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_plane_data.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		return out_plane_data;
	}


	const PlaneData_SpectralComplex& operator/=(
			const std::complex<double> &i_value
	)	const
	{

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			spectral_space_data[idx] /= i_value;

		return *this;
	}



	PlaneData_SpectralComplex operator+(
			double i_value
	)	const
	{
		PlaneData_SpectralComplex out_plane_data(*this);

		out_plane_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		return out_plane_data;
	}



	PlaneData_SpectralComplex& operator+=(
			double i_value
	)
	{
		spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);
		return *this;
	}



	PlaneData_SpectralComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		PlaneData_SpectralComplex out_plane_data(*this);

		out_plane_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

		return out_plane_data;
	}



	PlaneData_SpectralComplex operator-(
			const std::complex<double> &i_value
	)	const
	{
		PlaneData_SpectralComplex out_plane_data(*this);

		out_plane_data.spectral_space_data[0] -= i_value*std::sqrt(4.0*M_PI);

		return out_plane_data;
	}



	PlaneData_SpectralComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		PlaneData_SpectralComplex out_plane_data(planeDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
			out_plane_data.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		return out_plane_data;
	}




public:
	void setup_if_required(
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		if (planeDataConfig != nullptr)
			return;

		setup(i_planeDataConfig);
	}


public:
	void setup(
			const PlaneDataConfig *i_planeConfig
	)
	{
		// assure that the initialization is not done twice!
		assert(planeDataConfig == nullptr);

		planeDataConfig = i_planeConfig;

		spectral_space_data = MemBlockAlloc::alloc<Tcomplex>(planeDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
	}


public:
	~PlaneData_SpectralComplex()
	{
		if (spectral_space_data != nullptr)
		{
			MemBlockAlloc::free(spectral_space_data, planeDataConfig->spectral_complex_array_data_number_of_elements * sizeof(Tcomplex));
		}
	}



	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	inline
	PlaneData_SpectralComplex spectral_solve_helmholtz(
			const std::complex<double> &i_a,
			const std::complex<double> &i_b,
			double r
	)	const
	{
		PlaneData_SpectralComplex out(*this);

		const std::complex<double> a = i_a;
		const std::complex<double> b = i_b;

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				io_data /= (a + (-b*(m * m + n * n));
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
		for (int n = 0; n <= planeDataConfig->spectral_complex_data_size[1]; n++)
		{
			int idx = planeDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				i_lambda(n, m, spectral_space_data[idx]);
				idx++;
			}
		}
	}


	inline
	void spectral_update_lambda_array_indices(
			std::function<void(int,int,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		PLANE_DATA_COMPLEX_SPECTRAL_FOR_IDX_ALL (
				i_lambda(idx, spectral_space_data[idx]);
		);
		spectral_zeroAliasingModes();
	}


	inline
	void spectral_update_lambda_modes(
			std::function<void(int,int,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		//		int modes_0 = planeDataConfig->spectral_complex_data_size[0];
		int modes_1 = planeDataConfig->spectral_complex_data_size[1];

		//		int half_modes_0 = modes_0/2;
		int half_modes_1 = modes_1/2;

		PLANE_DATA_SPECTRAL_FOR_IDX(
				{
			int k0 = ii;
			//					if (k0 > half_modes_0)
			//						k0 -= modes_0;

			int k1 = jj;
			if (k1 > half_modes_1)
				k1 = k1 - modes_1;

			i_lambda(k0, k1, spectral_space_data[idx]);
				}
		);

		spectral_zeroAliasingModes();
	}


	inline
	const std::complex<double>& spectral_get_(
			int in,
			int im
	)	const
	{
		assert(in >= 0);
		assert(in <= planeDataConfig->spectral_complex_data_size[1]);
		assert(std::abs(im) <= planeDataConfig->spectral_complex_data_size[0]);
		assert(std::abs(im) <= in);

		return spectral_space_data[planeDataConfig->getArrayIndexByModes_Complex(in, im)];
	}



	/*
	 * Set all values to zero in spectral space
	 */
	void spectral_set_zero()
	{

	SWEET_THREADING_SPACE_PARALLEL_FOR

		for (int n = 0; n <= planeDataConfig->spectral_complex_data_size[1]; n++)
		{
			int idx = planeDataConfig->getArrayIndexByModes_Complex(n, -n);
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
		for (int m = 0; m <= planeDataConfig->spectral_complex_data_size[0]; m++)
		{
			for (int n = m; n <= planeDataConfig->spectral_complex_data_size[1]; n++)
			{
				std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
				std::cout << spectral_space_data[idx] << "\t";
			}
			std::cout << std::endl;
		}
	}

};



inline
static
PlaneData_SpectralComplex operator*(
		const double i_value,
		const PlaneData_SpectralComplex &i_array_data
)
{
	return i_array_data*i_value;
}


inline
static
PlaneData_SpectralComplex operator*(
		const std::complex<double> &i_value,
		const PlaneData_SpectralComplex &i_array_data
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
PlaneData_SpectralComplex operator+(
		const double i_value,
		const PlaneData_SpectralComplex &i_array_data
)
{
	return ((PlaneData_SpectralComplex&)i_array_data)+i_value;
}

inline
static
PlaneData_SpectralComplex operator+(
		const std::complex<double> &i_value,
		const PlaneData_SpectralComplex &i_array_data
)
{
	return i_array_data+i_value;
}

inline
static
PlaneData_SpectralComplex operator-(
		const std::complex<double> &i_value,
		const PlaneData_SpectralComplex &i_array_data
)
{
	PlaneData_SpectralComplex out_plane_data(i_array_data.planeDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (std::size_t idx = 0; idx < i_array_data.planeDataConfig->spectral_complex_array_data_number_of_elements; idx++)
		out_plane_data.spectral_space_data[idx] = -i_array_data.spectral_space_data[idx];

	out_plane_data.spectral_space_data[0] += i_value*std::sqrt(4.0*M_PI);

	return out_plane_data;

}

#endif
