/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *  * 14 Mar 2022, Joao Steinstraesser <joao.steinstraesser@usp.br>
 *    Split into physical and spectral classes
 */

#ifndef SWEET_PLANE_DATA_SPECTRAL_HPP_
#define SWEET_PLANE_DATA_SPECTRAL_HPP_

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
#include <sweet/core/plane/PlaneDataConfig.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/core/SWEETError.hpp>


#define PLANE_DATA_SPECTRAL_FOR_IDX(CORE)					\
		SWEET_THREADING_SPACE_PARALLEL_FOR			\
		for (int r = 0; r < 2; r++)								\
		{														\
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2		\
			for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
			{				\
				for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
				{			\
					std::size_t idx = jj*planeDataConfig->spectral_data_size[0]+ii;	\
					CORE	\
				}			\
			}				\
		}

#if SWEET_USE_LIBFFT
#	include <fftw3.h>
#endif

#if SWEET_THREADING_SPACE
#	include <omp.h>
#endif

namespace sweet
{

class PlaneData_Spectral
{
	typedef std::complex<double> Tcomplex;

public:
	const PlaneDataConfig *planeDataConfig = nullptr;

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
			PlaneData_Spectral &i_planeData
	)
	{
		assert(planeDataConfig == i_planeData.planeDataConfig);

		std::swap(spectral_space_data, i_planeData.spectral_space_data);
	}


public:
	PlaneData_Spectral(
			const PlaneDataConfig *i_planeDataConfig
	)	:
		spectral_space_data(nullptr)
	{
		assert(i_planeDataConfig != 0);

		setup(i_planeDataConfig);
	}


public:
	PlaneData_Spectral(
			const PlaneDataConfig *i_planeDataConfig,
			const std::complex<double> &i_value
	)	:
		planeDataConfig(i_planeDataConfig),
		spectral_space_data(nullptr)
	{
		assert(i_planeDataConfig != 0);

		setup(i_planeDataConfig);
		spectral_set_value(i_value);
	}



public:
	PlaneData_Spectral(
			const PlaneDataConfig *i_planeDataConfig,
			double &i_value
	)	:
		planeDataConfig(i_planeDataConfig),
		spectral_space_data(nullptr)
	{
		assert(i_planeDataConfig != 0);

		setup(i_planeDataConfig);
		spectral_set_value(i_value);
	}



	/**
	 * Without setup where we need to call setup(...) later on
	 */
public:
	PlaneData_Spectral()	:
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
	}


	/**
	 * dummy initialization by handing over an unused integer
	 */
public:
	PlaneData_Spectral(int i)	:
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
{
}


public:
	PlaneData_Spectral(
			const PlaneData_Spectral &i_plane_data
	)	:
		planeDataConfig(i_plane_data.planeDataConfig),
		spectral_space_data(nullptr)
	{
		if (i_plane_data.planeDataConfig == nullptr)
			return;

		alloc_data();

		operator=(i_plane_data);
	}



public:
	PlaneData_Spectral(
			PlaneData_Spectral &&i_plane_data
	)	:
		planeDataConfig(i_plane_data.planeDataConfig),
		spectral_space_data(nullptr)
	{
		if (i_plane_data.planeDataConfig == nullptr)
			return;

		std::swap(spectral_space_data, i_plane_data.spectral_space_data);
	}




	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void _validateRes(
			const PlaneDataConfig *i_planeDataConfig
	)	const
	{
		assert(planeDataConfig->physical_res[0] == i_planeDataConfig->physical_res[0]);
		assert(planeDataConfig->physical_res[1] == i_planeDataConfig->physical_res[1]);

		assert(planeDataConfig->spectral_data_size[0] == i_planeDataConfig->spectral_data_size[0]);
		assert(planeDataConfig->spectral_data_size[1] == i_planeDataConfig->spectral_data_size[1]);
	}


public:
	PlaneData_Spectral& operator=(
			const PlaneData_Spectral &i_plane_data
	)
	{
		if (i_plane_data.planeDataConfig == nullptr)
			return *this;

		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		parmemcpy(spectral_space_data, i_plane_data.spectral_space_data, sizeof(Tcomplex)*planeDataConfig->spectral_array_data_number_of_elements);

		return *this;
	}

public:
	/**
	 * assignment operator
	 */
	PlaneData_Spectral &operator=(double i_value)
	{
		spectral_set_value(i_value);

		return *this;
	}

public:
	/**
	 * assignment operator
	 */
	PlaneData_Spectral &operator=(int i_value)
	{
		spectral_set_value(i_value);

		return *this;
	}


public:
	void spectral_zeroAliasingModes()
	{
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING || 1	/// ALWAYS run this to eliminate Nyquist Frequency even without dealiasing activated

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int k = 0; k < 2; k++)
		{
			if (k == 0)
			{
				/*
				 * First process part between top and bottom spectral data blocks
				 */
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[0][1][1]; jj < planeDataConfig->spectral_data_iteration_ranges[1][1][0]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii++)
					{
						//spectral_space_data[jj*planeDataConfig->spectral_data_size[0]+ii] = 0;
						spectral_space_data[planeDataConfig->getArrayIndexByModes(jj, ii)] = 0;
					}
			}
			else
			{
				/*
				 * Then process the aliasing block on the right side
				 */
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (std::size_t jj = 0; jj < planeDataConfig->spectral_data_size[1]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii < planeDataConfig->spectral_data_size[0]; ii++)
					{
						//spectral_space_data[jj*planeDataConfig->spectral_data_size[0]+ii] = 0;
						spectral_space_data[planeDataConfig->getArrayIndexByModes(jj, ii)] = 0;
					}
			}
		}
#else

#endif
	}

public:
	void spectral_debugCheckForZeroAliasingModes()	const
	{
#if SWEET_DEBUG

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int k = 0; k < 2; k++)
		{
			if (k == 0)
			{
				/*
				 * First process part between top and bottom spectral data blocks
				 */
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[0][1][1]; jj < planeDataConfig->spectral_data_iteration_ranges[1][1][0]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii++)
					{
						std::complex<double> &data = spectral_space_data[planeDataConfig->getArrayIndexByModes(jj, ii)];

						double error = std::sqrt(data.real()*data.real() + data.imag()*data.imag());
						if (error >= 1e-9)
						{
							print_spectralData_zeroNumZero();
							std::cout << "Value at spectral coordinate " << jj << ", " << ii << " should be zero, but is " << data << std::endl;
							SWEETError("EXIT");
						}
					}
			}
			else
			{
				/*
				 * Then process the aliasing block on the right side
				 */
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (std::size_t jj = 0; jj < planeDataConfig->spectral_data_size[1]; jj++)
					for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[0][0][1]; ii < planeDataConfig->spectral_data_size[0]; ii++)
					{
						std::complex<double> &data = spectral_space_data[planeDataConfig->getArrayIndexByModes(jj, ii)];

						double error = std::sqrt(data.real()*data.real() + data.imag()*data.imag());
						if (error >= 1e-9)
						{
							print_spectralData_zeroNumZero();
							std::cout << "Value at spectral coordinate " << jj << ", " << ii << " should be zero, but is " << data << std::endl;
							SWEETError("EXIT");
						}
					}
			}
		}
#endif
	}




	/**
	 * This function implements copying the spectral data only.
	 *
	 * This becomes handy if coping with data which should be only transformed without dealiasing.
	 */
public:
	PlaneData_Spectral& load_nodealiasing(
			const PlaneData_Spectral &i_plane_data		///< data to be converted to planeDataConfig_nodealiasing
	)
	{
		if (planeDataConfig == nullptr)
			SWEETError("planeDataConfig not initialized");

		parmemcpy(spectral_space_data, i_plane_data.spectral_space_data, sizeof(Tcomplex)*planeDataConfig->spectral_array_data_number_of_elements);

		return *this;
	}


public:
	PlaneData_Spectral& operator=(
			PlaneData_Spectral &&i_plane_data
	)
	{
		if (planeDataConfig == nullptr)
			setup(i_plane_data.planeDataConfig);

		std::swap(spectral_space_data, i_plane_data.spectral_space_data);

		return *this;
	}


public:
	PlaneData_Spectral spectral_returnWithDifferentModes(
			const PlaneDataConfig *i_planeDataConfig
	)	const
	{
		PlaneData_Spectral out(i_planeDataConfig);

		/*
		 *  0 = invalid
		 * -1 = scale down
		 *  1 = scale up
		 */
		int scaling_mode = 0;

		if (planeDataConfig->spectral_modes[0] < out.planeDataConfig->spectral_modes[0])
		{
			scaling_mode = 1;
		}
		else if (planeDataConfig->spectral_modes[0] > out.planeDataConfig->spectral_modes[0])
		{
			scaling_mode = -1;
		}

		if (planeDataConfig->spectral_modes[1] < out.planeDataConfig->spectral_modes[1])
		{
			assert(scaling_mode != -1);
			scaling_mode = 1;
		}
		else if (planeDataConfig->spectral_modes[1] > out.planeDataConfig->spectral_modes[1])
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

		double rescale =
				(double)(out.planeDataConfig->physical_array_data_number_of_elements)
				/
				(double)(planeDataConfig->physical_array_data_number_of_elements);

		//rescale = 1.0;

		{
			if (scaling_mode == -1)
			{
				/*
				 * more modes -> less modes
				 */

				/*
				 * Region #1
				 *
				 * 00000000 7
				 * 00000000 6
				 * 00000000 5
				 * 00000000 4
				 * 00000000 3
				 * XXXXX000 2
				 * XXXXX000 1
				 * XXXXX000 0
				 */
				{
					const std::size_t* src_range_dim0 = &(planeDataConfig->spectral_data_iteration_ranges[0][0][0]);
					const std::size_t* src_range_dim1 = &(planeDataConfig->spectral_data_iteration_ranges[0][1][0]);

					const std::size_t* dst_range_dim0 = &(out.planeDataConfig->spectral_data_iteration_ranges[0][0][0]);
					const std::size_t* dst_range_dim1 = &(out.planeDataConfig->spectral_data_iteration_ranges[0][1][0]);

					assert(src_range_dim0[0] == 0);
					assert(dst_range_dim0[0] == 0);
					assert(src_range_dim1[0] == 0);
					assert(dst_range_dim1[0] == 0);

					std::size_t dst_size = dst_range_dim0[1];//-dst_range_dim0[0];

					SWEET_THREADING_SPACE_PARALLEL_FOR
					for (std::size_t j = 0; j < dst_range_dim1[1]; j++)
					{
						std::complex<double> *src = &spectral_space_data[planeDataConfig->spectral_data_size[0]*j];
						std::complex<double> *dst = &out.spectral_space_data[out.planeDataConfig->spectral_data_size[0]*j];

						for (std::size_t i = 0; i < dst_size; i++)
							dst[i] = src[i]*rescale;
					}
				}


				/*
				 * Region #2
				 *
				 * XXXXX000 7
				 * XXXXX000 6
				 * XXXXX000 5
				 * 00000000 4
				 * 00000000 3
				 * 00000000 2
				 * 00000000 1
				 * 00000000 0
				 */
				{
					const std::size_t* src_range_dim0 = &(planeDataConfig->spectral_data_iteration_ranges[1][0][0]);
					const std::size_t* src_range_dim1 = &(planeDataConfig->spectral_data_iteration_ranges[1][1][0]);

					const std::size_t* dst_range_dim0 = &(out.planeDataConfig->spectral_data_iteration_ranges[1][0][0]);
					const std::size_t* dst_range_dim1 = &(out.planeDataConfig->spectral_data_iteration_ranges[1][1][0]);

					assert(src_range_dim0[0] == 0);
					assert(dst_range_dim0[0] == 0);

					std::size_t dst_size = dst_range_dim0[1];//-dst_range_dim0[0];

					SWEET_THREADING_SPACE_PARALLEL_FOR
					for (std::size_t j = dst_range_dim1[0]; j < dst_range_dim1[1]; j++)
					{
						std::complex<double> *src = &spectral_space_data[planeDataConfig->spectral_data_size[0]*(src_range_dim1[1]-(dst_range_dim1[1]-j))];
						std::complex<double> *dst = &out.spectral_space_data[out.planeDataConfig->spectral_data_size[0]*j];

						for (std::size_t i = 0; i < dst_size; i++)
							dst[i] = src[i]*rescale;
					}
				}

				out.spectral_zeroAliasingModes();
			}
			else
			{
				/*
				 * less modes -> more modes
				 */

				/*
				 * Region #1
				 *
				 * 00000000 7
				 * 00000000 6
				 * 00000000 5
				 * 00000000 4
				 * 00000000 3
				 * XXXXX000 2
				 * XXXXX000 1
				 * XXXXX000 0
				 */
				out.spectral_set_zero();

				{
					int r = 0;

					const std::size_t* src_range_dim0 = &(planeDataConfig->spectral_data_iteration_ranges[r][0][0]);
					const std::size_t* src_range_dim1 = &(planeDataConfig->spectral_data_iteration_ranges[r][1][0]);

					const std::size_t* dst_range_dim0 = &(out.planeDataConfig->spectral_data_iteration_ranges[r][0][0]);
					const std::size_t* dst_range_dim1 = &(out.planeDataConfig->spectral_data_iteration_ranges[r][1][0]);

					assert(src_range_dim0[0] == 0);
					assert(dst_range_dim0[0] == 0);
					assert(src_range_dim1[0] == 0);

					std::size_t src_size = src_range_dim0[1];//-dst_range_dim0[0];


					SWEET_THREADING_SPACE_PARALLEL_FOR
					for (std::size_t j = 0; j < src_range_dim1[1]; j++)
					{
						std::complex<double> *src = &spectral_space_data[planeDataConfig->spectral_data_size[0]*(j-src_range_dim1[0]+dst_range_dim1[0])];
						std::complex<double> *dst = &out.spectral_space_data[out.planeDataConfig->spectral_data_size[0]*j];

						for (std::size_t i = 0; i < src_size; i++)
							dst[i] = src[i]*rescale;
					}
				}


				/*
				 * Region #2
				 *
				 * XXXXX000 7
				 * XXXXX000 6
				 * XXXXX000 5
				 * 00000000 4
				 * 00000000 3
				 * 00000000 2
				 * 00000000 1
				 * 00000000 0
				 */
				{
					int r = 1;

					const std::size_t* src_range_dim0 = &(planeDataConfig->spectral_data_iteration_ranges[r][0][0]);
					const std::size_t* src_range_dim1 = &(planeDataConfig->spectral_data_iteration_ranges[r][1][0]);

					const std::size_t* dst_range_dim0 = &(out.planeDataConfig->spectral_data_iteration_ranges[r][0][0]);
					const std::size_t* dst_range_dim1 = &(out.planeDataConfig->spectral_data_iteration_ranges[r][1][0]);

					assert(src_range_dim0[0] == 0);
					assert(dst_range_dim0[0] == 0);

					std::size_t src_size0 = src_range_dim0[1];//-dst_range_dim0[0];
					std::size_t src_size1 = src_range_dim1[1]-src_range_dim1[0];


					SWEET_THREADING_SPACE_PARALLEL_FOR
					for (std::size_t j = src_range_dim1[0]; j < src_range_dim1[1]; j++)
					{
						std::complex<double> *src = &spectral_space_data[planeDataConfig->spectral_data_size[0]*j];
						std::complex<double> *dst = &out.spectral_space_data[out.planeDataConfig->spectral_data_size[0]*(dst_range_dim1[1]-src_size1+(j-src_range_dim1[0]))];

						for (std::size_t i = 0; i < src_size0; i++)
						{
#if SWEET_DEBUG
							assert((int)(out.spectral_space_data - &dst[i]) < (int)out.planeDataConfig->spectral_array_data_number_of_elements);
#endif
							dst[i] = src[i]*rescale;
						}
					}

				}
			}

		}

		return out;
	}


	/**
	 * Return Plane Array with all spectral coefficients a+bi --> 1/(a+bi)
	 */
	inline
	PlaneData_Spectral spectral_invert()	const
	{
		PlaneData_Spectral out(planeDataConfig);

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = 1.0/spectral_space_data[idx];
		);

		out.spectral_zeroAliasingModes();

		return out;
	}



	/*
	 * Setup spectral sphere data based on data in physical space
	 */
	void loadPlaneDataPhysical(
		const PlaneData_Physical &i_planeDataPhysical
	)
	{
		/**
		 * Warning: The fftw functions are in-situ operations.
		 * Therefore, the data in the source array will be destroyed.
		 * Hence, we create a copy
		 */
		PlaneData_Physical tmp(i_planeDataPhysical);
		planeDataConfig->fft_physical_to_spectral(tmp.physical_space_data, this->spectral_space_data);
		// ALWAYS zero aliasing modes after doing transformation to spectral space
		this->spectral_zeroAliasingModes();
	}



	/*
	 * Return the data converted to physical space
	 */
	PlaneData_Physical getPlaneDataPhysical()	const
	{
		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		PlaneData_Spectral tmp(*this);
		PlaneData_Physical retval(planeDataConfig);

		planeDataConfig->fft_spectral_to_physical(tmp.spectral_space_data, retval.physical_space_data);
		return retval;
	}


	/*
	 * Return the data converted to physical space
	 *
	 * alias for "getPlaneDataPhysical"
	 */
	PlaneData_Physical toPhys()	const
	{
		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		PlaneData_Spectral tmp(*this);
		PlaneData_Physical retval(planeDataConfig);
		planeDataConfig->fft_spectral_to_physical(tmp.spectral_space_data, retval.physical_space_data);

		return retval;
	}


	PlaneData_PhysicalComplex getPlaneDataPhysicalComplex()	const
	{
		PlaneData_PhysicalComplex out(planeDataConfig);

		/*
		 * WARNING:
		 * We have to use a temporary array here because of destructive FFTW transformations
		 */
		PlaneData_Spectral tmp_spectral(*this);
		PlaneData_Physical tmp_physical(planeDataConfig);
		planeDataConfig->fft_spectral_to_physical(tmp_spectral.spectral_space_data, tmp_physical.physical_space_data);

		parmemcpy(out.physical_space_data, tmp_physical.physical_space_data, sizeof(double)*planeDataConfig->physical_array_data_number_of_elements);

		return out;
	}



	PlaneData_Spectral(
			const PlaneData_Physical &i_plane_data_physical
	)	:
		planeDataConfig(nullptr),
		spectral_space_data(nullptr)
	{
		setup(i_plane_data_physical.planeDataConfig);

		loadPlaneDataPhysical(i_plane_data_physical);
	}



private:
	void _main_setup(
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(planeDataConfig == nullptr);

		planeDataConfig = i_planeDataConfig;
		alloc_data();
	}



public:
	void setup(
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		_main_setup(i_planeDataConfig);
	}



	/*
	 * Wrapper for main setup
	 */
public:
	void setup(
		const PlaneDataConfig &i_planeDataConfig
	)
	{
		_main_setup(&i_planeDataConfig);
	}

public:
	void setup(
		const PlaneDataConfig *i_planeDataConfig,
		double i_value
	)
	{
		_main_setup(i_planeDataConfig);

		spectral_set_value(i_value);
	}


private:
	void alloc_data()
	{
		assert(spectral_space_data == nullptr);
		spectral_space_data = MemBlockAlloc::alloc<Tcomplex>(planeDataConfig->spectral_array_data_number_of_elements * sizeof(Tcomplex));
	}


public:
	void setup_if_required(
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		if (planeDataConfig != nullptr)
			return;

		_main_setup(i_planeDataConfig);
	}


public:
	void clear()
	{
		if (spectral_space_data != nullptr)
		{
			MemBlockAlloc::free(spectral_space_data, planeDataConfig->spectral_array_data_number_of_elements * sizeof(Tcomplex));

			spectral_space_data = nullptr;
			planeDataConfig = nullptr;
		}
	}


public:
	~PlaneData_Spectral()
	{
		clear();
	}



	PlaneData_Spectral operator+(
			const PlaneData_Spectral &i_plane_data
	)	const
	{
		_validateRes(i_plane_data.planeDataConfig);

		PlaneData_Spectral out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx] + i_plane_data.spectral_space_data[idx];

		out.spectral_zeroAliasingModes();

		return out;
	}



	PlaneData_Spectral& operator+=(
			const PlaneData_Spectral &i_plane_data
	)
	{
		_validateRes(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] += i_plane_data.spectral_space_data[idx];

		this->spectral_zeroAliasingModes();

		return *this;
	}


	PlaneData_Spectral& operator-=(
			const PlaneData_Spectral &i_plane_data
	)
	{
		_validateRes(i_plane_data.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] -= i_plane_data.spectral_space_data[idx];

		this->spectral_zeroAliasingModes();

		return *this;
	}



	PlaneData_Spectral operator-(
			const PlaneData_Spectral &i_plane_data
	)	const
	{
		_validateRes(i_plane_data.planeDataConfig);

		PlaneData_Spectral out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx] - i_plane_data.spectral_space_data[idx];

		out.spectral_zeroAliasingModes();

		return out;
	}


	PlaneData_Spectral operator-(
			const PlaneData_Physical &i_plane_data_physical
	)	const
	{
		_validateRes(i_plane_data_physical.planeDataConfig);

		PlaneData_Spectral i_plane_data_spectral(planeDataConfig);
		PlaneData_Spectral out(planeDataConfig);

		i_plane_data_spectral.loadPlaneDataPhysical(i_plane_data_physical);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx] - i_plane_data_spectral.spectral_space_data[idx];

		out.spectral_zeroAliasingModes();

		return out;
	}



	PlaneData_Spectral operator-()	const
	{
		PlaneData_Spectral out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = -spectral_space_data[idx];

		out.spectral_zeroAliasingModes();

		return out;
	}



	PlaneData_Spectral operator*(
			const PlaneData_Spectral &i_plane_data
	)	const
	{
		_validateRes(i_plane_data.planeDataConfig);

		PlaneData_Physical a = this->toPhys();
		PlaneData_Physical b = i_plane_data.toPhys();

		PlaneData_Physical mul(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			mul.physical_space_data[i] = a.physical_space_data[i]*b.physical_space_data[i];
		}

		// Dealiasing is performed inside the following call
		PlaneData_Spectral out(mul);

		return out;
	}


	PlaneData_Physical multiplication_physical_space(
				const PlaneData_Physical &i_a,
				const PlaneData_Physical &i_b
	) const
	{
		_validateRes(i_a.planeDataConfig);
		_validateRes(i_b.planeDataConfig);

		PlaneData_Physical mul(planeDataConfig);
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			mul.physical_space_data[i] = i_a.physical_space_data[i]*i_b.physical_space_data[i];
		}

		PlaneData_Spectral out_spec(mul);
		return out_spec.toPhys();
	}


	PlaneData_Spectral operator/(
			const PlaneData_Spectral &i_plane_data
	)	const
	{
		_validateRes(i_plane_data.planeDataConfig);

		PlaneData_Physical a = this->toPhys();
		PlaneData_Physical b = i_plane_data.toPhys();

		PlaneData_Physical div(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			div.physical_space_data[i] = a.physical_space_data[i]/b.physical_space_data[i];

		PlaneData_Spectral out(div);

		return out;
	}


	PlaneData_Spectral operator*(
			double i_value
	)	const
	{
		PlaneData_Spectral out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx]*i_value;

		out.spectral_zeroAliasingModes();

		return out;
	}


	const PlaneData_Spectral& operator*=(
			double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}


	const PlaneData_Spectral& operator/=(
			double i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] /= i_value;

		return *this;
	}


	const PlaneData_Spectral& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			spectral_space_data[idx] *= i_value;

		return *this;
	}


	PlaneData_Spectral operator/(
			double i_value
	)	const
	{
		PlaneData_Spectral out(planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = spectral_space_data[idx]/i_value;

		out.spectral_zeroAliasingModes();

		return out;
	}


	PlaneData_Spectral operator+(
			double i_value
	)	const
	{
		PlaneData_Spectral out(*this);

		out.spectral_space_data[0] += i_value * (double)planeDataConfig->physical_array_data_number_of_elements;
		out.spectral_zeroAliasingModes();

		return out;
	}


	PlaneData_Spectral operator_scalar_sub_this(
			double i_value
	)	const
	{
		PlaneData_Spectral out(planeDataConfig);
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
			out.spectral_space_data[idx] = -spectral_space_data[idx];

		out.spectral_space_data[0] = i_value*std::sqrt(4.0*M_PI) + out.spectral_space_data[0];
		return out;
	}


	const PlaneData_Spectral& operator+=(
			double i_value
	)	const
	{
		spectral_space_data[0] += i_value * (double)planeDataConfig->physical_array_data_number_of_elements;

		return *this;
	}


	PlaneData_Spectral operator-(
			double i_value
	)	const
	{
		PlaneData_Spectral out(*this);
		out.spectral_space_data[0] -= i_value * (double)planeDataConfig->physical_array_data_number_of_elements;

		out.spectral_zeroAliasingModes();

		return out;
	}


	PlaneData_Spectral& operator-=(
			double i_value
	)
	{
		//SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		//for (std::size_t idx = 0; idx < planeDataConfig->spectral_array_data_number_of_elements; idx++)
		//	spectral_space_data[idx] -= i_value * (double)planeDataConfig->physical_array_data_number_of_elements;
		spectral_space_data[0] -= i_value * (double)planeDataConfig->physical_array_data_number_of_elements;

		this->spectral_zeroAliasingModes();

		return *this;
	}

public:
	/**
	 * Add scalar to all spectral modes
	 */
	inline
	PlaneData_Spectral spectral_addScalarAll(
			const double &i_value
	)	const
	{
		PlaneData_Spectral out(planeDataConfig);

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = spectral_space_data[idx] + i_value;
		);

		out.spectral_zeroAliasingModes();

		return out;
	}


	/**
	 * Invert the application of a linear operator in spectral space.
	 * The operator is given in i_array_data
	 *
	 * Solving for phi in
	 *
	 * 		\Delta^2 phi = delta_phi
	 *
	 * can be accomplished by
	 *
	 * 		phi = delta_phi.spectral_div_element_wise(\Delta^2)
	 */
	inline
	PlaneData_Spectral spectral_div_element_wise(
			const PlaneData_Spectral &i_array_data	///< operator
	)	const
	{
		PlaneData_Spectral out(planeDataConfig);

		PLANE_DATA_SPECTRAL_FOR_IDX(

				double ar = spectral_space_data[idx].real();
				double ai = spectral_space_data[idx].imag();
				double br = i_array_data.spectral_space_data[idx].real();
				double bi = i_array_data.spectral_space_data[idx].imag();

				double den = br*br+bi*bi;

				if (den == 0)
				{
					out.spectral_space_data[idx].real(0);
					out.spectral_space_data[idx].imag(0);
				}
				else
				{
					out.spectral_space_data[idx].real((ar*br + ai*bi)/den);
					out.spectral_space_data[idx].imag((ai*br - ar*bi)/den);
				}
		);

		out.spectral_zeroAliasingModes();

		return out;
	}



	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	inline
	PlaneData_Spectral spectral_solve_helmholtz(
			const double &i_a,
			const double &i_b,
			double r
	)
	{
		PlaneData_Spectral out(*this);

		const double a = i_a;
		const double b = i_b;

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				io_data /= (a + (-b* (m * m + n * n) )); //???
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
	PlaneData_Spectral spectral_solve_laplace(
			double r
	)
	{
		PlaneData_Spectral out(*this);

		const double b = 1.0;

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				if (n == 0)
					io_data = 0;
				else
					io_data /= (-b*(m * m + n * n));  // ????
			}
		);

		return out;
	}


	/**
	 * Truncate modes which are not representable in spectral space
	 */
	const PlaneData_Spectral& spectral_truncate()	const
	{
		PlaneData_Physical tmp(planeDataConfig);

		planeDataConfig->fft_spectral_to_physical(this->spectral_space_data, tmp.physical_space_data);
		planeDataConfig->fft_physical_to_spectral(tmp.physical_space_data, this->spectral_space_data);

		return *this;
	}


	void spectral_update_lambda(
			std::function<void(int,int,Tcomplex&)> i_lambda
	)
	{
		PLANE_DATA_SPECTRAL_FOR_IDX(
					i_lambda(jj, ii, spectral_space_data[idx]);
				)
	}


	const std::complex<double>& spectral_get(
			int i_n,
			int i_m
	)	const
	{
		assert(i_n >= 0 && i_m >= 0);
		assert(i_n <= (int)planeDataConfig->spectral_data_size[1]);
		assert(i_m <= (int)planeDataConfig->spectral_data_size[0]);
		///assert(i_m <= i_n);

		return spectral_space_data[planeDataConfig->getArrayIndexByModes(i_n, i_m)];
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

		if (i_m > (int)planeDataConfig->spectral_data_size[0])
			SWEETError("Out of boundary b");

		if (i_n > (int)planeDataConfig->spectral_data_size[1])
			SWEETError("Out of boundary c");
#endif

		spectral_space_data[planeDataConfig->getArrayIndexByModes(i_n, i_m)] = i_data;
	}


	/*
	 * Set all values to zero
	 */
	void spectral_set_zero()
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
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
		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
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

		for (std::size_t n = 0; n <= planeDataConfig->spectral_data_size[1]; n++)
			for (std::size_t m = 0; m <= planeDataConfig->spectral_data_size[0]; m++)
			{
				std::size_t idx = planeDataConfig->getArrayIndexByModes(n, m);
				sum += spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
			}

		//////m=0 case - weight 1
		////std::size_t idx = planeDataConfig->getArrayIndexByModes(0, 0);
		////for (std::size_t n = 0; n <= planeDataConfig->spectral_data_size[1]; n++)
		////{
		////	sum += spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
		////	idx++;
		////}
		////
		//////m>0 case - weight 2, as they appear twice
		////for (std::size_t m = 0; m <= planeDataConfig->spectral_data_size[0]; m++)
		////{
		////	std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
		////	for (std::size_t n = m; n <= planeDataConfig->spectral_data_size[1]; n++)
		////	{
		////		sum += 2.0*spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
		////		idx++;
		////	}
		////}

		return sum.real();
	}

	/**
	 * reduce to sum square in spectrum
	 */
	double spectral_reduce_sum_sq()
	{
		double rms = 0.0;

		std::complex<double> sum = 0;

		for (int r = 0; r < 2; r++)								
		{														
			for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
			{	
				for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
				{	
					std::size_t idx = jj*planeDataConfig->spectral_data_size[0]+ii;	
					sum += (spectral_space_data[idx]*std::conj(spectral_space_data[idx]));
				}
			}
		}

		//sum = std::__complex_sqrt (sum/(double)(planeDataConfig->spectral_array_data_number_of_elements));

		if(sum.imag()>DBL_EPSILON)
			SWEETError("Reduce operation of complex values (rms) error");

		rms = sum.real()/(double)(planeDataConfig->spectral_array_data_number_of_elements); 
		//rms = (sum.real()*sum.real()+sum.imag()*sum.imag()); 
		return rms;

	}


	/**
	 * Return the minimum value
	 */
	std::complex<double> spectral_reduce_min()	const
	{
		std::complex<double> error = std::numeric_limits<double>::infinity();

		for (std::size_t j = 0; j < planeDataConfig->spectral_array_data_number_of_elements; j++)
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

		for (std::size_t j = 0; j < planeDataConfig->spectral_array_data_number_of_elements; j++)
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
		for (std::size_t j = 0; j < planeDataConfig->spectral_array_data_number_of_elements; j++)
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

		assert (rnorm <= planeDataConfig->spectral_data_size[0]);
		assert (rnorm <= planeDataConfig->spectral_data_size[1]);

		double error = -std::numeric_limits<double>::infinity();
		std::complex<double> w = {0,0};

		for (std::size_t n = 0; n <= rnorm; n++)
			for (std::size_t m = 0; m <= rnorm; m++)
			{
				std::size_t idx = planeDataConfig->getArrayIndexByModes(n, m);
				w = spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
				error = std::max(std::abs(w), error);
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
		for (std::size_t j = 0; j < planeDataConfig->spectral_array_data_number_of_elements; j++)
		{
			w = spectral_space_data[j]*std::conj(spectral_space_data[j]);
			error = std::min(std::abs(w), error);
		}

		return error;
	}

	/**
	 * reduce to root mean square in spectrum
	 */
	double spectral_reduce_rms()
	{
		double rms = 0.0;

		std::complex<double> sum = 0;

		for (int r = 0; r < 2; r++)								
		{														
			for (std::size_t jj = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; jj < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; jj++)		\
			{	
				for (std::size_t ii = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; ii < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; ii++)	\
				{	
					std::size_t idx = jj*planeDataConfig->spectral_data_size[0]+ii;	
					sum += spectral_space_data[idx]*std::conj(spectral_space_data[idx]);
				}
			}
		}

		//sum = std::__complex_sqrt (sum/(double)(planeDataConfig->spectral_array_data_number_of_elements));

		if(sum.imag()>DBL_EPSILON)
			SWEETError("Reduce operation of complex values (rms) error");

		rms = std::sqrt(sum.real()/(double)(planeDataConfig->spectral_array_data_number_of_elements)); 
		return rms;

	}


	bool spectral_reduce_is_any_nan_or_inf()	const
	{
		bool retval = false;

#if SWEET_THREADING_SPACE
		#pragma omp parallel for simd reduction(|:retval)
#endif
		for (std::size_t j = 0; j < planeDataConfig->spectral_array_data_number_of_elements; j++)
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
		for (std::size_t m = 0; m < planeDataConfig->spectral_data_size[0]; m++)
		{
			///std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
			for (std::size_t n = 0; n < planeDataConfig->spectral_data_size[1]; n++)
			{
				std::size_t idx = planeDataConfig->getArrayIndexByModes(n, m);
				if (std::abs(spectral_space_data[idx]) < i_abs_threshold)
					std::cout << 0 << "\t";
				else
					std::cout << spectral_space_data[idx] << "\t";
				idx++;
			}
			std::cout << std::endl;
		}
	}

	/**
	* Call file_physical_loadData from PlaneData_Physical
	 */
	bool file_physical_loadData(
			const char *i_filename,		///< Name of file to load data from
			bool i_binary_data = false	///< load as binary data (disabled per default)
	)
	{
		PlaneData_Physical phys(this->planeDataConfig);
		bool status = phys.file_physical_loadData(i_filename, i_binary_data);
		if (status)
			this->loadPlaneDataPhysical(phys);
		return status;
	}

	inline
	void print_spectralData()	const
	{
		PlaneData_Spectral &rw_array_data = (PlaneData_Spectral&)*this;

		//for (std::size_t y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--) // https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
		for (std::size_t y = planeDataConfig->spectral_data_size[1]-1; y < planeDataConfig->spectral_data_size[1]; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = rw_array_data.spectral_get(y, x);
				std::cout << "(" << value.real() << ", " << value.imag() << ")\t";
			}
			std::cout << std::endl;
		}
	}


	/**
	 * print spectral data and zero out values which are numerically close to zero
	 */
	inline
	void print_spectralData_zeroNumZero(double i_zero_threshold = 1e-13)	const
	{
		PlaneData_Spectral &rw_array_data = (PlaneData_Spectral&)*this;

		//for (std::size_t y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--) //  https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
		for (std::size_t y = planeDataConfig->spectral_data_size[1]-1; y < planeDataConfig->spectral_data_size[1]; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = rw_array_data.spectral_get(y, x);

				double re = value.real();
				double im = value.imag();

				if (std::abs(re) < i_zero_threshold)	re = 0.0;
				if (std::abs(im) < i_zero_threshold)	im = 0.0;

				std::cout << "(" << re << ", " << im << ")\t";
			}
			std::cout << std::endl;
		}
	}

	inline
	void print_spectralIndex()	const
	{
		PlaneData_Spectral &rw_array_data = (PlaneData_Spectral&)*this;

		//for (std::size_t y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--) //  https://stackoverflow.com/questions/3623263/reverse-iteration-with-an-unsigned-loop-variable
		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = rw_array_data.spectral_get(y, x);
				std::cout << "(" << x << ", "<< y << ", "<< value.real() << ", " << value.imag() << ")\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

	}

	inline
	void print_spectralNonZero()	const
	{
		PlaneData_Spectral &rw_array_data = (PlaneData_Spectral&)*this;

		for (int y = planeDataConfig->spectral_data_size[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = rw_array_data.spectral_get(y, x);
				if (value.real()*value.real()+value.imag()*value.imag() > 1.0e-13)
					std::cout << "(" << x << ", "<< y << ", "<< value.real() << ", " << value.imag() << ")" <<std::endl;;
			}
		}
	}


	void spectral_structure_print(
			int i_precision = 16,
			double i_abs_threshold = -1
	)	const
	{
		std::cout << std::setprecision(i_precision);

		std::cout << "m \\ n ----->" << std::endl;
		for (std::size_t m = 0; m < planeDataConfig->spectral_data_size[0]; m++)
		{
			//std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
			for (std::size_t n = 0; n < planeDataConfig->spectral_data_size[1]; n++)
			{
				std::cout << "(" << m <<"," << n <<")" << "\t";
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
		file << "#SWEET_PLANE_SPECTRAL_DATA_ASCII" << std::endl;

  		file << "#TI " << i_title << std::endl;

  		// Use 0 to make it processable by python
  		file << "0\t";

  		std::complex<double> w = {0,0};
  		std::vector<double> sum(planeDataConfig->spectral_data_size[1]+1,0);
  		std::vector<double> sum_squared(planeDataConfig->spectral_data_size[1]+1,0);
  		std::vector<double> max(planeDataConfig->spectral_data_size[1]+1,0);

  		for (std::size_t m = 0; m < planeDataConfig->spectral_data_size[0]; m++)
  		{
  			///std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
  			for (std::size_t n = 0; n < planeDataConfig->spectral_data_size[1]; n++)
  			{
  				std::size_t idx = planeDataConfig->getArrayIndexByModes(n, m);
  				w = spectral_space_data[idx];
  				idx++;

  				sum[n]         += std::abs(w);
  				sum_squared[n] += std::abs(w * w);
  				if (std::abs(w) >= max[n]) max[n] = std::abs(w);
  			}
  		}

  		file << planeDataConfig->spectral_data_size[0] << " "
  				<< planeDataConfig->spectral_data_size[1] << std::endl;
  		for (std::size_t n = 0; n <= planeDataConfig->spectral_data_size[1]; n++)
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

		if(i_time == 0.0){
			file.open(i_filename, std::ios_base::trunc);
			file << std::setprecision(i_precision);
			file << "#SWEET_PLANE_SPECTRAL_ABS_EVOL_ASCII" << std::endl;
  			file << "#TI " << i_title << std::endl;
			file << "0\t"<< std::endl; // Use 0 to make it processable by python
			file << "(n_max="<<planeDataConfig->spectral_data_size[1] << " m_max="
					<< planeDataConfig->spectral_data_size[1] << ")" << std::endl;
			file << "timestamp\t" ; 
			for (std::size_t m = 0; m < planeDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
			{
				//std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
				for (std::size_t n = 0; n < planeDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
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
  		for (std::size_t m = 0; m < planeDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
  		{
  			///std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
  			for (std::size_t n = 0; n < planeDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
  			{
  				std::size_t idx = planeDataConfig->getArrayIndexByModes(n, m);
  				w = spectral_space_data[idx];
				wabs = std::abs(w * std::conj(w));
				if ( m > 0 ) sum += 2*wabs; //term appears twice in the spectrum
				else sum += wabs;      // term appears only once
  				
				if ( wabs < i_abs_threshold){
					//file << "(" << n << "," << m << ")\t"<<std::endl;
					file <<  0 << "\t"; //<<std::endl;
				}
				else{
					//file << "(" << n << "," << m << ")\t"<<std::endl;
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

		if(i_time == 0.0){
			file.open(i_filename, std::ios_base::trunc);
			file << std::setprecision(i_precision);
			file << "#SWEET_PLANE_SPECTRAL_PHASE_EVOL_ASCII" << std::endl;
  			file << "#TI " << i_title << std::endl;
			file << "0\t"<< std::endl; // Use 0 to make it processable by python
			file << "(n_max="<<planeDataConfig->spectral_data_size[1] << " m_max="
					<< planeDataConfig->spectral_data_size[1] << ")" << std::endl;
			file << "timestamp\t" ; 
			for (std::size_t m = 0; m < planeDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
			{
				//std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
				for (std::size_t n = 0; n < planeDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
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
		
		file << i_time << "\t";
  		for (std::size_t m = 0; m < planeDataConfig->spectral_data_size[0]/i_reduce_mode_factor; m++)
  		{
  			///std::size_t idx = planeDataConfig->getArrayIndexByModes(m, m);
  			for (std::size_t n = 0; n < planeDataConfig->spectral_data_size[1]/i_reduce_mode_factor; n++)
  			{
  				std::size_t idx = planeDataConfig->getArrayIndexByModes(n, m);
  				w = spectral_space_data[idx];
				wphase = std::arg(w); // std::abs(w * std::conj(w));
				
				file <<  wphase << "\t"; //<<std::endl;;
				
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
  		file << "DATA_TYPE MODES_DATA" << std::endl;
  		file << "MODES_M_MAX " << planeDataConfig->spectral_data_size[0] << std::endl;
  		file << "MODES_N_MAX " << planeDataConfig->spectral_data_size[1] << std::endl;
  		file << "GRID_TYPE GAUSSIAN" << std::endl;
  		file << "NUM_ELEMENTS " << planeDataConfig->spectral_array_data_number_of_elements << std::endl;
  		file << "FIN" << std::endl;

  		file.write((const char*)spectral_space_data, sizeof(std::complex<double>)*planeDataConfig->spectral_array_data_number_of_elements);

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
  		int num_x = -1;
  		int num_y = -1;
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

  			if (buf == "NUM_X")
  			{
  				file >> buf;
  				num_x = std::stoi(buf);
  				std::cout << num_x << std::endl;
  				continue;
  			}

  			if (buf == "NUM_Y")
  			{
  				file >> buf;
  				num_y = std::stoi(buf);
  				std::cout << num_y << std::endl;
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

  		if (data_type != "MODES_DATA")
  			SWEETError("Unknown data type '"+data_type+"'");

  		if (num_x != (int)planeDataConfig->spectral_data_size[0])
  			SWEETError("NUM_X "+std::to_string(num_x)+" doesn't match planeDataConfig");

  		if (num_y != (int)planeDataConfig->spectral_data_size[1])
  			SWEETError("NUM_Y "+std::to_string(num_y)+" doesn't match planeDataConfig");

  		file.read((char*)spectral_space_data, sizeof(std::complex<double>)*planeDataConfig->spectral_array_data_number_of_elements);

  		file.close();
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

		file << "#SWEET_PLANE_SPECTRAL_CPLX_DATA_ASCII" << std::endl;

		size_t ymax = 0;
		if (dimension == 2)
			ymax = planeDataConfig->spectral_data_size[1];
		else
			ymax = 1;

		for (std::size_t y = 0; y < ymax; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				const std::complex<double> &value = spectral_get(y, x);
				file << "(" << value.real() << ", " << value.imag() << ")";

				if (x < planeDataConfig->spectral_data_size[0]-1)
					file << i_separator;
				else
					file << std::endl;
			}
		}

		return true;
	}

	/**
	 * Write amplitude of spectral data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_spectral_abs_saveData_ascii(
			const char *i_filename,		///< Name of file to store data to
			char i_separator = '\t',	///< separator to use for each line
			int i_precision = 12,		///< number of floating point digits
			int dimension = 2			///< store 1D or 2D
	)	const
	{

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		file << "#SWEET_PLANE_SPECTRAL_ABS_DATA_ASCII" << std::endl;

		size_t ymax = 0;
		if (dimension == 2)
			ymax = planeDataConfig->spectral_data_size[1];
		else
			ymax = 1;

		for (std::size_t y = 0; y < ymax; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				file << spectral_return_amplitude(y, x);

				if (x < planeDataConfig->spectral_data_size[0]-1)
					file << i_separator;
				else
					file << std::endl;
			}
		}

		return true;
	}

	/**
	 * Write spectral data to ASCII file
	 *
	 * Each array row is stored to a line.
	 * Per default, a tab separator is used in each line to separate the values.
	 */
	bool file_spectral_arg_saveData_ascii(
			const char *i_filename,		///< Name of file to store data to
			char i_separator = '\t',	///< separator to use for each line
			int i_precision = 12,		///< number of floating point digits
			int dimension = 2			///< store 1D or 2D
	)	const
	{

		std::ofstream file(i_filename, std::ios_base::trunc);
		file << std::setprecision(i_precision);

		size_t ymax = 0;
		if (dimension == 2)
			ymax = planeDataConfig->spectral_data_size[1];
		else
			ymax = 1;

		for (std::size_t y = 0; y < ymax; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{

				file << spectral_return_phase(y, x);

				if (x < planeDataConfig->spectral_data_size[0]-1)
					file << i_separator;
				else
					file << std::endl;
			}
		}

		return true;
	}


	/**
	 * Return average which is given by the first mode
	 */
	inline
	double get_average()	const
	{
		return spectral_space_data[0].real()/(double)planeDataConfig->physical_array_data_number_of_elements;
	}

	inline
	double spectral_return_amplitude(
			std::size_t j,
			std::size_t i
	) const
	{
		std::complex<double> val = spectral_get(j,i);
		val.real(val.real()*2/(planeDataConfig->physical_array_data_number_of_elements));
		val.imag(val.imag()*2/(planeDataConfig->physical_array_data_number_of_elements));

		return std::abs(val);
	}

	inline
	double spectral_return_phase(
			std::size_t j,
			std::size_t i
	) const
	{
		std::complex<double> val = spectral_get(j,i);
		val.real(val.real()*2/(planeDataConfig->physical_array_data_number_of_elements));
		val.imag(val.imag()*2/(planeDataConfig->physical_array_data_number_of_elements));

		return std::arg(val);
	}



	void normalize(
			const std::string &normalization = ""
	)
	{
		if (normalization == "avg_zero")
		{
			// move average value to 0
			double phi_min = getPlaneDataPhysical().physical_reduce_min();
			double phi_max = getPlaneDataPhysical().physical_reduce_max();

			double avg = 0.5*(phi_max+phi_min);

			operator-=(avg);
		}
		else if (normalization == "min_zero")
		{
			// move minimum value to zero
			double phi_min = getPlaneDataPhysical().physical_reduce_min();
			operator-=(phi_min);
		}
		else if (normalization == "max_zero")
		{
			// move maximum value to zero
			double phi_max = getPlaneDataPhysical().physical_reduce_max();
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
	 * Apply a linear operator given by this class to the input data array.
	 */
	inline
	PlaneData_Spectral operator()(
			const PlaneData_Spectral &i_array_data
	)	const
	{
		PlaneData_Spectral out(planeDataConfig);

		PLANE_DATA_SPECTRAL_FOR_IDX(
				out.spectral_space_data[idx] = spectral_space_data[idx]*i_array_data.spectral_space_data[idx];
		);

		out.spectral_zeroAliasingModes();

		return out;
	}


	/**
	 * Interpolate from a finer mesh. Remove highest frequency modes.
	 */
	PlaneData_Spectral restrict(
			const PlaneData_Spectral &i_array_data
	)
	{

		PlaneData_Spectral out = *this;
		out = i_array_data.spectral_returnWithDifferentModes(out.planeDataConfig);

		//////std::size_t M_fine = i_array_data.planeDataConfig->spectral_data_size[0];
		//////std::size_t N_fine = i_array_data.planeDataConfig->spectral_data_size[1];
		//////std::size_t M_coarse = out.planeDataConfig->spectral_data_size[0];
		//////std::size_t N_coarse = out.planeDataConfig->spectral_data_size[1];

		//////assert(M_fine >= M_coarse);
		//////assert(N_fine >= N_coarse);

		//////double rescale =
		//////		(double)(out.planeDataConfig->physical_array_data_number_of_elements)
		//////		/
		//////		(double)(i_array_data.planeDataConfig->physical_array_data_number_of_elements);

		//////// just copy data
		//////if (M_fine == M_coarse && N_fine == N_coarse)
		//////	out = i_array_data;
		//////else
		//////	for (std::size_t m = 0; m < M_coarse; m++)
		//////		for (std::size_t n = 0; n < N_coarse; n++)
		//////		{
		//////			std::size_t idx_coarse = out.planeDataConfig->getArrayIndexByModes(n, m);
		//////			std::size_t idx_fine = i_array_data.planeDataConfig->getArrayIndexByModes(n, m);
		//////			out.spectral_space_data[idx_coarse] = i_array_data.spectral_space_data[idx_fine] * rescale;
		//////		}



		return out;
	}


	/**
	 * Interpolate from a coarser mesh. Pad zeros corresponding to highest frequency modes.
	 */
	PlaneData_Spectral pad_zeros(
			const PlaneData_Spectral &i_array_data
	)
	{

		PlaneData_Spectral out = *this;
		out = i_array_data.spectral_returnWithDifferentModes(out.planeDataConfig);

		///////std::size_t M_coarse = i_array_data.planeDataConfig->spectral_data_size[0];
		///////std::size_t N_coarse = i_array_data.planeDataConfig->spectral_data_size[1];
		///////std::size_t M_fine = out.planeDataConfig->spectral_data_size[0];
		///////std::size_t N_fine = out.planeDataConfig->spectral_data_size[1];

		///////assert(M_fine >= M_coarse);
		///////assert(N_fine >= N_coarse);

		///////double rescale =
		///////		(double)(out.planeDataConfig->physical_array_data_number_of_elements)
		///////		/
		///////		(double)(i_array_data.planeDataConfig->physical_array_data_number_of_elements);

		///////// just copy data
		///////if (M_fine == M_coarse && N_fine == N_coarse)
		///////	out = i_array_data;
		///////else
		///////{
		///////	out.spectral_set_zero();
		///////	for (std::size_t m = 0; m < M_coarse; m++)
		///////		for (std::size_t n = 0; n < N_coarse; n++)
		///////		{
		///////			std::size_t idx_coarse = i_array_data.planeDataConfig->getArrayIndexByModes(n, m);
		///////			std::size_t idx_fine = out.planeDataConfig->getArrayIndexByModes(n, m);
		///////			out.spectral_space_data[idx_fine] = i_array_data.spectral_space_data[idx_coarse] * rescale;
		///////		}
		///////}

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
PlaneData_Spectral operator*(
		double i_value,
		const PlaneData_Spectral &i_array_data
)
{
	return ((PlaneData_Spectral&)i_array_data)*i_value;
}






/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData
 *
 */
inline
static
PlaneData_Spectral operator+(
		double i_value,
		const PlaneData_Spectral &i_array_data
)
{
	return ((PlaneData_Spectral&)i_array_data)+i_value;
}



/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData
 *
 */
inline
static
PlaneData_Spectral operator-(
		double i_value,
		const PlaneData_Spectral &i_array_data
)
{
	return i_array_data.operator_scalar_sub_this(i_value);
}


/**
 * operator to support operations such as:
 *
 * array_data_physical - arrayData
 *
 */
inline
static
PlaneData_Spectral operator-(
		const PlaneData_Physical &i_plane_data_physical,
		const PlaneData_Spectral &i_array_data
)
{
	return - (i_array_data - i_plane_data_physical);
}

}

#endif
