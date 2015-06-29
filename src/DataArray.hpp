/*
 * DataArray.hpp
 *
 *  Created on: 28 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_DATAARRAY_HPP_
#define SRC_DATAARRAY_HPP_



#include <cassert>
#include <cstddef>
#include <cmath>
#include <memory>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <utility>

#include <fftw3.h>
#include <omp.h>

/**
 * allocator which allocated memory blocks aligned at 128 byte boundaries
 */
template <typename T=void>
T *alloc_aligned_mem(
		std::size_t i_size
)
{
	T *data;
	int retval = posix_memalign((void**)&data, 128, i_size);
	if (retval != 0)
	{
		std::cerr << "Unable to allocate memory" << std::endl;
		assert(false);
		exit(-1);
	}
	return data;
}



/**
 * Data array and operator support.
 *
 * Here, we assume the Cartesian coordinate system given similar to the following sketch:
 *
 *  Y ^
 *    |
 *    |
 *    |
 *    +-------->
 *           X
 *
 * Also the arrays are stored in this way:
 * 		A[Y0...YN-1][X0...XN-1]
 *
 */
template <int D>
class DataArray
{
public:
	/**
	 * global size of allocated array
	 * (x,y[,z])
	 */
	std::size_t resolution[D];

	/**
	 * local data in cartesian space
	 */
	std::size_t array_data_cartesian_length;
	double *array_data_cartesian_space;
	bool array_data_cartesian_space_valid;

	/**
	 * local data in spectral space
	 */
	std::size_t array_data_spectral_length;
	double *array_data_spectral_space;
	bool array_data_spectral_space_valid;

	/**
	 * local ranges
	 */
	std::size_t range_start[D];
	std::size_t range_end[D];
	std::size_t range_size[D];

	/**
	 * temporary data?
	 *
	 * Temporary data can be created if e.g. operators are evaluated:
	 * h = hu/h
	 *
	 * This first creates a temporary DataArray to compute hu/h.
	 *
	 * This is then followed by an assignment of this data to h.
	 */
	bool temporary_data;

/*

private:
	// http://stackoverflow.com/questions/124856/how-do-i-prevent-a-class-from-being-allocated-via-the-new-operator-id-like
	// Prevent heap allocation
	void *operator new(std::size_t);
	void *operator new[](std::size_t);
	void operator delete(void*);
	void operator delete[](void*);
*/

private:
	/**
	 * prohobit empty initialization
	 */
	DataArray()	{}



public:
	/**
	 * copy constructor, used e.g. in
	 * 	DataArray<2> tmp_h = h;
	 * 	DataArray<2> tmp_h2(h);
	 *
	 * Duplicate all data
	 */
	DataArray(
			const DataArray<D> &i_dataArray
	)	:
		temporary_data(false)
	{
//		std::cout << "Copy constructor called" << std::endl;

		for (int i = 0; i < D; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;
		array_data_cartesian_space_valid = i_dataArray.array_data_cartesian_space_valid;
		array_data_cartesian_space = alloc_aligned_mem<double>(array_data_cartesian_length*sizeof(double));
		if (array_data_cartesian_space_valid)
			memcpy(array_data_cartesian_space, i_dataArray.array_data_cartesian_space, array_data_cartesian_length*sizeof(double));

		array_data_spectral_length = i_dataArray.array_data_spectral_length;
		array_data_spectral_space_valid = i_dataArray.array_data_spectral_space_valid;
		array_data_spectral_space = alloc_aligned_mem<double>(array_data_spectral_length*sizeof(double));
		if (array_data_spectral_space_valid)
			memcpy(array_data_spectral_space, i_dataArray.array_data_spectral_space, array_data_cartesian_length*sizeof(double));


		auto fft_ptr = *fftGetSingletonPtr();
		fft_ptr->ref_counter++;
	}


public:
	/**
	 * Move constructor
	 */
	DataArray(
			DataArray<D> &&i_dataArray
	)	:
		temporary_data(false)
	{
//		std::cout << "Move constructor called" << std::endl;

		for (int i = 0; i < D; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;
		array_data_cartesian_space = i_dataArray.array_data_cartesian_space;
		i_dataArray.array_data_cartesian_space = nullptr;
		array_data_cartesian_space_valid = i_dataArray.array_data_cartesian_space_valid;
		i_dataArray.array_data_cartesian_space_valid = false;

		array_data_spectral_length = i_dataArray.array_data_spectral_length;
		array_data_spectral_space = i_dataArray.array_data_spectral_space;
		i_dataArray.array_data_spectral_space = nullptr;
		array_data_spectral_space_valid = i_dataArray.array_data_spectral_space_valid;
		i_dataArray.array_data_spectral_space_valid = false;

		auto fft_ptr = *fftGetSingletonPtr();
		fft_ptr->ref_counter++;
	}



	/**
	 * default constructor
	 */
public:
	DataArray(
		const std::size_t i_resolution[D]	///< size of array
	)	:
		temporary_data(false)
	{
//		std::cout << "Std constructor called" << std::endl;

		array_data_cartesian_length = 1;
		for (int i = 0; i < D; i++)
		{
			array_data_cartesian_length *= i_resolution[i];

			resolution[i] = i_resolution[i];
			range_start[i] = 0;
			range_end[i] = i_resolution[i];
			range_size[i] = range_end[i]-range_start[i];
		}


		array_data_cartesian_space = alloc_aligned_mem<double>(array_data_cartesian_length*sizeof(double));
		array_data_cartesian_space_valid = false;

		array_data_spectral_length = array_data_cartesian_length/i_resolution[D-1];	/// see FFTW documentation for allocation of memory buffers
		array_data_spectral_length *= (i_resolution[D-1]/2+1)*2;
		array_data_spectral_space = alloc_aligned_mem<double>(array_data_spectral_length*sizeof(double));
		array_data_spectral_space_valid = false;

		// initialize fft if not yet done
		fftTestAndInit(resolution);

		auto fft_ptr = *fftGetSingletonPtr();
		fft_ptr->ref_counter++;
	}

	~DataArray()
	{
//		std::cout << "Deconstructor called" << std::endl;
		free(array_data_cartesian_space);
		free(array_data_spectral_space);

		auto fft_ptr = *fftGetSingletonPtr();
		fft_ptr->ref_counter--;

		assert(fft_ptr->ref_counter >= 0);
		if (fft_ptr->ref_counter == 0)
		{
			delete *fftGetSingletonPtr();
			*fftGetSingletonPtr() = nullptr;
		}
	}


	inline
	double &getDataRef(
			std::size_t j,
			std::size_t i
	)
	{
		assert(i >= range_start[0] && i < range_end[0]);
		assert(j >= range_start[1] && j < range_end[1]);

		return array_data_cartesian_space[
							(j-range_start[1])*range_size[0]+
							(i-range_start[0])
						];
	}


	inline
	double &getDataRef(
			std::size_t k,
			std::size_t j,
			std::size_t i
	)
	{
		assert(i >= range_start[0] && i < range_end[0]);
		assert(j >= range_start[1] && j < range_end[1]);
		assert(k >= range_start[2] && k < range_end[2]);

		return array_data_cartesian_space[
						   (k-range_start[2])*range_size[1]*range_size[0]+
						   (j-range_start[1])*range_size[0]+
						   (i-range_start[0])
					  ];
	}


	inline
	void data_setall(
			double i_value
	)
	{
#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] = i_value;

		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;
	}



	class FFTWSingletonClass
	{
	public:
		int ref_counter;

	private:
		fftw_plan	plan_forward;
		std::size_t plan_forward_output_length;

		fftw_plan	plan_backward;
		std::size_t plan_backward_output_length;

	public:
		FFTWSingletonClass(
				DataArray<D> &i_dataArray
		)	:
			ref_counter(0)
		{
			// support threading
			fftw_init_threads();

		    fftw_plan_with_nthreads(omp_get_max_threads());

			plan_backward_output_length = i_dataArray.array_data_cartesian_length;
			plan_forward_output_length = i_dataArray.array_data_spectral_length;

			double *data_cartesian = alloc_aligned_mem<double>(i_dataArray.array_data_cartesian_length*sizeof(double));
			double *data_spectral = alloc_aligned_mem<double>(i_dataArray.array_data_spectral_length*sizeof(double));

			plan_forward =
					fftw_plan_dft_r2c_2d(
						i_dataArray.resolution[1],
						i_dataArray.resolution[0],
						data_cartesian,
						(fftw_complex*)data_spectral,
						FFTW_PRESERVE_INPUT
					);

			if (plan_forward == nullptr)
			{
				std::cerr << "Failed to create forward plan for fftw" << std::endl;
				exit(-1);
			}

			plan_backward =
					fftw_plan_dft_c2r_2d(
						i_dataArray.resolution[1],
						i_dataArray.resolution[0],
						(fftw_complex*)data_spectral,
						data_cartesian,
						0
					);

			if (plan_backward == nullptr)
			{
				std::cerr << "Failed to create backward plan for fftw" << std::endl;
				exit(-1);
			}

			free(data_cartesian);
			free(data_spectral);
		}

		void fft_forward(
				DataArray<D> &io_dataArray
		)
		{
			fftw_execute_dft_r2c(plan_forward, io_dataArray.array_data_cartesian_space, (fftw_complex*)io_dataArray.array_data_spectral_space);
		}

		void fft_backward(
				DataArray<D> &io_dataArray
		)
		{
			fftw_execute_dft_c2r(plan_backward, (fftw_complex*)io_dataArray.array_data_spectral_space, io_dataArray.array_data_cartesian_space);
			// spectral data is not valid anymore, since c2r is destructive!
			io_dataArray.array_data_spectral_space_valid = false;

			double scale = (1.0/(double)plan_backward_output_length);
#pragma omp parallel for simd
			for (std::size_t i = 0; i < plan_backward_output_length; i++)
				io_dataArray.array_data_cartesian_space[i] *= scale;
		}

		~FFTWSingletonClass()
		{
			fftw_destroy_plan(plan_forward);
			fftw_destroy_plan(plan_backward);

			fftw_cleanup();
			fftw_cleanup_threads();
		}
	};


private:
	FFTWSingletonClass** fftGetSingletonPtr()
	{
		static FFTWSingletonClass *fftw_singleton_data = nullptr;
		return &fftw_singleton_data;
	}


private:
	FFTWSingletonClass* fftTestAndInit(
		std::size_t i_size[D]
	)
	{
		FFTWSingletonClass **fftw_singleton_data = fftGetSingletonPtr();

		if (*fftw_singleton_data != nullptr)
			return *fftw_singleton_data;

		*fftw_singleton_data = new FFTWSingletonClass(*this);

		return *fftw_singleton_data;
	}


public:
	void requestDataInSpectralSpace()
	{
		if (array_data_spectral_space_valid)
			return;		// nothing to do

		assert(array_data_cartesian_space_valid == true);

		(*fftGetSingletonPtr())->fft_forward(*this);

		array_data_spectral_space_valid = true;
	}


	void requestDataInCartesianSpace()
	{
		if (array_data_cartesian_space_valid)
			return;		// nothing to do

		assert(array_data_spectral_space_valid == true);

		(*fftGetSingletonPtr())->fft_backward(*this);

		array_data_cartesian_space_valid = true;
	}


	/**
	 * return the maximum of all absolute values
	 */
	double get_maxAbs()
	{
		requestDataInCartesianSpace();

		double maxabs = -1;
#pragma omp parallel for simd reduction(max:maxabs)
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			maxabs = std::max(maxabs, std::abs(array_data_cartesian_space[i]));

		return maxabs;
	}


public:
	template <int S>
	void setup_kernel(
			const double i_kernel_array[S][S]
	)
	{
		double inv_kernel_array[S][S];

		for (int j = 0; j < S; j++)
			for (int i = 0; i < S; i++)
				inv_kernel_array[j][i] = i_kernel_array[S-j-1][S-i-1];

		assert(D == 2);

		// assure symmetric kernel
		assert((S & 1) == 1);

		// radius of kernel (half size)
		std::size_t R = S>>1;

		// zero data in Cartesian space
#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] = 0;

		// left lower corner
		//kernel_cart[    0:R[0]+1,       0:R[0]+1    ] = conv_kernel[    R[0]:,  R[0]:   ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = ky-R;
				std::size_t cx = kx-R;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				getDataRef(cy, cx) = inv_kernel_array[ky][kx];
			}


		// right bottom corner
		//kernel_cart[    0:R[0]+1,       res-R[0]:   ] = conv_kernel[    R[0]:,  0:R[0]  ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = ky-R;
				std::size_t cx = resolution[0] - R + kx;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				getDataRef(cy, cx) = inv_kernel_array[ky][kx];
			}


		// left top corner
		//kernel_cart[    res-R[0]:,      0:R[0]+1   ] = conv_kernel[    0:R[0],  R[0]:  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = resolution[1] - R + ky;
				std::size_t cx = kx-R;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				getDataRef(cy, cx) = inv_kernel_array[ky][kx];
			}


		// right top corner
		//kernel_cart[    res-R[0]:,      res-R[0]:   ] = conv_kernel[    0:R[0], 0:R[0]  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				std::size_t cy = resolution[1] - R + ky;
				std::size_t cx = resolution[0] - R + kx;

				if (range_start[0] > cx || range_end[0] <= cx)
					continue;
				if (range_start[1] > cy || range_end[1] <= cy)
					continue;

				getDataRef(cy, cx) = inv_kernel_array[ky][kx];
			}

		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;

		requestDataInSpectralSpace();	/// convert kernel to spectral space
	}


public:
	template <int S>
	void setup_kernel(
			const double i_kernel_array[S][S][S]
	)
	{
		assert(D == 3);
		assert(false);//TODO
	}


public:
	/**
	 * assignment operator
	 *
	 * hasdfasdf = h*hasdf;
	 */
	DataArray<D> &operator=(
			const DataArray<D> &i_dataArray
	)
	{
//		std::cout << "Assignment operator called" << std::endl;

		for (int i = 0; i < D; i++)
		{
			resolution[i] = i_dataArray.resolution[i];

			range_start[i] = i_dataArray.range_start[i];
			range_end[i] = i_dataArray.range_end[i];
			range_size[i] = i_dataArray.range_size[i];
		}

		array_data_cartesian_length = i_dataArray.array_data_cartesian_length;
		array_data_spectral_length = i_dataArray.array_data_spectral_length;

		if (i_dataArray.array_data_cartesian_space_valid)
		{
			array_data_cartesian_space_valid = true;

			/**
			 * If this data was generated based on temporary data sets (e.g. via h = hu/u), then only swap pointers.
			 */
			if (i_dataArray.temporary_data)
				std::swap(array_data_cartesian_space, ((DataArray<D> &)i_dataArray).array_data_cartesian_space);
			else
				memcpy(array_data_cartesian_space, i_dataArray.array_data_cartesian_space, array_data_cartesian_length*sizeof(double));
		}
		else
		{
			array_data_cartesian_space_valid = false;
		}

		if (i_dataArray.array_data_spectral_space_valid)
		{
			array_data_spectral_space_valid = true;
			if (i_dataArray.temporary_data)
				std::swap(array_data_spectral_space, ((DataArray<D> &)i_dataArray).array_data_spectral_space);
			else
				memcpy(array_data_spectral_space, i_dataArray.array_data_spectral_space, array_data_spectral_length*sizeof(double));
		}
		else
		{
			array_data_spectral_space_valid = false;
		}

		return *this;
	}


	/**
	 * Apply a linear operator given by this class to the input data array.
	 */
	inline
	DataArray<D> operator()(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		DataArray<D> out  = DataArray<D>(this->resolution);
		out.temporary_data = true;

		requestDataInSpectralSpace();
		rw_array_data.requestDataInSpectralSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_spectral_length; i+=2)
		{
			double ar = array_data_spectral_space[i];
			double ai = array_data_spectral_space[i+1];
			double br = i_array_data.array_data_spectral_space[i];
			double bi = i_array_data.array_data_spectral_space[i+1];

			out.array_data_spectral_space[i] = ar*br - ai*bi;
			out.array_data_spectral_space[i+1] = ar*bi + ai*br;
		}

		out.array_data_spectral_space_valid = true;
		out.array_data_cartesian_space_valid = false;

		return out;
	}



	/**
	 * Compute element-wise addition
	 */
	inline
	DataArray<D> operator+(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		auto out  = DataArray<D>(this->resolution);
		out.temporary_data = true;

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]+
					i_array_data.array_data_cartesian_space[i];

		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
		return out;
	}




	/**
	 * Compute element-wise addition
	 */
	inline
	DataArray<D> operator+(
			const double i_value
	)
	{
		auto out  = DataArray<D>(this->resolution);
		out.temporary_data = true;

		requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]+i_value;

		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
		return out;
	}




	/**
	 * Compute element-wise addition
	 */
	inline
	DataArray<D>& operator+=(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] +=
					i_array_data.array_data_cartesian_space[i];

		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;
		return *this;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	DataArray<D>& operator-=(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] -=
					i_array_data.array_data_cartesian_space[i];

		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;
		return *this;
	}


	/**
	 * Compute element-wise subtraction
	 */
	inline
	DataArray<D> operator-(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		auto out  = DataArray<D>(this->resolution);
		out.temporary_data = true;

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]-
					i_array_data.array_data_cartesian_space[i];

		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
		return out;
	}



	/**
	 * Compute element-wise subtraction
	 */
	inline
	DataArray<D> operator-(
			const double i_value
	)
	{
		auto out  = DataArray<D>(this->resolution);
		out.temporary_data = true;

		requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]-i_value;

		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
		return out;
	}


	/**
	 * Invert sign
	 */
	inline
	DataArray<D>& operator-()
	{
		requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			array_data_cartesian_space[i] = -array_data_cartesian_space[i];

		array_data_cartesian_space_valid = true;
		array_data_spectral_space_valid = false;

		return *this;
	}


	/**
	 * Compute element-wise multiplication
	 */
	inline
	DataArray<D> operator*(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		DataArray<D> out  = DataArray<D>(i_array_data.resolution);
		out.temporary_data = true;

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]*
					i_array_data.array_data_cartesian_space[i];

		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;

		return out;
	}

	/**
	 * Compute element-wise multiplication with a scalar
	 */
	inline
	DataArray<D> operator*(
			const double i_value
	)
	{
		DataArray<D> out  = DataArray<D>(resolution);
		out.temporary_data = true;

		if (array_data_spectral_space_valid)
		{
#pragma omp parallel for simd
			for (std::size_t i = 0; i < array_data_spectral_length; i++)
				out.array_data_spectral_space[i] =
						array_data_spectral_space[i]*i_value;

			out.array_data_cartesian_space_valid = false;
			out.array_data_spectral_space_valid = true;
			return out;
		}


		if (array_data_cartesian_space_valid)
		{
#pragma omp parallel for simd
			for (std::size_t i = 0; i < array_data_cartesian_length; i++)
				out.array_data_cartesian_space[i] =
						array_data_cartesian_space[i]*i_value;

			out.array_data_spectral_space_valid = false;
			out.array_data_cartesian_space_valid = true;
			return out;
		}

		assert(false);
		return out;
	}


	/**
	 * Compute element-wise division
	 */
	inline
	DataArray<D> operator/(
			const DataArray<D> &i_array_data
	)
	{
		DataArray<D> &rw_array_data = (DataArray<D>&)i_array_data;

		auto out  = DataArray<D>(this->resolution);
		out.temporary_data = true;

		requestDataInCartesianSpace();
		rw_array_data.requestDataInCartesianSpace();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < array_data_cartesian_length; i++)
			out.array_data_cartesian_space[i] =
					array_data_cartesian_space[i]/
					i_array_data.array_data_cartesian_space[i];

		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
		return out;
	}



	friend
	inline
	std::ostream& operator<<(std::ostream &o_ostream, DataArray<D> &i_dataArray)
	{
		i_dataArray.requestDataInCartesianSpace();

		assert(D == 2);
		if (D == 2)
		{
			for (int y = i_dataArray.resolution[1]-1; y >= 0; y--)
			{
				for (std::size_t x = 0; x < i_dataArray.resolution[0]; x++)
				{
					double value = i_dataArray.getDataRef(y, x);
//					if (std::abs(value) < 1e-13)
//						value = 0;
					std::cout << value << "\t";
				}
				std::cout << std::endl;
			}
		}
		return o_ostream;
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
DataArray<2> operator*(
		const double i_value,
		const DataArray<2> &i_array_data
)
{
	return ((DataArray<2>&)i_array_data)*i_value;
}

#endif /* SRC_DATAARRAY_HPP_ */
