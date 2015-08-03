/*
 * Complex2DArrayFFT.hpp
 *
 *  Created on: 11 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_COMPLEX2DARRAYFFT_HPP_
#define SRC_INCLUDE_SWEET_COMPLEX2DARRAYFFT_HPP_

#if SWEET_USE_SPECTRAL_SPACE
#	include <fftw3.h>
#endif

#include <sweet/DataArray.hpp>
#include <cstddef>
#include <complex>


/**
 * this class is a convenient handler for arrays with complex numbers
 * to convert them to Fourier space forward and backward.
 *
 * In contrast to the DataArray class, this class also supports
 * complex numbers in Cartesian space.
 */
class Complex2DArrayFFT
{
public:
	std::size_t resolution[2];

	fftw_plan plan_to_cart;
	fftw_plan plan_to_spec;

	double *data;

	void setup_fftw()
	{
		// create dummy array for plan creation
		// IMPORTANT! if we use the same array for input/output,
		// a plan will be created with does not support out-of-place
		// FFTs, see http://www.fftw.org/doc/New_002darray-Execute-Functions.html
		double *dummy_data = DataArray<2>::alloc_aligned_mem<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		plan_to_spec =
				fftw_plan_dft_2d(
					resolution[1],	// n0 = ny
					resolution[0],	// n1 = nx
					(fftw_complex*)data,
					(fftw_complex*)dummy_data,
					FFTW_FORWARD,
					0
				);


		plan_to_cart =
				fftw_plan_dft_2d(
					resolution[1],	// n0 = ny
					resolution[0],	// n1 = nx
					(fftw_complex*)data,
					(fftw_complex*)dummy_data,
					FFTW_BACKWARD,
					0
				);

		free(dummy_data);

		if (plan_to_spec == nullptr)
		{
			std::cerr << "Failed to create plan_backward for fftw" << std::endl;
			exit(-1);
		}
	}

public:
	Complex2DArrayFFT(
			const std::size_t i_res[2]
	)	:
		plan_to_cart(nullptr),
		plan_to_spec(nullptr)
	{
		resolution[0] = i_res[0];
		resolution[1] = i_res[1];

		data = DataArray<2>::alloc_aligned_mem<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		setup_fftw();
	}


public:
	Complex2DArrayFFT(
			const Complex2DArrayFFT &i_testArray
	)
	:
		plan_to_cart(nullptr),
		plan_to_spec(nullptr)
	{
		resolution[0] = i_testArray.resolution[0];
		resolution[1] = i_testArray.resolution[1];

		data = DataArray<2>::alloc_aligned_mem<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		setup_fftw();

		memcpy(data, i_testArray.data, sizeof(double)*resolution[0]*resolution[1]*2);
	}


	~Complex2DArrayFFT()
	{
		fftw_free(plan_to_spec);
		fftw_free(plan_to_cart);

		free(data);
	}


public:
	Complex2DArrayFFT& operator=(
			const Complex2DArrayFFT &i_testArray
	)
	{
		resolution[0] = i_testArray.resolution[0];
		resolution[1] = i_testArray.resolution[1];

		memcpy(data, i_testArray.data, sizeof(double)*resolution[0]*resolution[1]*2);
		return *this;
	}


	Complex2DArrayFFT toSpec()
	{
		Complex2DArrayFFT o_testArray(resolution);

		fftw_execute_dft(
				plan_to_spec,
				(fftw_complex*)this->data,
				(fftw_complex*)o_testArray.data
			);

		return o_testArray;
	}


	Complex2DArrayFFT toCart()
	{
		Complex2DArrayFFT o_testArray(resolution);

		fftw_execute_dft(
				plan_to_cart,
				(fftw_complex*)this->data,
				(fftw_complex*)o_testArray.data
			);

		return o_testArray;
	}


	void set(int y, int x, double re, double im)
	{
		data[(y*resolution[0]+x)*2+0] = re;
		data[(y*resolution[0]+x)*2+1] = im;
	}


	void setRe(int y, int x, double re)
	{
		data[(y*resolution[0]+x)*2+0] = re;
	}

	void setIm(int y, int x, double im)
	{
		data[(y*resolution[0]+x)*2+1] = im;
	}

	double getRe(int y, int x)	const
	{
		return data[(y*resolution[0]+x)*2+0];
	}

	double getIm(int y, int x)	const
	{
		return data[(y*resolution[0]+x)*2+1];
	}

	void setAll(double re, double im)
	{
		for (std::size_t y = 0; y < resolution[1]; y++)
			for (std::size_t x = 0; x < resolution[0]; x++)
				set(y, x, re, im);
	}

	void setAllRe(double re)
	{
		for (std::size_t y = 0; y < resolution[1]; y++)
			for (std::size_t x = 0; x < resolution[0]; x++)
				setRe(y, x, re);
	}

	void setAllIm(double im)
	{
		for (std::size_t y = 0; y < resolution[1]; y++)
			for (std::size_t x = 0; x < resolution[0]; x++)
				setIm(y, x, im);
	}


	/**
	 * return |u|^2 with
	 */
	double getNormSquared(int y, int x)	const
	{
		double re = data[(y*resolution[0]+x)*2+0];
		double im = data[(y*resolution[0]+x)*2+1];

		return std::abs(re*re+im*im);
	}


	/**
	 * compute
	 *
	 * int(x*value(x), x=0..inf)/int(x, x=0..inf)
	 */
	double return_centroidFromEnergy()
	{
		double nominator = 0;
		double denominator = 0;

		for (std::size_t kb = 0; kb < resolution[1]; kb++)
		{
//			for (std::size_t ka = 0; ka < resolution[0]; ka++)
//			{
				double k = std::sqrt((double)(kb*kb));
//				double k = std::sqrt((double)(ka*ka)+(double)(kb*kb));
//				double k = ka+kb;

//				double energy_ampl = std::sqrt(getNormSquared(kb, ka));
				double energy_ampl = std::sqrt(getNormSquared(kb, 0));

				nominator += k*energy_ampl;
				denominator += energy_ampl;
//			}
		}

		std::cout << nominator << std::endl;
		return nominator/denominator;
//		return nominator/(double)(resolution[0]*resolution[1]);
	}

	/**
	 * return the energy centroid with the velocity components given
	 * by this class and the parameter
	 */
	double return_energyCentroid(
			Complex2DArrayFFT &v
	)
	{
		Complex2DArrayFFT &u = *this;

		double nominator = 0;
		double denominator = 0;

		for (std::size_t kb = 0; kb < resolution[1]; kb++)
			for (std::size_t ka = 0; ka < resolution[0]; ka++)
			{
				double k = std::sqrt((double)(ka*ka)+(double)(kb*kb));

				double u2 = u.getNormSquared(kb, ka);
				double v2 = v.getNormSquared(kb, ka);

				double ampl = u2+v2;
				nominator += k*ampl;
				denominator += ampl;
			}

		return nominator/denominator;
	}


	/**
	 * return the maximum of all absolute values
	 */
	double reduce_sumAbs()	const
	{
		double sum = 0;
		for (std::size_t i = 0; i < resolution[0]*resolution[1]; i++)
			sum += std::abs(data[i*2])+std::abs(data[i*2+1]);

		return sum;
	}


	friend
	inline
	std::ostream& operator<<(std::ostream &o_ostream, const Complex2DArrayFFT &i_testArray)
	{
		for (int y = i_testArray.resolution[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < i_testArray.resolution[0]; x++)
			{
				double value_re = i_testArray.getRe(y, x);
				double value_im = i_testArray.getIm(y, x);
				o_ostream << "(" << value_re << ", " << value_im << ")\t";
			}
			o_ostream << std::endl;
		}
		return o_ostream;
	}



	Complex2DArrayFFT spec_div_element_wise(
			Complex2DArrayFFT &i_d
	)
	{
		Complex2DArrayFFT out(resolution);

		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double ar = data[i];
			double ai = data[i+1];
			double br = i_d.data[i];
			double bi = i_d.data[i+1];

			double den = (br*br+bi*bi);

			if (std::abs(den) == 0)
			{
				// For Laplace solution, this is the integration constant C
				out.data[i] = 0;
				out.data[i+1] = 0;
			}
			else
			{
				out.data[i] = (ar*br + ai*bi)/den;
				out.data[i+1] = (ai*br - ar*bi)/den;
			}
		}

		return out;
	}



	DataArray<2> getRealWithDataArray()
	{
		DataArray<2> out(resolution);

		for (std::size_t j = 0; j < resolution[1]; j++)
			for (std::size_t i = 0; i < resolution[0]; i++)
				out.set(j, i, getRe(j, i));

		return out;
	}



	DataArray<2> getImagWithDataArray()
	{
		DataArray<2> out(resolution);

		for (std::size_t j = 0; j < resolution[1]; j++)
			for (std::size_t i = 0; i < resolution[0]; i++)
				out.set(j, i, getIm(j, i));

		return out;
	}

	Complex2DArrayFFT &loadRealFromDataArray(
			const DataArray<2> &i_dataArray_Real
	)
	{
		i_dataArray_Real.requestDataInCartesianSpace();

		for (std::size_t j = 0; j < resolution[1]; j++)
			for (std::size_t i = 0; i < resolution[0]; i++)
				set(j, i, i_dataArray_Real.get(j, i), 0);

		return *this;
	}


	Complex2DArrayFFT &loadRealAndImagFromDataArrays(
			const DataArray<2> &i_dataArray_Real,
			const DataArray<2> &i_dataArray_Imag
	)
	{
		i_dataArray_Real.requestDataInCartesianSpace();
		i_dataArray_Imag.requestDataInCartesianSpace();

		for (std::size_t j = 0; j < resolution[1]; j++)
			for (std::size_t i = 0; i < resolution[0]; i++)
				set(	j, i,
						i_dataArray_Real.get(j, i),
						i_dataArray_Imag.get(j, i)
				);

		return *this;
	}



	/**
	 * Compute multiplication with a complex scalar
	 */
	inline
	Complex2DArrayFFT operator*(
			const std::complex<double> &i_value
	)	const
	{
		Complex2DArrayFFT out(resolution);

		double br = i_value.real();
		double bi = i_value.imag();

		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double ar = data[i];
			double ai = data[i+1];

			out.data[i] = ar*br - ai*bi;
			out.data[i+1] = ar*bi + ai*br;
		}
		return out;
	}



	/**
	 * Compute element-wise addition
	 */
	inline
	Complex2DArrayFFT operator+(
			const Complex2DArrayFFT &i_array_data
	)	const
	{
		Complex2DArrayFFT out(this->resolution);

		#pragma omp parallel for OPENMP_SIMD
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i++)
		{
			out.data[i] = data[i] + i_array_data.data[i];
			out.data[i+1] = data[i+1] + i_array_data.data[i+1];
		}

		return out;
	}

};


inline
static
Complex2DArrayFFT operator*(
		const std::complex<double> &i_value,
		const Complex2DArrayFFT &i_array_data
)
{
	return ((Complex2DArrayFFT&)i_array_data)*i_value;
}

#endif /* SRC_INCLUDE_SWEET_COMPLEX2DARRAYFFT_HPP_ */
