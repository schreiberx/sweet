/*
 * Complex2DArrayFFT.hpp
 *
 *  Created on: 11 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_COMPLEX2DARRAYFFT_HPP_
#define SRC_INCLUDE_SWEET_COMPLEX2DARRAYFFT_HPP_


#include <fftw3.h>
#include <sweet/DataArray.hpp>
#include <sweet/NUMABlockAlloc.hpp>
#include <cstddef>
#include <complex>
#include <sweet/openmp_helper.hpp>



/**
 * this class is a convenient handler for arrays with complex numbers
 * to convert them to Fourier space forward and backward.
 *
 * In contrast to the DataArray class, this class also supports
 * complex numbers in Cartesian space.
 */
class Complex2DArrayFFT
{
	typedef std::complex<double> complex;

public:
	/**
	 * resolution of data array
	 */
	std::size_t resolution[2];

	/**
	 * Is this an aliased version?
	 * Then this uses different FFTW plans
	 */
	bool aliased_scaled = false;

	/**
	 * Data associated to this complex array
	 */
	double *data;

	/**
	 * Flag which is true if this is the initializer class of the FFTW library
	 */
	bool is_fft_data_initialized;


	/**
	 * container for different plans supported by this class
	 */
	class Plans
	{
	public:
		fftw_plan to_cart;
		fftw_plan to_spec;
		fftw_plan to_cart_aliasing;
		fftw_plan to_spec_aliasing;

		std::size_t resolution[2];
		std::size_t resolution_aliasing[2];
	};


	/**
	 * Singleton of plans
	 */
	Plans &fft_getSingleton_Plans()	const
	{
		static Plans plans;
		return plans;
	}


	/**
	 * Singleton of reference counter
	 */
	int &fft_getSingleton_RefCounter()
	{
		static int ref_counter = 0;
		return ref_counter;
	}



	/**
	 * Setup the FFTW
	 *
	 * This may be only called once
	 */
	void fft_setup()
	{
		assert(!is_fft_data_initialized);
		is_fft_data_initialized = true;

		assert(fft_getSingleton_RefCounter() >= 0);

		int &ref_counter = fft_getSingleton_RefCounter();

#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp atomic
#endif

		ref_counter++;

		if (ref_counter != 1)
			return;

		bool wisdom_loaded = DataArray<2>::FFTWSingletonClass::loadWisdom();

		{
			// create dummy array for plan creation
			// IMPORTANT! if we use the same array for input/output,
			// a plan will be created with does not support out-of-place
			// FFTs, see http://www.fftw.org/doc/New_002darray-Execute-Functions.html
			double *dummy_data = NUMABlockAlloc::alloc<double>(sizeof(double)*resolution[0]*resolution[1]*2);

			fft_getSingleton_Plans().to_spec =
					fftw_plan_dft_2d(
						resolution[1],	// n0 = ny
						resolution[0],	// n1 = nx
						(fftw_complex*)data,
						(fftw_complex*)dummy_data,
						FFTW_FORWARD,
						(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
					);

			if (fft_getSingleton_Plans().to_spec == nullptr)
			{
				std::cerr << "Failed to create plan_forward for fftw" << std::endl;
				std::cerr << "complex forward preverse_input forward " << resolution[0] << " x " << resolution[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << resolution[0] << "x" << resolution[1] << std::endl;
				exit(-1);
			}

			fft_getSingleton_Plans().to_cart =
					fftw_plan_dft_2d(
						resolution[1],	// n0 = ny
						resolution[0],	// n1 = nx
						(fftw_complex*)data,
						(fftw_complex*)dummy_data,
						FFTW_BACKWARD,
						(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
					);

			fft_getSingleton_Plans().resolution[0] = resolution[0];
			fft_getSingleton_Plans().resolution[1] = resolution[1];

			if (fft_getSingleton_Plans().to_cart == nullptr)
			{
				std::cerr << "Failed to create plan_backward for fftw" << std::endl;
				std::cerr << "complex backward preverse_input forward " << resolution[0] << " x " << resolution[1] << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << resolution[0] << "x" << resolution[1] << std::endl;
				exit(-1);
			}

			NUMABlockAlloc::free(dummy_data, sizeof(double)*resolution[0]*resolution[1]*2);
		}

		{
			double *dummy_data_aliasing_in = NUMABlockAlloc::alloc<double>(sizeof(double)*resolution[0]*resolution[1]*2*4);
			double *dummy_data_aliasing_out = NUMABlockAlloc::alloc<double>(sizeof(double)*resolution[0]*resolution[1]*2*4);

			fft_getSingleton_Plans().to_spec_aliasing =
					fftw_plan_dft_2d(
						resolution[1]*2,	// n0 = ny
						resolution[0]*2,	// n1 = nx
						(fftw_complex*)dummy_data_aliasing_in,
						(fftw_complex*)dummy_data_aliasing_out,
						FFTW_FORWARD,
						(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
					);

			if (fft_getSingleton_Plans().to_spec_aliasing == nullptr)
			{
				std::cerr << "Failed to create plan_forward for aliasing fftw" << std::endl;
				std::cerr << "complex backward preverse_input forward " << resolution[0]*2 << " x " << resolution[1]*2 << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << resolution[0]*2 << "x" << resolution[1]*2 << std::endl;
				exit(-1);
			}


			fft_getSingleton_Plans().to_cart_aliasing =
					fftw_plan_dft_2d(
						resolution[1]*2,	// n0 = ny
						resolution[0]*2,	// n1 = nx
						(fftw_complex*)dummy_data_aliasing_out,
						(fftw_complex*)dummy_data_aliasing_in,
						FFTW_BACKWARD,
						(!wisdom_loaded ? FFTW_PRESERVE_INPUT : FFTW_PRESERVE_INPUT | FFTW_WISDOM_ONLY)
					);

			fft_getSingleton_Plans().resolution_aliasing[0] = resolution[0]*2;
			fft_getSingleton_Plans().resolution_aliasing[1] = resolution[1]*2;


			if (fft_getSingleton_Plans().to_cart_aliasing == nullptr)
			{
				std::cerr << "Failed to create plan_backward for aliasing fftw" << std::endl;
				std::cerr << "complex backward preverse_input forward " << resolution[0]*2 << " x " << resolution[1]*2 << std::endl;
				std::cerr << "fftw-wisdom plan: cf" << resolution[0]*2 << "x" << resolution[1]*2 << std::endl;
				exit(-1);
			}


			NUMABlockAlloc::free(dummy_data_aliasing_out, sizeof(double)*resolution[0]*resolution[1]*2*4);
			NUMABlockAlloc::free(dummy_data_aliasing_in, sizeof(double)*resolution[0]*resolution[1]*2*4);
		}
	}

	void fftw_shutdown()
	{
		if (!is_fft_data_initialized)
			return;

		is_fft_data_initialized = false;

		int &ref_counter = fft_getSingleton_RefCounter();

#if SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp atomic
#endif
		ref_counter--;

		if (ref_counter > 0)
			return;

		assert(ref_counter >= 0);

		fftw_destroy_plan(fft_getSingleton_Plans().to_spec);
		fftw_destroy_plan(fft_getSingleton_Plans().to_cart);

		fftw_destroy_plan(fft_getSingleton_Plans().to_spec_aliasing);
		fftw_destroy_plan(fft_getSingleton_Plans().to_cart_aliasing);

		fft_getSingleton_Plans().resolution[0] = 0;
		fft_getSingleton_Plans().resolution[1] = 0;

		fft_getSingleton_Plans().resolution_aliasing[0] = 0;
		fft_getSingleton_Plans().resolution_aliasing[1] = 0;
	}



	Complex2DArrayFFT()	:
		data(nullptr),
		is_fft_data_initialized(false)
	{

#if !SWEET_USE_LIBFFT
		std::cerr << "This class only makes sense with FFT" << std::endl;
		exit(1);
#endif

	}



public:
	Complex2DArrayFFT(
			const std::size_t i_res[2],
			bool i_aliased_scaled = false
	)	:
		is_fft_data_initialized(false)
	{
#if !SWEET_USE_LIBFFT
		std::cerr << "This class only makes sense with FFT" << std::endl;
		exit(1);
#endif

		aliased_scaled = i_aliased_scaled;

		resolution[0] = i_res[0];
		resolution[1] = i_res[1];

		data = NUMABlockAlloc::alloc<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		fft_setup();
	}



public:
	void setup(
			const std::size_t i_res[2],
			bool i_aliased_scaled = false
	)
	{
		aliased_scaled = i_aliased_scaled;

		resolution[0] = i_res[0];
		resolution[1] = i_res[1];

		if (data)
			cleanup();

		data = NUMABlockAlloc::alloc<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		if (!is_fft_data_initialized)
			fft_setup();
	}



	void cleanup()
	{
		NUMABlockAlloc::free(data, sizeof(double)*resolution[0]*resolution[1]*2);
		data = nullptr;
	}



public:
	Complex2DArrayFFT(
			const Complex2DArrayFFT &i_testArray
	)	:
		is_fft_data_initialized(false)
	{
#if !SWEET_USE_LIBFFT
		std::cerr << "This class only makes sense with FFT" << std::endl;
		exit(1);
#endif

		resolution[0] = i_testArray.resolution[0];
		resolution[1] = i_testArray.resolution[1];

		aliased_scaled = i_testArray.aliased_scaled;

		data = NUMABlockAlloc::alloc<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		fft_setup();

		par_doublecopy(data, i_testArray.data, resolution[0]*resolution[1]*2);
	}


	inline
	void par_doublecopy(double *o_data, double *i_data, std::size_t i_size)
	{
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < i_size; i++)
			o_data[i] = i_data[i];
	}

	~Complex2DArrayFFT()
	{
		fftw_shutdown();

		cleanup();
	}



public:
	Complex2DArrayFFT& operator=(
			const Complex2DArrayFFT &i_testArray
	)
	{
		if (	data == nullptr ||
				resolution[0] != i_testArray.resolution[0] ||
				resolution[1] != i_testArray.resolution[1]
		)
		{
			setup(i_testArray.resolution, aliased_scaled);
		}

		aliased_scaled = i_testArray.aliased_scaled;
		resolution[0] = i_testArray.resolution[0];
		resolution[1] = i_testArray.resolution[1];

		assert(resolution[0] == i_testArray.resolution[0]);
		assert(resolution[1] == i_testArray.resolution[1]);

		par_doublecopy(data, i_testArray.data, resolution[0]*resolution[1]*2);
		return *this;
	}



	Complex2DArrayFFT toSpec()	const
	{
		Complex2DArrayFFT o_testArray(resolution, aliased_scaled);

		if (aliased_scaled)
		{
			assert(resolution[0] == fft_getSingleton_Plans().resolution_aliasing[0]);
			assert(resolution[1] == fft_getSingleton_Plans().resolution_aliasing[1]);

			fftw_execute_dft(
					fft_getSingleton_Plans().to_spec_aliasing,
					(fftw_complex*)this->data,
					(fftw_complex*)o_testArray.data
				);
		}
		else
		{
			assert(resolution[0] == fft_getSingleton_Plans().resolution[0]);
			assert(resolution[1] == fft_getSingleton_Plans().resolution[1]);

			fftw_execute_dft(
					fft_getSingleton_Plans().to_spec,
					(fftw_complex*)this->data,
					(fftw_complex*)o_testArray.data
				);
		}

		return o_testArray;
	}


	Complex2DArrayFFT toCart()
	{
		Complex2DArrayFFT o_testArray(resolution, aliased_scaled);

		if (aliased_scaled)
		{
			fftw_execute_dft(
					fft_getSingleton_Plans().to_cart_aliasing,
					(fftw_complex*)this->data,
					(fftw_complex*)o_testArray.data
				);
		}
		else
		{
			fftw_execute_dft(
					fft_getSingleton_Plans().to_cart,
					(fftw_complex*)this->data,
					(fftw_complex*)o_testArray.data
				);
		}

		/*
		 * do the scaling only if we convert the data back to cartesian space
		 */
		double scale = (1.0/((double)resolution[0]*(double)resolution[1]));
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			o_testArray.data[i] *= scale;
			o_testArray.data[i+1] *= scale;
		}

		return o_testArray;
	}


	inline
	void set(int y, int x, double re, double im)
	{
		data[(y*resolution[0]+x)*2+0] = re;
		data[(y*resolution[0]+x)*2+1] = im;
	}


	inline
	void set(int y, int x, const std::complex<double> &i_value)
	{
		data[(y*resolution[0]+x)*2+0] = i_value.real();
		data[(y*resolution[0]+x)*2+1] = i_value.imag();
	}


	inline
	void setRe(int y, int x, double re)
	{
		data[(y*resolution[0]+x)*2+0] = re;
	}

	inline
	void setIm(int y, int x, double im)
	{
		data[(y*resolution[0]+x)*2+1] = im;
	}

	inline
	double getRe(int y, int x)	const
	{
		return data[(y*resolution[0]+x)*2+0];
	}

	inline
	double getIm(int y, int x)	const
	{
		return data[(y*resolution[0]+x)*2+1];
	}

	inline
	complex get(int y, int x)	const
	{
		std::size_t idx = (y*resolution[0]+x)*2;

		return complex(
				data[idx+0],
				data[idx+1]
				);
	}

	void setAll(double re, double im)
	{
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			data[i] = re;
			data[i+1] = im;
		}
	}


	void zero()
	{
#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			data[i] = 0;
			data[i+1] = 0;
		}
	}


	void setAll(
			const std::complex<double> &i_value
	)
	{
		double re = i_value.real();
		double im = i_value.imag();

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			data[i] = re;
			data[i+1] = im;
		}
	}

	void setAllRe(double re)
	{
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			data[i] = re;
		}
	}

	void setAllIm(double im)
	{
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			data[i+1] = im;
		}
	}

#if 0
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
#endif



	/**
	 * apply a 3x3 stencil
	 */
	Complex2DArrayFFT op_stencil_3x3(
			const complex *i_kernel_data
	)
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

		int res_x = resolution[0];
		int res_y = resolution[1];

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y, out, i_kernel_data)
#endif
		for (int y = 0; y < res_y; y++)
		{
			for (int x = 0; x < res_x; x++)
			{
				double *data_out = &out.data[(y*res_x+x)*2];
				data_out[0] = 0;
				data_out[1] = 0;

				for (int j = -1; j <= 1; j++)
				{
					int pos_y = y+j;

					pos_y -= (pos_y >= res_y ? res_y : 0);
					pos_y += (pos_y < 0 ? res_y : 0);

					assert(pos_y >= 0 && pos_y < res_y);

					for (int i = -1; i <= 1; i++)
					{
						int pos_x = x+i;

						pos_x -= (pos_x >= res_x ? res_x : 0);
						pos_x += (pos_x < 0 ? res_x : 0);

						assert(pos_x >= 0 && pos_x < res_x);

						int idx = (j+1)*3+(i+1);
						assert(idx >= 0 && idx < 9);

						double kre = i_kernel_data[idx].real();
						double kim = i_kernel_data[idx].imag();

						double dre = data[(pos_y*res_x+pos_x)*2];
						double dim = data[(pos_y*res_x+pos_x)*2+1];

						data_out[0] += kre*dre - kim*dim;
						data_out[1] += kre*dim + kim*dre;
					}
				}
			}
		}

		return out;
	}


	/**
	 * apply a 3x3 stencil
	 */
	Complex2DArrayFFT op_stencil_Re_3x3(
			const double *i_kernel_data
	)
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

		int res_x = resolution[0];
		int res_y = resolution[1];

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y, out, i_kernel_data)
#endif
		for (int y = 0; y < res_y; y++)
		{
			for (int x = 0; x < res_x; x++)
			{
				double *data_out = &out.data[(y*res_x+x)*2];
				data_out[0] = 0;
				data_out[1] = 0;

				for (int j = -1; j <= 1; j++)
				{
					int pos_y = y+j;

					pos_y -= (pos_y >= res_y ? res_y : 0);
					pos_y += (pos_y < 0 ? res_y : 0);

					assert(pos_y >= 0 && pos_y < res_y);

					for (int i = -1; i <= 1; i++)
					{
						int pos_x = x+i;

						pos_x -= (pos_x >= res_x ? res_x : 0);
						pos_x += (pos_x < 0 ? res_x : 0);

						assert(pos_x >= 0 && pos_x < res_x);

						int idx = (j+1)*3+(i+1);
						assert(idx >= 0 && idx < 9);

						double kre = i_kernel_data[idx];

						double dre = data[(pos_y*res_x+pos_x)*2];
						double dim = data[(pos_y*res_x+pos_x)*2+1];

						data_out[0] += dre*kre;
						data_out[1] += dim*kre;
					}
				}
			}
		}

		return out;
	}



	/**
	 * apply a cross-shaped stencil by value real value 'A' given
	 *
	 * 0  Dy 0
	 * Dx 0  Dx
	 * 0  Dy 0
	 */
	Complex2DArrayFFT op_stencil_Re_X(
			const double i_scalar_Re_Dx,
			const double i_scalar_Re_Dy
	)
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

		int res_x = resolution[0];
		int res_y = resolution[1];


		// from right
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[y*res_x*2];
			double *i_data = &data[(y*res_x+1)*2];

			for (int x = 0; x < res_x-1; x++)
			{
				o_data[0] = i_data[0]*i_scalar_Re_Dx;
				o_data[1] = i_data[1]*i_scalar_Re_Dx;

				o_data += 2;
				i_data += 2;
			}

			i_data -= res_x*2;
			o_data[0] = i_data[0]*i_scalar_Re_Dx;
			o_data[1] = i_data[1]*i_scalar_Re_Dx;
		}


		// from left
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[(y*res_x+1)*2];
			double *i_data = &data[y*res_x*2];

			for (int x = 1; x < res_x; x++)
			{
				o_data[0] += i_data[0]*i_scalar_Re_Dx;
				o_data[1] += i_data[1]*i_scalar_Re_Dx;

				o_data += 2;
				i_data += 2;
			}

			o_data -= res_x*2;
			o_data[0] += i_data[0]*i_scalar_Re_Dx;
			o_data[1] += i_data[1]*i_scalar_Re_Dx;
		}


		// from upper
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[y*res_x*2];
			double *i_data = &data[(y == res_y-1 ? 0 : y+1)*res_x*2];

			for (int x = 0; x < res_x; x++)
			{
				o_data[0] += i_data[0]*i_scalar_Re_Dy;
				o_data[1] += i_data[1]*i_scalar_Re_Dy;

				o_data += 2;
				i_data += 2;
			}
		}


		// from lower
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[y*res_x*2];
			double *i_data = &data[(y == 0 ? res_y-1 : y-1)*res_x*2];

			for (int x = 0; x < res_x; x++)
			{
				o_data[0] += i_data[0]*i_scalar_Re_Dy;
				o_data[1] += i_data[1]*i_scalar_Re_Dy;

				o_data += 2;
				i_data += 2;
			}
		}

		return out;
	}

	/**
	 * apply a cross-shaped stencil by value real value 'A' given
	 *
	 * 0  Dy 0
	 * Dx C  Dx
	 * 0  Dy 0
	 */
	Complex2DArrayFFT op_stencil_Re_X_C(
			const double i_scalar_Re_Dx,
			const double i_scalar_Re_Dy,
			const double i_scalar_Re_C
	)
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

		int res_x = resolution[0];
		int res_y = resolution[1];

		// from right
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[y*res_x*2];
			double *i_data = &data[(y*res_x+1)*2];

			for (int x = 0; x < res_x-1; x++)
			{
				o_data[0] = i_data[0]*i_scalar_Re_Dx;
				o_data[1] = i_data[1]*i_scalar_Re_Dx;

				o_data += 2;
				i_data += 2;
			}

			i_data -= res_x*2;
			o_data[0] = i_data[0]*i_scalar_Re_Dx;
			o_data[1] = i_data[1]*i_scalar_Re_Dx;
		}


		// from left
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[(y*res_x+1)*2];
			double *i_data = &data[y*res_x*2];

			for (int x = 1; x < res_x; x++)
			{
				o_data[0] += i_data[0]*i_scalar_Re_Dx;
				o_data[1] += i_data[1]*i_scalar_Re_Dx;

				o_data += 2;
				i_data += 2;
			}

			o_data -= res_x*2;
			o_data[0] += i_data[0]*i_scalar_Re_Dx;
			o_data[1] += i_data[1]*i_scalar_Re_Dx;
		}


		// from upper
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[y*res_x*2];
			double *i_data = &data[(y == res_y-1 ? 0 : y+1)*res_x*2];

			for (int x = 0; x < res_x; x++)
			{
				o_data[0] += i_data[0]*i_scalar_Re_Dy;
				o_data[1] += i_data[1]*i_scalar_Re_Dy;

				o_data += 2;
				i_data += 2;
			}
		}


		// from lower
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[y*res_x*2];
			double *i_data = &data[(y == 0 ? res_y-1 : y-1)*res_x*2];

			for (int x = 0; x < res_x; x++)
			{
				o_data[0] += i_data[0]*i_scalar_Re_Dy;
				o_data[1] += i_data[1]*i_scalar_Re_Dy;

				o_data += 2;
				i_data += 2;
			}
		}


		// center
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp parallel for OPENMP_PAR_SIMD shared(res_x, res_y)
#endif
		for (int y = 0; y < res_y; y++)
		{
			double *o_data = &out.data[y*res_x*2];
			double *i_data = &data[(y*res_x)*2];

			for (int x = 0; x < res_x; x++)
			{
				o_data[0] += i_data[0]*i_scalar_Re_C;
				o_data[1] += i_data[1]*i_scalar_Re_C;

				o_data += 2;
				i_data += 2;
			}
		}

		return out;
	}


	/**
	 * setup the data in this class as a differential operator (d/dx)
	 */
public:
	void op_setup_diff_x(
			const double i_domain_size[2],
			bool i_use_complex_spectral_diffs = true
	)
	{
		setAll(0,0);

		if (i_use_complex_spectral_diffs)
		{
			double scale = 2.0*M_PIl/i_domain_size[0];

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < resolution[1]/2; j++)
			{
				for (std::size_t i = 1; i < resolution[0]/2; i++)
				{
					set(	j,
							i,
							0,
							(double)i*scale
						);
					set(
							resolution[1]-1-j,
							i,
							0,
							(double)i*scale
						);

					set(	j,
							resolution[0]-i,
							0,
							-(double)i*scale
						);
					set(
							resolution[1]-1-j,
							resolution[0]-i,
							0,
							-(double)i*scale
						);
				}
			}
		}
		else
		{
			double h[2] = {(double)i_domain_size[0] / (double)resolution[0], (double)i_domain_size[1] / (double)resolution[1]};

			/*
			 * setup FD operator
			 */
			set(0, 1, -1.0/(2.0*h[0]), 0);
			set(0, resolution[0]-1, 1.0/(2.0*h[0]), 0);

			*this = this->toSpec();
			// TODO: maybe set highest modes to zero?
		}
	}



public:
	void op_setup_diff2_x(
			const double i_domain_size[2],
			bool i_use_complex_spectral_diffs = true
	)
	{
		setAll(0,0);

		if (i_use_complex_spectral_diffs)
		{

			double scale = 2.0*M_PIl/i_domain_size[0];

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < resolution[1]/2; j++)
			{
				for (std::size_t i = 1; i < resolution[0]/2; i++)
				{
					set(	j,
							i,
							-(double)i*scale*(double)i*scale,
							0
						);
					set(
							resolution[1]-1-j,
							i,
							-(double)i*scale*(double)i*scale,
							0
						);

					set(	j,
							resolution[0]-i,
							-(double)i*scale*(double)i*scale,
							0
						);
					set(
							resolution[1]-1-j,
							resolution[0]-i,
							-(double)i*scale*(double)i*scale,
							0
						);
				}
			}
		}
		else
		{

			double h[2] = {
					(double)i_domain_size[0] / (double)resolution[0],
					(double)i_domain_size[1] / (double)resolution[1]
			};

			/*
			 * setup FD operator
			 */
			set(0, 1,
					1.0/(h[0]*h[0]), 0);
			set(0, 0,
					-2.0/(h[0]*h[0]), 0);
			set(0, resolution[0]-1,
					1.0/(h[0]*h[0]), 0);

			*this = this->toSpec();
			// TODO: maybe set highest modes to zero? Shouldn't be a big problem

		}
	}



public:
	void op_setup_diff_y(
			const double i_domain_size[2],
			bool i_use_complex_spectral_diffs = true
	)
	{
		setAll(0,0);

		if (i_use_complex_spectral_diffs)
		{

			double scale = 2.0*M_PIl/i_domain_size[1];

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (int j = 1; j < (int)resolution[1]/2; j++)
			{
				for (int i = 0; i < (int)resolution[0]/2; i++)
				{
					set(
							j,
							i,
							0,
							(double)j*scale
						);
					set(
							resolution[1]-j,
							i,
							0,
							-(double)j*scale
						);

					set(
							j,
							resolution[0]-i-1,
							0,
							(double)j*scale
						);
					set(
							resolution[1]-j,
							resolution[0]-i-1,
							0,
							-(double)j*scale
						);
				}
			}
		}
		else
		{
			double h[2] = {(double)i_domain_size[0] / (double)resolution[0], (double)i_domain_size[1] / (double)resolution[1]};

			/*
			 * setup FD operator
			 */
			set(1, 0, -1.0/(2.0*h[1]), 0);
			set(resolution[1]-1, 0, 1.0/(2.0*h[1]), 0);

			*this = this->toSpec();
			// TODO: maybe set highest modes to zero?

		}
	}


public:
	void op_setup_diff2_y(
			const double i_domain_size[2],
			bool i_use_complex_spectral_diffs = true
	)
	{
		setAll(0,0);


		if (i_use_complex_spectral_diffs)
		{

			double scale = 2.0*M_PIl/i_domain_size[1];

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (int j = 1; j < (int)resolution[1]/2; j++)
			{
				for (int i = 0; i < (int)resolution[0]/2; i++)
				{
					set(
							j,
							i,
							-(double)j*scale*(double)j*scale,
							0
						);
					set(
							resolution[1]-j,
							i,
							-(double)j*scale*(double)j*scale,
							0
						);

					set(
							j,
							resolution[0]-i-1,
							-(double)j*scale*(double)j*scale,
							0
						);
					set(
							resolution[1]-j,
							resolution[0]-i-1,
							-(double)j*scale*(double)j*scale,
							0
						);
				}
			}
		}
		else
		{

			double h[2] = {(double)i_domain_size[0] / (double)resolution[0], (double)i_domain_size[1] / (double)resolution[1]};

			/*
			 * setup FD operator
			 */
			set(1, 0, 1.0/(h[1]*h[1]), 0);
			set(0, 0, -2.0/(h[1]*h[1]), 0);
			set(resolution[1]-1, 0, 1.0/(h[1]*h[1]), 0);

			*this = this->toSpec();
			// TODO: maybe set highest modes to zero?

		}
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
#if 0
				if (std::abs(value_re) < 10e-10)
					value_re = 0;
				if (std::abs(value_im) < 10e-10)
					value_im = 0;
#endif
				o_ostream << "(" << value_re << ", " << value_im << ")\t";
			}
			o_ostream << std::endl;
		}
		return o_ostream;
	}



	inline
	Complex2DArrayFFT spec_div_element_wise(
			const Complex2DArrayFFT &i_d
//			double i_instability_threshold = 0
	)	const
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double ar = data[i];
			double ai = data[i+1];
			double br = i_d.data[i];
			double bi = i_d.data[i+1];

			double den = (br*br+bi*bi);

			//if (std::abs(den) == 0)
			// special handling merged with bit trick
			// For Laplace solution, this is the integration constant C
			out.data[i] = (std::abs(den) == 0 ? 0 : (ar*br + ai*bi)/den);
			out.data[i+1] = (std::abs(den) == 0 ? 0 : (ai*br - ar*bi)/den);
		}

		return out;
	}



	inline
	void spec_div_element_wise(
			const Complex2DArrayFFT &i_d,
			const Complex2DArrayFFT &o_out
//			double i_instability_threshold = 0
	)	const
	{
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double ar = data[i];
			double ai = data[i+1];
			double br = i_d.data[i];
			double bi = i_d.data[i+1];

			double den = (br*br+bi*bi);

			//if (std::abs(den) == 0)
			// special handling merged with bit trick
			// For Laplace solution, this is the integration constant C
			o_out.data[i] = (std::abs(den) == 0 ? 0 : (ar*br + ai*bi)/den);
			o_out.data[i+1] = (std::abs(den) == 0 ? 0 : (ai*br - ar*bi)/den);
		}
	}



	DataArray<2> getRealWithDataArray()
	{
		DataArray<2> out(resolution);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for
#endif

		for (std::size_t j = 0; j < resolution[1]; j++)
		{
#pragma omp OPENMP_SIMD
			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				out.array_data_cartesian_space[
									(j-out.range_start[1])*out.range_size[0]+
									(i-out.range_start[0])
								]
								= data[(j*resolution[0]+i)*2+0];
			}
		}

#if SWEET_USE_SPECTRAL_SPACE
		out.array_data_cartesian_space_valid = true;
		out.array_data_spectral_space_valid = false;
#endif

		return out;
	}
#if 0
	DataArray<2> getImagWithDataArray()
	{
		DataArray<2> out(resolution);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution[1]; j++)
		{
			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				out.set(j, i, getIm(j, i));
			}
		}


		return out;
	}
#endif


#if 1
	Complex2DArrayFFT &loadRealFromDataArray(
			const DataArray<2> &i_dataArray_Real
	)
	{
		i_dataArray_Real.requestDataInCartesianSpace();

// TODO: Make this SIMD
#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for //OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution[1]; j++)
		{
			#pragma omp OPENMP_SIMD
			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				data[(j*resolution[0]+i)*2+0] =
						i_dataArray_Real.array_data_cartesian_space[
								(j-i_dataArray_Real.range_start[1])*i_dataArray_Real.range_size[0]+
								(i-i_dataArray_Real.range_start[0])
							];

				data[(j*resolution[0]+i)*2+1] = 0;
			}
		}

		return *this;
	}
#endif

#if 0
	Complex2DArrayFFT &loadRealAndImagFromDataArrays(
			const DataArray<2> &i_dataArray_Real,
			const DataArray<2> &i_dataArray_Imag
	)
	{
		i_dataArray_Real.requestDataInCartesianSpace();
		i_dataArray_Imag.requestDataInCartesianSpace();

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution[1]; j++)
		{
			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				set(	j, i,
						i_dataArray_Real.get(j, i),
						i_dataArray_Imag.get(j, i)
				);
			}
		}

		return *this;
	}
#endif


#if 0
	DataArray<2> toDataArrays_Real()	const
	{
		DataArray<2> out(resolution);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution[1]; j++)
		{
			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				out.set(
					j, i,
					getRe(j, i)
				);
			}
		}

		return out;
	}
#endif

	void toDataArrays_Real(
			DataArray<2> &o_out
	)	const
	{

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution[1]; j++)
		{
			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				o_out.array_data_cartesian_space[
									(j-o_out.range_start[1])*o_out.range_size[0]+
									(i-o_out.range_start[0])
								] =
					data[(j*resolution[0]+i)*2+0];
			}
		}

#if SWEET_USE_SPECTRAL_SPACE
		o_out.array_data_cartesian_space_valid = true;
		o_out.array_data_spectral_space_valid = false;
#endif
	}


#if 0
	DataArray<2> toDataArrays_Imag()	const
	{
		DataArray<2> out(resolution);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution[1]; j++)
		{
			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				out.set(
					j, i,
					getIm(j, i)
				);
			}
		}

		return out;
	}
#endif


	/**
	 * Apply a linear operator given by this class to the input data array.
	 */
	inline
	Complex2DArrayFFT operator()(
			const Complex2DArrayFFT &i_array_data
	)	const
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double ar = data[i];
			double ai = data[i+1];
			double br = i_array_data.data[i];
			double bi = i_array_data.data[i+1];

			out.data[i] = ar*br - ai*bi;
			out.data[i+1] = ar*bi + ai*br;
		}

		return out;
	}


	/**
	 * Compute multiplication with a complex scalar
	 */
	inline
	Complex2DArrayFFT operator*(
			const std::complex<double> &i_value
	)	const
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

		double br = i_value.real();
		double bi = i_value.imag();

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
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
	 * Compute multiplication with real and imaginary components separately
	 */
	inline
	Complex2DArrayFFT multiply_real_imag(
			double i_real,
			double i_imag
	)	const
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i] = data[i] * i_real;
			out.data[i+1] = data[i+1] * i_imag;
		}
		return out;
	}



	/**
	 * Compute multiplication with a real value
	 */
	inline
	Complex2DArrayFFT operator*(
			double &i_value_real
	)	const
	{
		Complex2DArrayFFT out(resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double ar = data[i];
			double ai = data[i+1];

			out.data[i] = ar*i_value_real;
			out.data[i+1] = ai*i_value_real;
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
		Complex2DArrayFFT out(resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i] = data[i] + i_array_data.data[i];
			out.data[i+1] = data[i+1] + i_array_data.data[i+1];
		}

		return out;
	}




	/**
	 * Compute element-wise addition
	 */
	inline
	Complex2DArrayFFT& operator+=(
			const Complex2DArrayFFT &i_array_data
	)
	{
#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			data[i] += i_array_data.data[i];
			data[i+1] += i_array_data.data[i+1];
		}

		return *this;
	}



	/**
	 * Compute element-wise subtraction
	 */
	inline
	Complex2DArrayFFT operator-(
			const Complex2DArrayFFT &i_array_data
	)	const
	{
		Complex2DArrayFFT out(this->resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i] = data[i] - i_array_data.data[i];
			out.data[i+1] = data[i+1] - i_array_data.data[i+1];
		}

		return out;
	}



	/**
	 * Compute addition with complex value
	 */
	inline
	Complex2DArrayFFT addScalar_Spec(
			const complex &i_value
	)	const
	{
		Complex2DArrayFFT out(this->resolution, aliased_scaled);


#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i] = data[i];
			out.data[i+1] = data[i+1];
		}

		double scale = resolution[0]*resolution[1];
//		double scale = 1.0;
		out.data[0] += i_value.real()*scale;
		out.data[1] += i_value.imag()*scale;

		return out;
	}



	/**
	 * Compute addition with complex value
	 */
	inline
	Complex2DArrayFFT addScalar_Cart(
			const complex &i_value
	)	const
	{
		double re = i_value.real();
		double im = i_value.imag();

		Complex2DArrayFFT out(this->resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i] = data[i]+re;
			out.data[i+1] = data[i+1]+im;
		}

		return out;
	}



	/**
	 * Compute subtraction with complex value
	 */
	inline
	Complex2DArrayFFT subScalar_Spec(
			const complex &i_value
	)	const
	{
		Complex2DArrayFFT out(this->resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i] = data[i];
			out.data[i+1] = data[i+1];
		}

		double scale = resolution[0]*resolution[1];
		out.data[0] -= i_value.real()*scale;
		out.data[1] -= i_value.imag()*scale;

		return out;
	}



	/**
	 * Compute addition with complex value
	 */
	inline
	Complex2DArrayFFT subScalar_Cart(
			const complex &i_value
	)	const
	{
		Complex2DArrayFFT out(this->resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i] = data[i]-i_value.real();
			out.data[i+1] = data[i+1]-i_value.imag();
		}

		return out;
	}



	/**
	 * scale down
	 */
	inline
	Complex2DArrayFFT scale_down()	const
	{
		std::size_t res[2];
		res[0] = resolution[0]>>1;
		res[1] = resolution[1]>>1;

		// check for power of 2
		assert(res[0]*2 == resolution[0]);
		assert(res[1]*2 == resolution[1]);

		Complex2DArrayFFT out(res, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < res[1]; j++)
		{
			double *o_data = &out.data[j*res[0]*2];

			// first line
			double *i_data1 = &data[j*2*resolution[0]*2];
			// second line
			double *i_data2 = &data[(j*2+1)*resolution[0]*2];

			for (std::size_t i = 0; i < res[0]; i++)
			{
				// real
				o_data[0] = (1.0/4.0)*(i_data1[0] + i_data1[2] + i_data2[0] + i_data2[2]);

				// imag
				o_data[1] = (1.0/4.0)*(i_data1[1] + i_data1[3] + i_data2[1] + i_data2[3]);

				o_data += 2;
				i_data1 += 4;
				i_data2 += 4;
			}
		}

		return out;
	}



	/**
	 * scale up
	 */
	inline
	Complex2DArrayFFT scale_up()	const
	{
		std::size_t res[2];
		res[0] = resolution[0]*2;
		res[1] = resolution[1]*2;

		Complex2DArrayFFT out(res, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t j = 0; j < resolution[1]; j++)
		{
			// first line
			double *o_data1 = &out.data[(j*2)*res[0]*2];
			// second line
			double *o_data2 = &out.data[(j*2+1)*res[0]*2];

			double *i_data = &data[j*resolution[0]*2];

			for (std::size_t i = 0; i < resolution[0]; i++)
			{
				o_data1[0] = i_data[0];
				o_data1[1] = i_data[1];
				o_data1[2] = i_data[0];
				o_data1[3] = i_data[1];

				o_data2[0] = i_data[0];
				o_data2[1] = i_data[1];
				o_data2[2] = i_data[0];
				o_data2[3] = i_data[1];

				o_data1 += 4;
				o_data2 += 4;
				i_data += 2;
			}
		}

		return out;
	}



	/**
	 * Compute element-wise multiplication
	 */
	inline
	Complex2DArrayFFT operator*(
			const Complex2DArrayFFT &i_array_data
	)	const
	{
		Complex2DArrayFFT out(i_array_data.resolution, aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			out.data[i+0] = data[i]*i_array_data.data[i] - data[i+1]*i_array_data.data[i+1];
			out.data[i+1] = data[i]*i_array_data.data[i+1] + data[i+1]*i_array_data.data[i];
		}

		return out;
	}


	/**
	 * Compute dot product with complex conjugate
	 */
	inline
	std::complex<double> dotProd(
			const Complex2DArrayFFT &i_array_data
	)	const
	{
		Complex2DArrayFFT out(i_array_data.resolution, aliased_scaled);

		double sum_re = 0;
		double sum_im = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for OPENMP_PAR_SIMD reduction (+:sum_re,sum_im)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			sum_re += data[i]*i_array_data.data[i] - data[i+1]*i_array_data.data[i+1];
			sum_im += data[i]*i_array_data.data[i+1] + data[i+1]*i_array_data.data[i];
		}

		return std::complex<double>(sum_re, sum_im);
	}



	/**
	 * return the sum of the absolute values
	 */
	double reduce_sumAbs()	const
	{
		double sum = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
			sum += std::abs(data[i])+std::abs(data[i+1]);

		return sum;
	}



	/**
	 * return the sum
	 */
	std::complex<double> reduce_sum()	const
	{
		double sum_re = 0;
		double sum_im = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum_re,sum_im)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			sum_re += data[i+0];
			sum_im += data[i+1];
		}

		return std::complex<double>(sum_re, sum_im);
	}


#if 0
	/**
	 * return the sum
	 */
	double reduce_sum_real_imag()	const
	{
		double sum_re = 0;
		double sum_im = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum_re,sum_im)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			sum_re += data[i+0];
			sum_im += data[i+1];
		}

		return sum_re+sum_im;
	}
#else
	/**
	 * return the sum
	 */
	double reduce_sum_real_imag()	const
	{
		double sum = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
			sum += data[i+0] + data[i+1];

		return sum;
	}
#endif


	/**
	 * reduce to root mean square
	 */
	double reduce_rms_quad()
	{
		double sum = 0;
		double c = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double radius2 = data[i]*data[i]+data[i+1]*data[i+1];
			//double value = std::sqrt(radius2);
			//value *= value;
			double value = radius2;

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		sum = std::sqrt(sum/double(resolution[0]*resolution[1]));
		return sum;
	}



	/**
	 * reduce to root mean square
	 */
	double reduce_rms()
	{
		double sum = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
			sum += data[i]*data[i]+data[i+1]*data[i+1];

		sum = std::sqrt(sum/double(resolution[0]*resolution[1]));
		return sum;
	}




	/**
	 * reduce to max value
	 */
	double reduce_max()
	{
		double max_abs = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(max:max_abs)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
			max_abs = std::max(std::max(max_abs, std::abs(data[i])), std::abs(data[i+1]));

		return max_abs;
	}




	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double reduce_sum_re_quad()	const
	{
		double sum = 0;
		double c = 0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double value = data[i];

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
	double reduce_norm2_quad()	const
	{
		double sum = 0.0;
		double c = 0.0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum,c)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
		{
			double value = data[i]*data[i]+data[i+1]*data[i+1];

			// Use Kahan summation
			double y = value - c;
			double t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		sum -= c;

		return std::sqrt(sum);
	}

	/**
	 * return the maximum of all absolute values, use quad precision for reduction
	 */
	double reduce_norm2_squared()	const
	{
		double sum = 0.0;

#if !SWEET_REXI_THREAD_PARALLEL_SUM
		#pragma omp parallel for reduction(+:sum)
#endif
		for (std::size_t i = 0; i < resolution[0]*resolution[1]*2; i+=2)
			sum += data[i]*data[i]+data[i+1]*data[i+1];

		return sum;
	}

};


inline
static
Complex2DArrayFFT operator*(
		const std::complex<double> &i_value,
		const Complex2DArrayFFT &i_array_data
)
{
	Complex2DArrayFFT out(i_array_data.resolution, i_array_data.aliased_scaled);

	double br = i_value.real();
	double bi = i_value.imag();

#if !SWEET_REXI_THREAD_PARALLEL_SUM
	#pragma omp parallel for OPENMP_PAR_SIMD
#endif
	for (std::size_t i = 0; i < i_array_data.resolution[0]*i_array_data.resolution[1]*2; i+=2)
	{
		double ar = i_array_data.data[i];
		double ai = i_array_data.data[i+1];

		out.data[i] = ar*br - ai*bi;
		out.data[i+1] = ar*bi + ai*br;
	}
	return out;
}


inline
static
Complex2DArrayFFT operator*(
		double &i_value_re,
		const Complex2DArrayFFT &i_array_data
)
{
	Complex2DArrayFFT out(i_array_data.resolution, i_array_data.aliased_scaled);

#if !SWEET_REXI_THREAD_PARALLEL_SUM
	#pragma omp parallel for OPENMP_PAR_SIMD
#endif
	for (std::size_t i = 0; i < i_array_data.resolution[0]*i_array_data.resolution[1]*2; i+=2)
	{
		double ar = i_array_data.data[i];
		double ai = i_array_data.data[i+1];

		out.data[i] = ar*i_value_re;
		out.data[i+1] = ai*i_value_re;
	}
	return out;
}

#endif /* SRC_INCLUDE_SWEET_COMPLEX2DARRAYFFT_HPP_ */
