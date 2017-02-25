/*
 * DoubleAndQuadPrecisionStuff.hpp
 *
 *  Created on: 14 Dec 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_REXI_DQSTUFF_HPP_
#define SRC_INCLUDE_REXI_DQSTUFF_HPP_

#include <quadmath.h>
#include <ccomplex>


/**
 * Double and quad precision stuff
 */
class DQStuff
{
public:
	static
	std::complex<double> exp(
			const std::complex<double> &i_value
	)
	{
		std::complex<double> ret_val(
				std::cos(i_value.imag()),
				std::sin(i_value.imag())
			);

		ret_val = ret_val*std::exp(i_value.real());

		return ret_val;
	}

	static
	std::complex<float> exp(
			const std::complex<float> &i_value
	)
	{
		std::complex<float> ret_val(
				std::cos(i_value.imag()),
				std::sin(i_value.imag())
			);

		ret_val = ret_val*std::exp(i_value.real());

		return ret_val;
	}

	static
	std::complex<__float128> exp(
			const std::complex<__float128> &i_value
	)
	{
		std::complex<__float128> ret_val(
				::cosq(i_value.imag()),
				::sinq(i_value.imag())
			);

		return ret_val*::expq(i_value.real());
	}



public:
	static
	std::complex<double> expIm(
			const double &i_value
	)
	{
		return std::complex<double>(::cos(i_value), ::sin(i_value));
	}

	static
	std::complex<float> expIm(
			const float &i_value
	)
	{
		return std::complex<float>(::cos(i_value), ::sin(i_value));
	}

	static
	std::complex<__float128> expIm(
			const __float128 &i_value
	)
	{
		return std::complex<__float128>(::cosq(i_value), ::sinq(i_value));
	}



public:
	static
	double exp(
			const double &i_value
	)
	{
		return ::exp(i_value);
	}

	static
	__float128 exp(
			const __float128 &i_value
	)
	{
		return ::expq(i_value);
	}


public:
	static
	double sqrt(
			const double &i_value
	)
	{
		return ::sqrt(i_value);
	}

	static
	float sqrt(
			const float &i_value
	)
	{
		return ::sqrt(i_value);
	}

	static
	__float128 sqrt(
			const __float128 &i_value
	)
	{
		return ::sqrtq(i_value);
	}



public:
	static
	double abs(
			const double &i_value
	)
	{
		return ::fabs(i_value);
	}

	static
	float abs(
			const float &i_value
	)
	{
		return ::fabs(i_value);
	}

	static
	__float128 abs(
			const __float128 &i_value
	)
	{
		return ::fabsq(i_value);
	}



public:
	static
	double max(
			const double &i_value1,
			const double &i_value2
	)
	{
		return std::max(i_value1, i_value2);
	}

	static
	float max(
			const float &i_value1,
			const float &i_value2
	)
	{
		return std::max(i_value1, i_value2);
	}

	static
	__float128 max(
			const __float128 &i_value1,
			const __float128 &i_value2
	)
	{
		return ::fmaxq(i_value1, i_value2);
	}



public:
	static
	double pow(
			const double &i_value,
			const double &i_pi
	)
	{
		return ::pow(i_value, i_pi);
	}

	static
	float pow(
			const float &i_value,
			const float &i_pi
	)
	{
		return ::pow(i_value, i_pi);
	}

	static
	__float128 pow(
			const __float128 &i_value,
			const __float128 &i_pi
	)
	{
		return ::powq(i_value, i_pi);
	}



public:
	template <typename T>
	static
	std::complex<T> convertComplex(
			const std::complex<__float128> &i_value
	)
	{
		return std::complex<T>((T)i_value.real(), (T)i_value.imag());
	}


	template <typename T>
	static
	std::complex<T> convertComplex(
			const std::complex<double> &i_value
	)
	{
		return std::complex<T>((T)i_value.real(), (T)i_value.imag());
	}



public:
	template <typename T>
	static
	T fromString(
			const char* i_value
	)
	{
		return (T)::strtoflt128(i_value, nullptr);
	}

};



#endif /* SRC_INCLUDE_REXI_DQSTUFF_HPP_ */
