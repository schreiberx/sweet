/*
 * GaussQuadrature.hpp
 *
 *  Created on: 17 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_GAUSSQUADRATURE_HPP_
#define SRC_INCLUDE_GAUSSQUADRATURE_HPP_

#include <cmath>
#include <limits>
#include <iostream>
#include <functional>
#include <sweet/libmath/DQStuff.hpp>
#include <sweet/core/SWEETError.hpp>


class GaussQuadrature
{
public:
	template <typename T = double>
	static T integrate5(
			T i_start,
			T i_end,
			std::function<T(T)> i_fun
	)
	{
		// see wikipedia
		// https://en.wikipedia.org/wiki/Gaussian_quadrature
		static T x[5] = {
				0,
				-(T)1.0/(T)3.0*DQStuff::sqrt((T)5.0-(T)2.0*DQStuff::sqrt((T)10.0/(T)7.0)),
				+(T)1.0/(T)3.0*DQStuff::sqrt((T)5.0-(T)2.0*DQStuff::sqrt((T)10.0/(T)7.0)),
				-(T)1.0/(T)3.0*DQStuff::sqrt((T)5.0+(T)2.0*DQStuff::sqrt((T)10.0/(T)7.0)),
				+(T)1.0/(T)3.0*DQStuff::sqrt((T)5.0+(T)2.0*DQStuff::sqrt((T)10.0/(T)7.0)),
		};

		static T w[5] = {
				(T)128.0/(T)225.0,
				((T)322.0+(T)13.0*DQStuff::sqrt((T)70.0))/(T)900.0,
				((T)322.0+(T)13.0*DQStuff::sqrt((T)70.0))/(T)900.0,
				((T)322.0-(T)13.0*DQStuff::sqrt((T)70.0))/(T)900.0,
				((T)322.0-(T)13.0*DQStuff::sqrt((T)70.0))/(T)900.0
		};

		T accum = 0;
		for (int i = 0; i < 5; i++)
		{
			accum += w[i]*i_fun((i_end-i_start)*(T)0.5*x[i]+ (i_end+i_start)*(T)0.5);
		}
		return accum*(i_end-i_start)*(T)0.5;
	}


	template <typename T = double>
	static T integrate5_intervals(
			T i_start,
			T i_end,
			std::function<T(T)> i_fun,
			int n_strips
	)
	{
		T sub_interval = (i_end - i_start)/(T)n_strips;

		T accum = 0;
		for (int i = 0; i < n_strips; i++)
			accum += integrate5(i_start + sub_interval*(T)i, i_start + sub_interval*(T)((T)i+(T)1.0), i_fun);

		return accum;
	}

	template <typename T = double>
	static T integrate5_intervals_adaptive_linear(
			T i_start,
			T i_end,
			std::function<T(T)> i_fun,
			T i_error_threshold = 10e-13
	)
	{
		T prev_value = std::numeric_limits<T>::infinity();

		int max_intervals = 1024*128;
		for (int num_intervals = 1; num_intervals < max_intervals; num_intervals++)
		{
			T value = integrate5_intervals(i_start, i_end, i_fun, num_intervals);
			T delta = std::abs(prev_value - value);

			if (delta <= i_error_threshold)
				return value;

			prev_value = value;
		}

		SWEETError("No convergence reached for integrate5_intervals_adaptive");
		return -1;
	}


	template <typename T = double>
	static T p_integrate5_intervals_adaptive_recursive(
			T i_start,
			T i_end,

			int i_current_depth,	///< current depth of recursion
			int i_min_depth,		///< minimum depth for recursion
			int i_max_depth,		///< maximum depth of recursion

			std::function<T(T)> i_fun,
			T i_rel_error_threshold = 10e-13,
			T i_prev_value = std::numeric_limits<T>::infinity()
	)
	{
		assert(i_end > i_start);

		T mid = (T)0.5*(i_end + i_start);

		if (i_current_depth > i_max_depth)
			return std::numeric_limits<T>::infinity();

		if (i_current_depth >= i_min_depth)
		{
			T left_value = integrate5(
					i_start,
					mid,
					i_fun
				);

			T right_value = integrate5(
					mid,
					i_end,
					i_fun
				);

			T sum = left_value + right_value;

			if (DQStuff::abs(sum-i_prev_value)/(i_end-i_start) < i_rel_error_threshold)
				return sum;


			T recursive_left_value = p_integrate5_intervals_adaptive_recursive(
					i_start,
					mid,
					i_current_depth+1,
					i_min_depth,
					i_max_depth,
					i_fun,
					i_rel_error_threshold,
					left_value
				);

			T recursive_right_value = p_integrate5_intervals_adaptive_recursive(
					mid,
					i_end,
					i_current_depth+1,
					i_min_depth,
					i_max_depth,
					i_fun,
					i_rel_error_threshold,
					right_value
				);

			return recursive_left_value + recursive_right_value;
		}
		else
		{

			T recursive_left_value = p_integrate5_intervals_adaptive_recursive(
					i_start,
					mid,
					i_current_depth+1,
					i_min_depth,
					i_max_depth,
					i_fun,
					i_rel_error_threshold
				);

			T recursive_right_value = p_integrate5_intervals_adaptive_recursive(
					mid,
					i_end,
					i_current_depth+1,
					i_min_depth,
					i_max_depth,
					i_fun,
					i_rel_error_threshold
				);

			return recursive_left_value + recursive_right_value;
		}
	}


	template <typename T = double>
	static T integrate5_intervals_adaptive_recursive(
			T i_start,
			T i_end,
			std::function<T(T)> i_fun,
			T i_rel_error_threshold = 10e-12
	)
	{
		T approx_integral = p_integrate5_intervals_adaptive_recursive(
				i_start,
				i_end,
				0,		/// current depth
				4,		/// min depth
				32,		/// max depth
				i_fun,
				i_rel_error_threshold
			);

		if (std::isinf((double)approx_integral))
			SWEETError("No convergence reached for integrate5_intervals_adaptive");

		return approx_integral;
	}
};



#endif
