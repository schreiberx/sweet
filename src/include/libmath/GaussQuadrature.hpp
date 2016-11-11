/*
 * GaussQuadrature.hpp
 *
 *  Created on: 17 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_GAUSSQUADRATURE_HPP_
#define SRC_INCLUDE_GAUSSQUADRATURE_HPP_

#include <functional>
#include <cmath>
#include <limits>
#include <iostream>
#include <sweet/FatalError.hpp>


class GaussQuadrature
{
public:
	static double integrate5(
			double i_start,
			double i_end,
			std::function<double(double)> i_fun
	)
	{
		// see wikipedia
		// https://en.wikipedia.org/wiki/Gaussian_quadrature
		static double x[5] = {
				0,
				-1.0/3.0*std::sqrt(5.0-2.0*std::sqrt(10.0/7.0)),
				+1.0/3.0*std::sqrt(5.0-2.0*std::sqrt(10.0/7.0)),
				-1.0/3.0*std::sqrt(5.0+2.0*std::sqrt(10.0/7.0)),
				+1.0/3.0*std::sqrt(5.0+2.0*std::sqrt(10.0/7.0)),
		};

		static double w[5] = {
				128.0/225.0,
				(322.0+13.0*std::sqrt(70.0))/900.0,
				(322.0+13.0*std::sqrt(70.0))/900.0,
				(322.0-13.0*std::sqrt(70.0))/900.0,
				(322.0-13.0*std::sqrt(70.0))/900.0
		};

		double accum = 0;
		for (int i = 0; i < 5; i++)
		{
			accum += w[i]*i_fun((i_end-i_start)*0.5*x[i]+ (i_end+i_start)*0.5);
		}
		return accum*(i_end-i_start)*0.5;
	}


	static double integrate5_intervals(
			double i_start,
			double i_end,
			std::function<double(double)> i_fun,
			int n_strips
	)
	{
		double sub_interval = (i_end - i_start)/(double)n_strips;

		double accum = 0;
		for (int i = 0; i < n_strips; i++)
			accum += integrate5(i_start + sub_interval*(double)i, i_start + sub_interval*(double)(i+1), i_fun);

		return accum;
	}

	static double integrate5_intervals_adaptive_linear(
			double i_start,
			double i_end,
			std::function<double(double)> i_fun,
			double i_error_threshold = 10e-13
	)
	{
		double prev_value = std::numeric_limits<double>::infinity();

		int max_intervals = 1024*128;
		for (int num_intervals = 1; num_intervals < max_intervals; num_intervals++)
		{
			double value = integrate5_intervals(i_start, i_end, i_fun, num_intervals);
			double delta = std::abs(prev_value - value);

			if (delta <= i_error_threshold)
				return value;

			prev_value = value;
		}

		FatalError("No convergence reached for integrate5_intervals_adaptive");
		return -1;
	}


	static double p_integrate5_intervals_adaptive_recursive(
			double i_start,
			double i_end,

			int i_current_depth,	///< current depth of recursion
			int i_min_depth,		///< minimum depth for recursion
			int i_max_depth,		///< maximum depth of recursion

			std::function<double(double)> i_fun,
			double i_rel_error_threshold = 10e-13,
			double i_prev_value = std::numeric_limits<double>::infinity()
	)
	{
		assert(i_end > i_start);

		double mid = 0.5*(i_end + i_start);

		if (i_current_depth > i_max_depth)
			return std::numeric_limits<double>::infinity();

		if (i_current_depth >= i_min_depth)
		{
			double left_value = integrate5(
					i_start,
					mid,
					i_fun
				);

			double right_value = integrate5(
					mid,
					i_end,
					i_fun
				);

			double sum = left_value + right_value;

			if (std::abs(sum-i_prev_value)/(i_end-i_start) < i_rel_error_threshold)
				return sum;


			double recursive_left_value = p_integrate5_intervals_adaptive_recursive(
					i_start,
					mid,
					i_current_depth+1,
					i_min_depth,
					i_max_depth,
					i_fun,
					i_rel_error_threshold,
					left_value
				);

			double recursive_right_value = p_integrate5_intervals_adaptive_recursive(
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

			double recursive_left_value = p_integrate5_intervals_adaptive_recursive(
					i_start,
					mid,
					i_current_depth+1,
					i_min_depth,
					i_max_depth,
					i_fun,
					i_rel_error_threshold
				);

			double recursive_right_value = p_integrate5_intervals_adaptive_recursive(
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


	static double integrate5_intervals_adaptive_recursive(
			double i_start,
			double i_end,
			std::function<double(double)> i_fun,
			double i_rel_error_threshold = 10e-12
	)
	{
		double approx_integral = p_integrate5_intervals_adaptive_recursive(
				i_start,
				i_end,
				0,		/// current depth
				4,		/// min depth
				32,		/// max depth
				i_fun,
				i_rel_error_threshold
			);

		if (std::isinf(approx_integral))
			FatalError("No convergence reached for integrate5_intervals_adaptive");

		return approx_integral;
	}
};



#endif /* SRC_INCLUDE_GAUSSQUADRATURE_HPP_ */
