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

	static double integrate5_intervals_adaptive(
			double i_start,
			double i_end,
			std::function<double(double)> i_fun,
			double i_error_threshold = 10e-13
	)
	{
//		double delta = std::numeric_limits<double>::infinity();
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
};



#endif /* SRC_INCLUDE_GAUSSQUADRATURE_HPP_ */
