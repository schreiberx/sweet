/*
 * test_quadrature.cpp
 *
 *  Created on: 17 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include <libmath/GaussQuadrature.hpp>
#include <sweet/SWEETError.hpp>
#include <iostream>
#include <iomanip>


double fun(double x)
{
	return x*(x+1.0);
};

double int_fun(double x)
{
	return 1.0/3.0*x*x*x + 1.0/2.0*x*x;
};


int main(
		int i_argc,
		char *const i_argv[]
)
{
	double error_threshold = 1e-11;

	double int_start = 0;
	double int_end = 1;

	double int_value = int_fun(int_end)-int_fun(int_start);

	std::cout << std::setprecision(20);

	{
		double quad_val = GaussQuadrature::integrate5_intervals_adaptive_linear<double>(int_start, int_end, fun, error_threshold);
		std::cout << "integrate5_intervals_adaptive_linear - QuadVal: " << quad_val << ", IntValue: " << int_value << std::endl;

		if (std::abs(quad_val - int_value) > error_threshold)
			SWEETError("ERROR threshold exceeded");
	}


	{
		double quad_val = GaussQuadrature::integrate5_intervals_adaptive_recursive<double>(int_start, int_end, fun, error_threshold);
		std::cout << "integrate5_intervals_adaptive_recursive - QuadVal: " << quad_val << ", IntValue: " << int_value << std::endl;

		if (std::abs(quad_val - int_value) > error_threshold)
			SWEETError("ERROR threshold exceeded");
	}

	{
		double quad_val = GaussQuadrature::integrate5_intervals<double>(int_start, int_end, fun, 3);
		std::cout << "integrate5_intervals - QuadVal: " << quad_val << ", IntValue: " << int_value << std::endl;

		if (std::abs(quad_val - int_value) > error_threshold)
			SWEETError("ERROR threshold exceeded");
	}
}
