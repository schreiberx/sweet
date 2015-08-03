/*
 * rexi.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <rexi/GaussianApproximation.hpp>
#include <rexi/ExponentialApproximation.hpp>

int main(int i_argc, char *i_argv[])
{
	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING GAUSSIAN APPROXIMATION" << std::endl;
		std::cout << "******************************************************" << std::endl;
		GaussianApproximation ga;

		for (double h = 8; h > 0.1; h *= 0.5)
		{
			double start = -100.0;
			double end = 100.0;
			double step_size = 0.001;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error = std::abs(ga.evalGaussian(x, h) - ga.approxGaussian(x, h));
				max_error = std::max(max_error, error);
			}

			std::cout << "max_error: " << max_error << " for h " << h << std::endl;
		}
	}


	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING EXPONENTIAL APPROXIMATION with analytical gaussian approximation" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (double h = 2; h > 0.001; h *= 0.5)
		{
			int M = 128/h;
			ExponentialApproximation ea(h, M);

			double start = -2.0;
			double end = 2.0;
			double step_size = 0.01;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error = std::abs(ea.evalExponential(x) - ea.approxExponentialEvalGaussian(x));
				max_error = std::max(max_error, error);
			}

			std::cout << "max_error: " << max_error << " for h " << h << " and M " << M << std::endl;
		}
	}

	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING EXPONENTIAL APPROXIMATION with approximated gaussian approximation" << std::endl;
		std::cout << "******************************************************" << std::endl;

	//	ea.print();

		for (double h = 2; h > 0.01; h *= 0.5)
		{
			int M = 128/h;

			ExponentialApproximation ea(h, M);

			double start = -2.0;
			double end = 2.0;
			double step_size = 0.01;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error = std::abs(ea.evalExponential(x) - ea.approxExponentialApproxGaussian(x));
				max_error = std::max(max_error, error);
			}

			std::cout << "max_error: " << max_error << " for h " << h << " and M " << M << std::endl;
		}
	}

	return 0;
}
