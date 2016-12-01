/*
 * test_rexi.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <rexi/GaussianApproximationNew.hpp>
#include <rexi/ExponentialApproximation.hpp>
#include <rexi/REXI.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>

int main(
		int i_argc,
		char *const i_argv[]
)
{
	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	double max_error_threshold = 1e-9;
	double max_error_threshold_machine = 1e-12;

	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "PHI 0 - EVALUATING GAUSSIAN APPROXIMATION (real)" << std::endl;
		std::cout << "******************************************************" << std::endl;
		GaussianApproximationNew ga;

		for (double h = 0.2; h > 0.001; h *= 0.5)
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

			if (std::abs(max_error) > max_error_threshold)
			{
				std::cerr << "MAX ERROR THRESHOLD EXCEEDED!" << std::endl;
				exit(-1);
			}
		}
	}

	return 0;
}
