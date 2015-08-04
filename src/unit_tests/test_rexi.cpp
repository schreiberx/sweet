/*
 * rexi.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <rexi/GaussianApproximation.hpp>
#include <rexi/ExponentialApproximation.hpp>
#include <rexi/REXI.hpp>

int main(int i_argc, char *i_argv[])
{
	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING GAUSSIAN APPROXIMATION (real)" << std::endl;
		std::cout << "******************************************************" << std::endl;
		GaussianApproximation ga;

//		ga.print();

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


	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING EXPONENTIAL APPROXIMATION with approx Gaussian (approx_e_ix)" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (double h = 2; h > 0.001; h *= 0.5)
		{
			int M = 32/h;
			ExponentialApproximation ea(h, M);

			double start = -2.0;
			double end = 2.0;
			double step_size = 0.01;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error = std::abs(ea.eval_e_ix(x) - ea.approx_e_ix(x));
				max_error = std::max(max_error, error);
			}

			std::cout << "max_error: " << max_error << " for h " << h << " and M " << M << std::endl;
		}
	}


#if 0
	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING REAL EXPONENTIAL APPROXIMATION with approximated Gaussian (approx_e_ix_returnReal)" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (double h = 2; h > 0.01; h *= 0.5)
		{
			int M = 32/h;

			ExponentialApproximation ea(h, M);

			double start = -2.0;
			double end = 2.0;
			double step_size = 0.01;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error = std::abs(ea.eval_e_ix(x).real() - ea.approx_e_ix_returnReal(x));
				max_error = std::max(max_error, error);
			}

			std::cout << "max_error: " << max_error << " for h " << h << " and M " << M << std::endl;
		}
	}
#endif

	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "REXI complex" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (double h = 2; h > 0.01; h *= 0.5)
		{
			int M = 256/h;

			REXI rexi(h, M);

			double start = -2.0;
			double end = 2.0;
			double step_size = 0.01;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error = std::abs(rexi.eval_e_ix(x) - rexi.approx_e_ix(x));
				max_error = std::max(max_error, error);
			}

			std::cout << "max_error: " << max_error << " for h " << h << " and M " << M << std::endl;
		}
	}

	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "REXI real" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (double h = 2; h > 0.01; h *= 0.5)
		{
			int M = 256/h;

			REXI rexi(h, M);

			double start = -2.0;
			double end = 2.0;
			double step_size = 0.01;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error_real = std::abs(rexi.eval_e_ix(x).real() - rexi.approx_e_ix_returnReal(x));
				max_error = std::max(max_error, error_real);
			}

			std::cout << "max_error: " << max_error << " for h " << h << " and M " << M << std::endl;
		}
	}

	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "REXI imag" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (double h = 2; h > 0.01; h *= 0.5)
		{
			int M = 256/h;

			REXI rexi(h, M);

			double start = -2.0;
			double end = 2.0;
			double step_size = 0.01;

			double max_error = 0;

			for (double x = start; x < end; x += step_size)
			{
				double error_imag = std::abs(rexi.eval_e_ix(x).imag() - rexi.approx_e_ix_returnImag(x));
				max_error = std::max(max_error, error_imag);
			}

			std::cout << "max_error: " << max_error << " for h " << h << " and M " << M << std::endl;
		}
	}

	return 0;
}
