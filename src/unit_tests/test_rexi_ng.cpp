/*
 * test_rexi_ng.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <rexi/RexiFile.hpp>
#include <iostream>
#include <sweet/SimulationVariables.hpp>
//#include <rexi/REXI.hpp>


typedef double T;

int main(
		int i_argc,
		char *const i_argv[]
)
{
	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr, false))
	{
		return -1;
	}

	double max_error_threshold = 1e-9;


	for (int fun_id = 0; fun_id <= 2; fun_id++)
	{
		std::ostringstream os;
		os << "phi" << fun_id;
		std::string function_name = os.str();

		std::cout << "******************************************************" << std::endl;
		std::cout << function_name << " - REXI real: Test for partition of unity and errors" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (T h = 0.2; h >= 0.01; h *= 0.5)
		{
//			for (int M = 32; M < 1024; M *= 2)
			for (int M = 32; M < 1024; M *= 2)
			{
				double test_min = -h*(T)M;
				double test_max = h*(T)M;

				RexiFile<T> rexiNG;
				bool retval = rexiNG.auto_load(
						function_name,
						0,
						rexiNG.None(),	/// max_error
						max_error_threshold,			/// max_error_double_precision
						test_min,		/// test_min
						test_max,		/// test_max
						rexiNG.None(),	/// basis_function_scaling
						rexiNG.None(),	/// basis_function_spacing
						rexiNG.None(),	/// basis_function_rat_shift

						false,			/// reduce to half
						"data/faf_data"	/// FAF data directory
				);

				std::cout << "Requesting REXI [" << test_min << ", " << test_max << "]" << std::endl;
				rexiNG.fafcoeffs.output();

				if (!retval)
					FatalError("No suitable coefficients found");

				//rexiNG.fafcoeffs.output();

				// REXI approximates the interval [-M*h;M*h] but gets inaccurate close to the interval boundaries
				double start = test_min;
				double end = test_max;

				double step_size = 0.001;

				T max_error = 0;
				T max_error_x = 0;

				for (double x = start; x < end; x += step_size)
				{
					T correct = rexiNG.eval_real(x);
					T approx = rexiNG.approx_real(x);

					//std::cout << x << "\t" << correct << "\t" << approx << std::endl;

					if (DQStuff::abs(approx) > 1.0)
					{
						std::cerr << "approx value " << (double)approx << " not bounded by unity (just a warning and not a problem) at x=" << x << std::endl;
						std::cerr << "correct: " << correct << std::endl;
						std::cerr << "approx: " << approx << std::endl;
					}
					T error_real = DQStuff::abs(correct - approx);

					if (error_real > max_error)
					{
						max_error = error_real;
						max_error_x = x;
					}
				}

				std::cout << "max_error: " << (double)max_error << " for h=" << (double)h << " and M=" << M << " at x=" << max_error_x << std::endl;

				if (DQStuff::abs(max_error) > max_error_threshold)
					FatalError("MAX ERROR THRESHOLD EXCEEDED!");
			}
		}
	}


	return 0;
}
