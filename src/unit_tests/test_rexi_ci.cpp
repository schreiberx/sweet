/*
 * test_rexi.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <rexi/REXI.hpp>
#include <rexi/REXIFunctions.hpp>
#include <sweet/SimulationVariables.hpp>



typedef double T;
typedef std::complex<T> TComplex;



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

	double max_error_threshold = 1e-8;


	for (int fun_id = 0; fun_id <= 5; fun_id++)
	{
		std::string function_name;
		switch(fun_id)
		{
		case 0:
			function_name = "phi0";
			break;

		case 1:
			function_name = "phi1";
			break;

		case 2:
			function_name = "phi2";
			break;

		case 3:
			function_name = "phi3";
			break;

		case 4:
			function_name = "phi4";
			break;

		case 5:
			function_name = "phi5";
			break;

		default:
			FatalError("this phi function is not implemented");
		}

		if (simVars.rexi.rexi_method != "ci")
			FatalError("This test is for rexi_method=='ci' only");

		if (simVars.rexi.use_half_poles)
			FatalError("Not yet supported");

		for (T r = 5.0; r < 25.0; r += 5.0)
		{
			{
				simVars.rexi.ci_n = (int)(32*r);
				simVars.rexi.ci_mu = 1.0;
				simVars.rexi.ci_r = r;

				std::cout << "******************************************************" << std::endl;
				std::cout << "PHI " << fun_id << " - REXI real: Test for partition of unity and accuracy" << std::endl;
				std::cout << "******************************************************" << std::endl;
				std::cout << "******************************************************" << std::endl;
				std::cout << "N: " << simVars.rexi.ci_n << std::endl;
				std::cout << "r: " << simVars.rexi.ci_r << std::endl;
				std::cout << "mu: " << simVars.rexi.ci_mu << std::endl;
				std::cout << "******************************************************" << std::endl;
				std::cout << "Setup coefficients... " << std::flush;

				std::vector<std::complex<double>> alpha;
				std::vector<std::complex<double>> beta;

				REXI rexi;
				rexi.load(&simVars.rexi, function_name, alpha, beta, simVars.misc.verbosity);
				std::cout << "OK" << std::endl;


				REXIFunctions<__float128> rexiFunctions(function_name);

				// REXI approximates the interval [-M*h;M*h] but gets inaccurate close to the interval boundaries
				T start = -r*0.7+abs(simVars.rexi.ci_mu);
				T end = -start;
				T step_size = 0.011;

				T max_error_real = 0.0;
				T max_error_imag = 0.0;

				double max_error_threshold_local = max_error_threshold*simVars.rexi.ci_n;

				for (T x = start; x < end; x += step_size)
				{
					std::complex<__float128> correct_ = rexiFunctions.eval(std::complex<__float128>(0, x));
					TComplex correct(correct_.real(), correct_.imag());

					TComplex approx = 0;
					for (std::size_t i = 0; i < alpha.size(); i++)
						approx += beta[i]/(TComplex(0, x) + alpha[i]);

#if 0
					std::cout << (double)x << ": " << correct << "\t" << approx << "\t" << (correct-approx) << std::endl;
					continue;
#endif

					if (DQStuff::abs(approx.real()) > 1.0)
						std::cerr << "approx value_real " << approx.real() << " not bounded by unity (just a warning and not a problem) at x=" << (double)x << std::endl;

					if (DQStuff::abs(approx.imag()) > 1.0)
						std::cerr << "approx value_imag " << approx.imag() << " not bounded by unity (just a warning and not a problem) at x=" << (double)x << std::endl;

					T error_real = DQStuff::abs(correct.real() - approx.real());
					T error_imag = DQStuff::abs(correct.imag() - approx.imag());

					max_error_real = DQStuff::max(max_error_real, error_real);
					max_error_imag = DQStuff::max(max_error_imag, error_imag);
				}

				std::cout << "max_error_real: " << (double)max_error_real << std::endl;

				if (max_error_real > max_error_threshold_local)
				{
					std::cout << "max_error_threshold: " << max_error_threshold_local << std::endl;
					FatalError("MAX ERROR THRESHOLD EXCEEDED for real values!");
				}

				std::cout << "max_error_imag: " << (double)max_error_imag << std::endl;

				if (max_error_imag > max_error_threshold_local)
				{
					std::cout << "max_error_threshold_imag: " << max_error_threshold_local << std::endl;
					FatalError("MAX ERROR THRESHOLD EXCEEDED for imaginary value!");
				}
			}
		}
	}


	return 0;
}
