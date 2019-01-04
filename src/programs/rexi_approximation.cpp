/*
 * rexi_approximation.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <rexi/REXI.hpp>
#include <rexi/REXICoefficients.hpp>
#include <rexi/REXICoefficientsSet.hpp>
#include <rexi/REXIFunctions.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/Timeloop.hpp>
#include <stdlib.h>


typedef double T;
typedef std::complex<T> TComplex;

std::string function_name;

int main(
		int i_argc,
		char *const i_argv[]
)
{
	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"function-name",		/// Frequency multipliers for special scenario setup
			"lambda-real",			/// Real part of lambda
			"lambda-imag",			/// Imaginary part of lambda
			"test-type",			/// Type of test
			nullptr
	};


	SimulationVariables simVars;

	simVars.bogus.var[0] = "";
	simVars.bogus.var[1] = "";
	simVars.bogus.var[2] = "";
	simVars.bogus.var[3] = "";

	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names, false))
	{
		std::cout << "User variables:" << std::endl;
		std::cout << std::endl;
		std::cout << "	--function-name=..." << std::endl;
		std::cout << "	--lambda-real=..." << std::endl;
		std::cout << "	--lambda-imag=..." << std::endl;
		std::cout << "  --test-mode=..." << std::endl;
		std::cout << "        0: use standard time stepping" << std::endl;
		std::cout << "        1: always start from u(0) with increasing time step sizes" << std::endl;
		return -1;
	}

	std::string function_name;
	if (simVars.bogus.var[0] != "")
		function_name = simVars.bogus.var[0];

	double lambda_real = 0.0;
	if (simVars.bogus.var[1] != "")
		lambda_real = atof(simVars.bogus.var[1].c_str());

	double lambda_imag = 0.0;
	if (simVars.bogus.var[2] != "")
		lambda_imag = atof(simVars.bogus.var[2].c_str());

	std::complex<double> lambda(lambda_real, lambda_imag);

	int test_mode = 0;
	if (simVars.bogus.var[3] != "")
		lambda_imag = atoi(simVars.bogus.var[3].c_str());

	if (simVars.timecontrol.current_timestep_size <= 0)
	{
		std::cerr << "Error: Specify time step size" << std::endl;
		return -1;
	}

	if (std::abs(lambda) == 0)
	{
		std::cerr << "Error: Specify \\lambda of linear operators" << std::endl;
		return -1;
	}


	/*
	 * Load analytical function
	 */
	REXIFunctions<double> rexiFunctions;
	rexiFunctions.setup(function_name);

	/*
	 * Load REXI coefficients from file
	 */
	REXICoefficientsSet<> rexiCoefficientsSet;
	rexiCoefficientsSet.setup_from_files(simVars.rexi.rexi_files);


	std::cout << "+ test_mode: " << test_mode << std::endl;

	for (std::size_t i = 0; i < rexiCoefficientsSet.rexiCoefficients.size(); i++)
	{
		const std::string &function_name = rexiCoefficientsSet.rexiCoefficients[i].function_name;

		std::cout << "Running tests for function " << function_name << std::endl;


		// Load coefficients for function
		REXICoefficients<> rexiCoeffs = rexiCoefficientsSet.find_by_function_name(function_name);

		// Initial condition
		std::complex<double> U0(1.0, 0.0);

		// Current solution
		std::complex<double> U = U0;

		auto computeAndOutputError = [&](std::complex<double> &i_U) -> double
		{
			TComplex analU = rexiFunctions.eval(lambda*simVars.timecontrol.current_simulation_time);

			double error = std::abs(analU-U);
			std::cout <<
					"t=" << simVars.timecontrol.current_simulation_time << "\t"
					"error=" << error
					<< std::endl;

			return error;
		};


		SWEET_TIMELOOP
		{
			computeAndOutputError(U);

			// REXI time integration
			{
				std::complex<double> approx = rexiCoeffs.gamma*U;

				for (std::size_t i = 0; i < rexiCoeffs.alphas.size(); i++)
					approx += rexiCoeffs.betas[i]/(lambda*simVars.timecontrol.current_timestep_size + rexiCoeffs.alphas[i])*U;

				U = approx;
			}
		}

		double error = computeAndOutputError(U);

		std::cout << "[MULE] error: " << error << std::endl;
	}

	simVars.rexi.outputConfig();
/*
	std::vector< std::complex<double> > alpha;
	std::vector< std::complex<double> > beta;

	if (simVars.rexi.rexi_method == "ci")
		if (simVars.timecontrol.current_timestep_size <= 0)
			FatalError("Please specify time step size with --dt=...");
*/
/*
	std::cout << "Loading REXI coefficients..." << std::flush;
	REXI rexi;
	rexi.load(
			&simVars.rexi,
			function_name,
			alpha,
			beta,
			simVars.timecontrol.current_timestep_size,
			simVars.misc.verbosity
		);
	std::cout << "OK" << std::endl;
*/

	return 0;
}
