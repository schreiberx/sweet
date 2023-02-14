/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include <rexi/EXPFunctions.hpp>
#include <iostream>
#include <rexi/REXI.hpp>
#include <rexi/REXICoefficients.hpp>
#include <rexi/REXICoefficientsSet.hpp>
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
	const char *user_defined_prog_params[] = {
			"function-name",		/// Frequency multipliers for special scenario setup
			"lambda-real",			/// Real part of lambda
			"lambda-imag",			/// Imaginary part of lambda
			"test-mode",			/// Type of test
			nullptr
	};


	SimulationVariables simVars;

	simVars.user_defined.var[0] = "";
	simVars.user_defined.var[1] = "";
	simVars.user_defined.var[2] = "";
	simVars.user_defined.var[3] = "";

	if (!simVars.setupFromMainParameters(i_argc, i_argv, user_defined_prog_params, false))
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
	if (simVars.user_defined.var[0] != "")
		function_name = simVars.user_defined.var[0];

	double lambda_real = std::numeric_limits<double>::infinity();
	if (simVars.user_defined.var[1] != "")
		lambda_real = atof(simVars.user_defined.var[1].c_str());

	double lambda_imag = std::numeric_limits<double>::infinity();
	if (simVars.user_defined.var[2] != "")
		lambda_imag = atof(simVars.user_defined.var[2].c_str());

	std::complex<double> lambda(lambda_real, lambda_imag);

	int test_mode = 0;
	if (simVars.user_defined.var[3] != "")
		lambda_imag = atoi(simVars.user_defined.var[3].c_str());

	if (simVars.timecontrol.current_timestep_size <= 0)
	{
		std::cerr << "Error: Specify time step size" << std::endl;
		return -1;
	}

	if (std::isinf(std::abs(lambda)))
	{
		std::cerr << "Error: Specify \\lambda of linear operators" << std::endl;
		return -1;
	}


	/*
	 * Load analytical function
	 */
	EXPFunctions<double> rexiFunctions;
	rexiFunctions.setup(function_name);

	/*
	 * Load REXI coefficients from file
	 */
	REXICoefficientsSet<> rexiCoefficientsSet;

	if (simVars.rexi.exp_method == "direct")
	{
		SWEETError("Direct REXI mode not supported");
	}
	else if (simVars.rexi.exp_method == "file")
	{
		rexiCoefficientsSet.setup_from_files(simVars.rexi.rexi_files);

		if (rexiCoefficientsSet.rexiCoefficientVector.size() == 0)
			SWEETError("No REXI coefficient loaded");
	}
	else if (simVars.rexi.exp_method == "terry" || simVars.rexi.exp_method == "ci")
	{
		REXICoefficients<> rexiCoefficients;

		REXI<> rexi;
		rexi.load(&simVars.rexi, function_name, rexiCoefficients, 0);

		rexiCoefficientsSet.rexiCoefficientVector.push_back(rexiCoefficients);
	}
	else
	{
		SWEETError("This REXI method is not supported");
	}

	std::cout << "+ test_mode: " << test_mode << std::endl;

	simVars.rexi.outputConfig();

	for (std::size_t i = 0; i < rexiCoefficientsSet.rexiCoefficientVector.size(); i++)
	{
		const std::string &function_name = rexiCoefficientsSet.rexiCoefficientVector[i].function_name;

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



		for (	simVars.timecontrol.current_simulation_time = 0;
				simVars.timecontrol.current_simulation_time < simVars.timecontrol.max_simulation_time*(1.0-1e-12);
				simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size
		)
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

/*
	std::vector< std::complex<double> > alpha;
	std::vector< std::complex<double> > beta;

	if (simVars.rexi.exp_method == "ci")
		if (simVars.timecontrol.current_timestep_size <= 0)
			SWEETError("Please specify time step size with --dt=...");
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
