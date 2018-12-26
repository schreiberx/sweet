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
			"function-name",		/// frequency multipliers for special scenario setup
			"lambda-real",		/// Real part of lambda
			"lambda-imag",		/// Imaginary part of lambda
			nullptr
	};

	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names, false))
	{
		std::cout << "User variables:" << std::endl;
		std::cout << std::endl;
		std::cout << "	--function-name=..." << std::endl;
		std::cout << "	--lambda-real=..." << std::endl;
		std::cout << "	--lambda-imag=..." << std::endl;
		return -1;
	}

	std::string function_name = simVars.bogus.var[0];
	double lambda_real = atof(simVars.bogus.var[1].c_str());
	double lambda_imag = atof(simVars.bogus.var[2].c_str());
	std::complex<double> lambda(lambda_real, lambda_imag);


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

	REXICoefficients<> rexiCoeffs = rexiCoefficientsSet.find_by_function_name(function_name);


	// Initial condition
	std::complex<double> U0(1.0, 0.0);

	std::complex<double> U = U0;

	for (	simVars.timecontrol.current_simulation_time = 0;
			simVars.timecontrol.current_simulation_time < simVars.timecontrol.max_simulation_time*(1.0-1e-12);
			simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size
	)
	{
		TComplex anal = rexiFunctions.eval(lambda*simVars.timecontrol.current_simulation_time);
		std::cout <<
				"t=" << simVars.timecontrol.current_simulation_time << "\t"
				"error=" << std::abs(anal-U)
		<< std::endl;

		std::complex<double> approx = rexiCoeffs.gamma*U;

		for (std::size_t i = 0; i < rexiCoeffs.alphas.size(); i++)
			approx += rexiCoeffs.betas[i]/(lambda*simVars.timecontrol.current_timestep_size + rexiCoeffs.alphas[i])*U;

		U = approx;

	}

	simVars.rexi.outputConfig();

	std::vector< std::complex<double> > alpha;
	std::vector< std::complex<double> > beta;

	if (simVars.rexi.rexi_method == "ci")
		if (simVars.timecontrol.current_timestep_size <= 0)
			FatalError("Please specify time step size with --dt=...");

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

	return 0;
}
