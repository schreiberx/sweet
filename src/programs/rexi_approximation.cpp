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
			"test-start",		/// start position
			"test-end",			/// end position
			"test-delta",		/// delta
			nullptr
	};

	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names, false))
	{
		std::cout << "User variables:" << std::endl;
		std::cout << std::endl;
		std::cout << "	--function-name=..." << std::endl;
		std::cout << "	--test-start=..." << std::endl;
		std::cout << "	--test-end=..." << std::endl;
		std::cout << "	--test-delta=..." << std::endl;
		return -1;
	}

	std::string function_name = simVars.bogus.var[0];
	double test_start = atof(simVars.bogus.var[1].c_str());
	double test_end = atof(simVars.bogus.var[2].c_str());
	double test_delta = atof(simVars.bogus.var[3].c_str());



	double max_error_threshold = 1e-8;

	simVars.rexi.outputConfig();

	std::vector<std::complex<double>> alpha;
	std::vector<std::complex<double>> beta;

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

	REXIFunctions<double> rexiFun;
	rexiFun.setup(function_name);

	TComplex I(0.0, 1.0);

	for (T x = test_start; x < test_end; x += test_delta)
	{
		TComplex approx = 0;
		for (std::size_t i = 0; i < alpha.size(); i++)
			approx += beta[i]/(TComplex(0, x) + alpha[i]);

		TComplex anal = rexiFun.eval(x*I);

		std::cout << x;
		std::cout << "\t" << approx.real();
		std::cout << "\t" << approx.imag();
		std::cout << "\t" << anal.real();
		std::cout << "\t" << anal.imag();
		std::cout << "\t" << std::abs(approx-anal);
		std::cout << std::endl;
	}

	return 0;
}
