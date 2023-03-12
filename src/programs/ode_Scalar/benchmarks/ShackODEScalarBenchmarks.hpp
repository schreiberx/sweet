/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_

#include <string>
#include <iostream>
#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/StringSplit.hpp>

/**
 * Values and parameters to setup benchmarks simulations
 */
class ShackODEScalarBenchmarks	:
		public sweet::ShackInterface
{
public:
	/// benchmark scenario
	std::string benchmark_name = "";

	/// May the benchmark setup overwrite the simulation variables
	bool benchmark_override_simvars = true;

	/**
	 * Parameters a and b for du/dt = a * f1(u, t) + b * f2(u, t)
	 */
	double ode_parameters[2] = {0, 0};

	/**
	 * Initial solution
	 */
	double u0 = 0;


	bool validateNonStationaryODE()
	{
		if (ode_parameters[0] == 0 && ode_parameters[1] == 0)
			return error.set("Both ODE parameters are 0, use --advection-velocity=...");

		return true;
	}

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << i_prefix << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << i_prefix << "	--ode-parameters=[float],[float]	ode parameters a and b" << std::endl;
		std::cout << i_prefix << "	--u0=[float]	initial solution for the ODE" << std::endl;
		std::cout << i_prefix << "	--benchmark-override-simvars [bool]	Allow overwriting simulation variables by benchmark (default: 1)" << std::endl;
		std::cout << i_prefix << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

		std::string tmp;
		if (i_pa.getArgumentValueByKey("--ode-parameters", tmp))
			StringSplit::split2double(tmp, &ode_parameters[0], &ode_parameters[1]);

		i_pa.getArgumentValueByKey("--u0", u0);

		i_pa.getArgumentValueByKey("--benchmark-override-simvars", benchmark_override_simvars);

		if (error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "BENCHMARK:" << std::endl;
		std::cout << i_prefix << " + benchmark_name: " << benchmark_name << std::endl;
		std::cout << i_prefix << " + ode_parameters (a, b): " << ode_parameters[0] << ", " << ode_parameters[1] << std::endl;
		std::cout << i_prefix << " + benchmark_override_simvars: " << benchmark_override_simvars << std::endl;
		std::cout << i_prefix << std::endl;
	}
};



#endif
