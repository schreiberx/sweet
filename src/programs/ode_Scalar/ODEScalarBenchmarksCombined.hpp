/*
 *  Created on: 07 Mar 2023
 *      Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_ODE_SCALAR_COMBINED_HPP
#define SRC_INCLUDE_BENCHMARKS_ODE_SCALAR_COMBINED_HPP

#include <iostream>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "benchmarks/ShackODEScalarBenchmarks.hpp"
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include "benchmarks/ODEScalarBench_Default.hpp"

class ODEScalarBenchmarksCombined
{
public:
	sweet::ErrorBase error;
	ShackODEScalarBenchmarks *shackBenchmarks;

	ODEScalarBenchmarksCombined()	:
		shackBenchmarks(nullptr)
	{

	}


	/*
	 * Special function to register shacks for benchmarks.
	 *
	 * This is in particular important for the --help output function to include all benchmarks.
	 */
	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		shackBenchmarks = io_shackDict.getAutoRegistration<ShackODEScalarBenchmarks>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}

	void clear()
	{
		shackBenchmarks = nullptr;
	}


public:
	template <typename T>
	bool setupInitialConditions(
			T &o_u,
			sweet::ShackDictionary &io_shackDict
	)
	{
		if (!shackBenchmarks->validateNonStationaryODE())
			return error.forwardWithPositiveReturn(shackBenchmarks->error);

		if (shackBenchmarks->benchmark_name == "")
			return error.set("ODEScalarBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");

		if (shackBenchmarks->benchmark_name == "default")
		{
			ODEScalarBenchDefault default_bench;

			default_bench.shackRegistration(&io_shackDict);

			default_bench.setup();

			default_bench.setupBenchmark(
					o_u
			);

			return true;
		}

		shackBenchmarks->printProgramArguments();

		return error.set(std::string("Benchmark '")+shackBenchmarks->benchmark_name+ "' not found (or not available)");
	}
};



#endif
