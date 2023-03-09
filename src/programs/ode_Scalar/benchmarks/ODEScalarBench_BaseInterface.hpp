/*
 *  Created on: 07 Mar 2023
 *      Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_BASEINTERFACE_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_BASEINTERFACE_HPP_


#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "ShackODEScalarBenchmarks.hpp"

class ODEScalarBench_BaseInterface
{
public:
	sweet::ErrorBase error;

	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackODEScalarBenchmarks *shackBenchmarks;

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackBenchmarks = io_shackDict->getAutoRegistration<ShackODEScalarBenchmarks>();

		ERROR_FORWARD_WITH_RETURN_BOOLEAN(*io_shackDict);
	}

	virtual bool setup(
	)
	{
		return true;
	}

	virtual bool setupBenchmark(
			double &o_u
	) = 0;
};


#endif
