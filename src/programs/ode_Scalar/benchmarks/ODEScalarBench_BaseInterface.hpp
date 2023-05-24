/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_BASEINTERFACE_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_BENCHMARKS_ODESCALARBENCH_BASEINTERFACE_HPP_


#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "ShackODEScalarBenchmarks.hpp"
#include "../ShackODEScalar.hpp"

class ODEScalarBench_BaseInterface
{
public:
	sweet::ErrorBase error;

	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackODEScalar *shackODEScalar;
	ShackODEScalarBenchmarks *shackBenchmarks;

	ODEScalarBench_BaseInterface() :
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackODEScalar(nullptr),
		shackBenchmarks(nullptr)
	{
	}

	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackODEScalar = io_shackDict->getAutoRegistration<ShackODEScalar>();
		shackBenchmarks = io_shackDict->getAutoRegistration<ShackODEScalarBenchmarks>();

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*io_shackDict);

		return true;
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

