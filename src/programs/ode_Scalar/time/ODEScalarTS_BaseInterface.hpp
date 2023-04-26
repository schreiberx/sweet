/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "ShackODEScalarTimeDiscretization.hpp"
#include "../benchmarks/ShackODEScalarBenchmarks.hpp"
#include "../ShackODEScalar.hpp"

class ODEScalarTS_BaseInterface
{
public:
	sweet::ErrorBase error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackODEScalar *shackODEScalar;
	ShackODEScalarTimeDiscretization *shackODEScalarTimeDisc;
	ShackODEScalarBenchmarks *shackODEScalarBenchmark;

	ODEScalarTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackODEScalar(nullptr),
		shackODEScalarTimeDisc(nullptr),
		shackODEScalarBenchmark(nullptr)
	{

	}

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackODEScalar = io_shackDict->getAutoRegistration<ShackODEScalar>();
		shackODEScalarTimeDisc = io_shackDict->getAutoRegistration<ShackODEScalarTimeDiscretization>();
		shackODEScalarBenchmark = io_shackDict->getAutoRegistration<ShackODEScalarBenchmarks>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool setup(
	)
	{
		return true;
	}


public:
	virtual void runTimestep(
			double &io_u,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	) = 0;

	~ODEScalarTS_BaseInterface() {}
};



#endif
