/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_ADV_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include "ShackPDEAdvectionPlaneTimeDiscretization.hpp"
#include "../benchmarks/ShackPDEAdvectionPlaneBenchmarks.hpp"


class PDEAdvectionPlaneTS_BaseInterface
{
public:
	sweet::ErrorBase error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	ShackPDEAdvectionPlaneTimeDiscretization *shackPDEAdvTimeDisc;
	ShackPDEAdvectionPlaneBenchmarks *shackPDEAdvBenchmark;

	sweet::PlaneOperators *ops;

	PDEAdvectionPlaneTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackPlaneDataOps(nullptr),
		shackPDEAdvTimeDisc(nullptr),
		shackPDEAdvBenchmark(nullptr),
		ops(nullptr)
	{

	}

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackPlaneDataOps = io_shackDict->getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackPDEAdvTimeDisc = io_shackDict->getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();
		shackPDEAdvBenchmark = io_shackDict->getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool setup(
			sweet::PlaneOperators *io_ops
	)
	{
		ops = io_ops;
		return true;
	}


public:
	virtual void runTimestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	) = 0;

	~PDEAdvectionPlaneTS_BaseInterface() {}
};



#endif
