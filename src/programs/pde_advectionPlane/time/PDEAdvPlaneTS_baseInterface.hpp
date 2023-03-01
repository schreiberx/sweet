/*
 * ADV_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
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


class PDEAdvPlaneTS_baseInterface
{
public:
	sweet::ErrorBase error;

	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	ShackPDEAdvectionPlaneTimeDiscretization *shackPDEAdvTimeDisc;
	ShackPDEAdvectionPlaneBenchmarks *shackPDEAdvBenchmark;

	bool registerShacks(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackPlaneDataOps = io_shackDict->getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackPDEAdvTimeDisc = io_shackDict->getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();
		shackPDEAdvBenchmark = io_shackDict->getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

public:
	virtual void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	) = 0;

#if 0
	virtual bool registerShacks(
			sweet::ShackDictionary *io_shackDict
	) = 0;
#endif

	~PDEAdvPlaneTS_baseInterface() {}
};



#endif
