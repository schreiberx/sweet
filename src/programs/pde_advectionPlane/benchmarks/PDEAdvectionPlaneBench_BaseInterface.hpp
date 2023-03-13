/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_ADVECTIONPLANE_BENCHMARKS_PDEADVECTIONPLANEBENCH_BASEINTERFACE_HPP_
#define SRC_PROGRAMS_PDE_ADVECTIONPLANE_BENCHMARKS_PDEADVECTIONPLANEBENCH_BASEINTERFACE_HPP_


#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "ShackPDEAdvectionPlaneBenchmarks.hpp"


class PDEAdvectionPlaneBench_BaseInterface
{
public:
	sweet::ErrorBase error;

	sweet::ShackDictionary *shackDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDEAdvectionPlaneBenchmarks *shackBenchmarks;

	sweet::PlaneOperators *ops;


	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackPlaneDataOps = io_shackDict->getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackBenchmarks = io_shackDict->getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*io_shackDict);
	}

	virtual bool setup(
		sweet::PlaneOperators *io_ops
	)
	{
		ops = io_ops;

		return true;
	}

	virtual bool setupBenchmark(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v
	) = 0;
};


#endif
