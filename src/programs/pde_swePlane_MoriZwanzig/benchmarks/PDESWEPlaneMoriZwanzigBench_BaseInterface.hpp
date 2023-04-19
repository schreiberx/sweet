/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_SWE_PLANE_BENCHMARKS_BASEINTERFACE_HPP_
#define SRC_PROGRAMS_PDE_SWE_PLANE_BENCHMARKS_BASEINTERFACE_HPP_


#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "ShackPDESWEPlaneMoriZwanzigBenchmarks.hpp"
///#include "ShackPDESWEPlaneBench_PolvaniBench.hpp"
#include "../ShackPDESWEPlaneMoriZwanzig.hpp"

class PDESWEPlaneMoriZwanzigBench_BaseInterface
{
public:
	sweet::ErrorBase error;

	sweet::ShackDictionary *shackDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWEPlaneMoriZwanzigBenchmarks *shackBenchmarks;
	ShackPDESWEPlaneMoriZwanzig *shackPDESWEPlane;
	//////ShackPDESWEPlaneBench_PolvaniBench *shackPDESWEPlaneBench_PolvaniBench;

	sweet::PlaneOperators *ops;
	sweet::PlaneData_Config *planeDataConfig;


	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackPlaneDataOps = io_shackDict->getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackBenchmarks = io_shackDict->getAutoRegistration<ShackPDESWEPlaneMoriZwanzigBenchmarks>();
		shackPDESWEPlane = io_shackDict->getAutoRegistration<ShackPDESWEPlaneMoriZwanzig>();
		////shackPDESWEPlaneBench_PolvaniBench = io_shackDict->getAutoRegistration<ShackPDESWEPlaneBench_PolvaniBench>();

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(*io_shackDict);
	}

	virtual bool setup(
		sweet::PlaneOperators *io_ops,
		sweet::PlaneData_Config *i_planeDataConfig
	)
	{
		ops = io_ops;
		planeDataConfig = i_planeDataConfig;

		return true;
	}

	virtual bool setupBenchmark(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v
	) = 0;
};


#endif
