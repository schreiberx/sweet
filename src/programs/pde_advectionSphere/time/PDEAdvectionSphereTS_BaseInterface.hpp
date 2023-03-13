/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SPHERE_ADVECTION_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_SPHERE_ADVECTION_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/core/sphere/Sphere.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackSphereDataOps.hpp>
#include <sweet/core/time/ShackTimesteppingSemiLagrangianSphereData.hpp>
#include "ShackPDEAdvectionSphereTimeDiscretization.hpp"
#include "../benchmarks/ShackPDEAdvectionSphereBenchmarks.hpp"


class PDEAdvectionSphereTS_BaseInterface
{
public:
	sweet::ErrorBase error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	ShackPDEAdvectionSphereTimeDiscretization *shackPDEAdvectionTimeDisc;
	ShackPDEAdvectionSphereBenchmarks *shackPDEAdvBenchmark;
	sweet::ShackTimesteppingSemiLagrangianSphereData *shackSemiLagrangian;

	sweet::SphereOperators *ops;

	PDEAdvectionSphereTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphereDataOps(nullptr),
		shackPDEAdvectionTimeDisc(nullptr),
		shackPDEAdvBenchmark(nullptr),
		ops(nullptr)
	{

	}

	virtual bool shackRegistration(
		sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackSphereDataOps = io_shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
		shackPDEAdvectionTimeDisc = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphereTimeDiscretization>();
		shackPDEAdvBenchmark = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphereBenchmarks>();
		shackSemiLagrangian = io_shackDict->getAutoRegistration<sweet::ShackTimesteppingSemiLagrangianSphereData>();

		ERROR_CHECK_WITH_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

	virtual bool setup(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;
		return true;
	}

	/*
	 * Timestepping interface used by main timestepping loop
	 */
	virtual void runTimestep(
			std::vector<sweet::SphereData_Spectral> &io_prognostic_fields,	///< prognostic variables
			sweet::SphereData_Physical &io_u,
			sweet::SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	) = 0;

	virtual bool testImplementsTimesteppingMethod(
			const std::string &i_timestepping_method
		) = 0;

	virtual std::string getStringId() = 0;

	virtual void printImplementedTimesteppingMethods(
			std::ostream &o_ostream = std::cout,
			const std::string &i_prefix = ""
		) = 0;

	virtual ~PDEAdvectionSphereTS_BaseInterface()
	{
	}
};

#endif
