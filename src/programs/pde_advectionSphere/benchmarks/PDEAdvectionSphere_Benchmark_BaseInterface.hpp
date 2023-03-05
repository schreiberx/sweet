/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_INTERFACE_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_INTERFACE_HPP_

#include <sweet/core/ScalarDataArray.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "../time/ShackPDEAdvectionSphereTimeDiscretization.hpp"
#include "../benchmarks/ShackPDEAdvectionSphereBenchmarks.hpp"


class PDEAdvectionSphere_Benchmark_BaseInterface
{
public:
	sweet::ErrorBase error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	ShackPDEAdvectionSphereTimeDiscretization *shackPDEAdvTimeDisc;
	ShackPDEAdvectionSphereBenchmarks *shackPDEAdvBenchmark;

	sweet::SphereOperators *ops;

	PDEAdvectionSphere_Benchmark_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphereDataOps(nullptr),
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
		shackSphereDataOps = io_shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
		shackPDEAdvTimeDisc = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphereTimeDiscretization>();
		shackPDEAdvBenchmark = io_shackDict->getAutoRegistration<ShackPDEAdvectionSphereBenchmarks>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}


public:
	virtual void setup(
			sweet::ShackDictionary *io_shackDict,
			sweet::SphereOperators *io_ops
		)
	{
		shackDict = io_shackDict;
		ops = io_ops;
	}


	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;


	virtual std::string get_help() = 0;

	virtual void getInitialState(
		std::vector<sweet::SphereData_Spectral> &o_phi,
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	) = 0;


#if 0
	virtual void get_reference_state(
		std::vector<sweet::SphereData_Spectral*> &o_phi_pert,
		double i_timestamp
	)
	{
		SWEETError("'get_reference_state' for multiple prognostic variables not implemented for this benchmark");
	}


	virtual void get_varying_velocities(
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v,
		double i_timestamp
	)
	{
		SWEETError("'get_varying_velocities' not implemented for this benchmark");
	}


	virtual bool has_time_varying_state()
	{
		return false;
	}
#endif


	/*
	 * Return number of prognostic fields to be used
	 */
	virtual int getNumPrognosticFields()
	{
		return 1;
	}


	virtual ~PDEAdvectionSphere_Benchmark_BaseInterface()
	{
	}
};

#endif
