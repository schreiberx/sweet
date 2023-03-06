/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_TS_INTERFACE_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_TS_INTERFACE_HPP_

#include <sweet/core/ScalarDataArray.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "../time/ShackPDESWESphereTimeDiscretization.hpp"
#include "../benchmarks/ShackPDESWESphereBenchmarks.hpp"
#include "../ShackPDESWESphere.hpp"


class PDESWESphereBenchmarks_BaseInterface
{
public:
	sweet::ErrorBase error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */
	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	ShackPDESWESphereTimeDiscretization *shackPDESWETimeDisc;
	ShackPDESWESphereBenchmarks *shackPDESWEBenchmark;
	ShackPDESWESphere *shackPDESWESphere;

	sweet::SphereOperators *ops;

	PDESWESphereBenchmarks_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphereDataOps(nullptr),
		shackPDESWETimeDisc(nullptr),
		shackPDESWEBenchmark(nullptr),
		shackPDESWESphere(nullptr),
		ops(nullptr)
	{
	}

	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackSphereDataOps = shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
		shackPDESWETimeDisc = shackDict->getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
		shackPDESWEBenchmark = shackDict->getAutoRegistration<ShackPDESWESphereBenchmarks>();
		shackPDESWESphere = shackDict->getAutoRegistration<ShackPDESWESphere>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*shackDict);

		return true;
	}

public:
	virtual void setup_1_shackData() = 0;

public:
	virtual void setup_2_withOps(
			sweet::SphereOperators *io_ops
	) = 0;

	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;


	virtual std::string printHelp() = 0;

	virtual void getInitialState(
		sweet::SphereData_Spectral &o_phi,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	) = 0;


	virtual void get_reference_state(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div,
		double i_timestamp
	)
	{
		SWEETError("Not implemented for this benchmark");
	}

	virtual bool has_time_varying_state()
	{
		return false;
	}

	virtual void clear() = 0;

	virtual ~PDESWESphereBenchmarks_BaseInterface()
	{
	}
};

#endif
