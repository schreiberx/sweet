/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDEAdvectionSphereBenchmarksCombined.hpp"

#include "benchmarks/PDEAdvectionSphereBenchmark_advection_vector_uv_velocities.hpp"
#include "benchmarks/PDEAdvectionSphereBenchmark_advection_vector_uv_gauss_bumps.hpp"
#include "benchmarks/PDEAdvectionSphereBenchmark_advection_vector_3d_normal_vectors.hpp"
#include "benchmarks/PDEAdvectionSphereBenchmark_nair_lauritzen_sl.hpp"
#include "benchmarks/PDEAdvectionSphereBenchmark_williamson_1_advection_cos_bell.hpp"
#include "benchmarks/PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump.hpp"
#include "benchmarks/PDEAdvectionSphereBenchmark_zero.hpp"

#include "time/PDEAdvectionSphereTS_BaseInterface.hpp"



PDEAdvectionSphereBenchmarksCombined::PDEAdvectionSphereBenchmarksCombined()	:
	shackDict(nullptr),
	ops(nullptr),
	benchmark(nullptr)
{
}


bool PDEAdvectionSphereBenchmarksCombined::setup_1_registerAllBenchmark()
{
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere_Benchmark_BaseInterface*>(new PDEAdvectionSphereBenchmark_zero));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere_Benchmark_BaseInterface*>(new PDEAdvectionSphereBenchmark_nair_lauritzen_sl));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere_Benchmark_BaseInterface*>(new PDEAdvectionSphereBenchmark_advection_vector_uv_velocities));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere_Benchmark_BaseInterface*>(new PDEAdvectionSphereBenchmark_advection_vector_uv_gauss_bumps));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere_Benchmark_BaseInterface*>(new PDEAdvectionSphereBenchmark_advection_vector_3d_normal_vectors));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere_Benchmark_BaseInterface*>(new PDEAdvectionSphereBenchmark_williamson_1_advection_cos_bell));
	_registered_benchmarks.push_back(static_cast<PDEAdvectionSphere_Benchmark_BaseInterface*>(new PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump));

	return true;
}



bool PDEAdvectionSphereBenchmarksCombined::setup_2_shackRegistration(
	sweet::ShackDictionary *io_shackDict
)
{
	shackDict = io_shackDict;

	shackBenchmarks = shackDict->getAutoRegistration<ShackPDEAdvectionSphereBenchmarks>();

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		_registered_benchmarks[i]->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*_registered_benchmarks[i]);
	}

	return true;
}



bool PDEAdvectionSphereBenchmarksCombined::setup_3_benchmarkDetection()
{
	assert(benchmark == nullptr);

	if (shackBenchmarks->benchmark_name == "")
	{
		printAvailableBenchmarks();
		return error.set("Please choose benchmark with --benchmark-name=...");
	}

	/*
	 * Find right one
	 */
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		PDEAdvectionSphere_Benchmark_BaseInterface *ts = _registered_benchmarks[i];

		if (ts->implements_benchmark(shackBenchmarks->benchmark_name))
		{
			if (benchmark != nullptr)
			{
				//std::cout << "Processing " << i+1 << "th element" << std::endl;
				error.set("Duplicate implementation for benchmark "+shackBenchmarks->benchmark_name);
			}

			//std::cout << "Benchmark detection: found match with benchmark id " << i+1 << std::endl;
			benchmark = ts;
		}
	}

	if (benchmark == nullptr)
	{
		printAvailableBenchmarks();
		return error.set("No valid --benchmark-name '"+shackBenchmarks->benchmark_name+"' provided");
	}

	return true;
}


bool PDEAdvectionSphereBenchmarksCombined::setup_4_benchmarkSetup_1_withoutOps()
{
	assert(benchmark != nullptr);

	benchmark->setup_1_shackBenchmarkData(shackDict);
	ERROR_FORWARD_WITH_RETURN_BOOLEAN(*benchmark);

	return true;
}


bool PDEAdvectionSphereBenchmarksCombined::setup_5_benchmarkSetup_2_withOps(
		sweet::SphereOperators* io_ops
)
{
	ops = io_ops;

	benchmark->setup_2_withOps(shackDict, ops);
	ERROR_FORWARD_WITH_RETURN_BOOLEAN(*benchmark);
}

void PDEAdvectionSphereBenchmarksCombined::clear()
{
	benchmark = nullptr;
	shackDict = nullptr;
	ops = nullptr;

	_benchmarksFreeAll();
}



void PDEAdvectionSphereBenchmarksCombined::_benchmarksFreeAll(
		PDEAdvectionSphere_Benchmark_BaseInterface *skip_this
)
{
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		PDEAdvectionSphere_Benchmark_BaseInterface *ts = _registered_benchmarks[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	_registered_benchmarks.clear();
}


void PDEAdvectionSphereBenchmarksCombined::printAvailableBenchmarks(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		std::cout << _registered_benchmarks[i]->get_help();
		std::cout << std::endl;
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

}


PDEAdvectionSphereBenchmarksCombined::~PDEAdvectionSphereBenchmarksCombined()
{
	_benchmarksFreeAll();
}
