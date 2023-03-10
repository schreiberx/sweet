/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere_BenchmarksCombined.hpp"
#include "benchmarks/PDESWESphereBenchmark_zero.hpp"
#include "benchmarks/PDESWESphereBenchmark_galewsky.hpp"
#include "benchmarks/PDESWESphereBenchmark_gaussian_bump.hpp"
#include "benchmarks/PDESWESphereBenchmark_gaussian_bumps_pvd.hpp"
#include "benchmarks/PDESWESphereBenchmark_gaussian_bumps_test_cases.hpp"
#include "benchmarks/PDESWESphereBenchmark_three_gaussian_bumps.hpp"
#include "benchmarks/PDESWESphereBenchmark_williamson_1_advection_cos_bell.hpp"
#include "benchmarks/PDESWESphereBenchmark_williamson_1_advection_gauss_bump.hpp"
#include "benchmarks/PDESWESphereBenchmark_williamson_2_geostrophic_balance.hpp"
#include "benchmarks/PDESWESphereBenchmark_williamson_2_geostrophic_balance_linear.hpp"
#include "benchmarks/PDESWESphereBenchmark_williamson_5_flow_over_mountain.hpp"
#include "benchmarks/PDESWESphereBenchmark_williamson_6_rossby_haurwitz_wave.hpp"
#include "benchmarks/PDESWESphereBenchmark_barotropic_vort_modes.hpp"



PDESWESphere_BenchmarksCombined::PDESWESphere_BenchmarksCombined()	:
	shackDict(nullptr),
	ops(nullptr),
	benchmark(nullptr)
{
}


bool PDESWESphere_BenchmarksCombined::setup_1_registerAllBenchmark()
{
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_gaussian_bumps_test_cases));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_gaussian_bump));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_gaussian_bumps_pvd));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_galewsky));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_three_gaussian_bumps));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_williamson_1_advection_cos_bell));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_williamson_1_advection_gauss_bump));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_williamson_2_geostrophic_balance));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_williamson_2_geostrophic_balance_linear));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_williamson_5_flow_over_mountain));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_williamson_6_rossby_haurwitz_wave));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_barotropic_vort_modes));
	_registered_benchmarks.push_back(static_cast<PDESWESphereBenchmarks_BaseInterface*>(new PDESWESphereBenchmark_zero));

	return true;
}



bool PDESWESphere_BenchmarksCombined::setup_2_shackRegistration(
	sweet::ShackDictionary *io_shackDict
)
{
	shackDict = io_shackDict;

	shackBenchmarks = shackDict->getAutoRegistration<ShackPDESWESphereBenchmarks>();

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		_registered_benchmarks[i]->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*_registered_benchmarks[i]);
	}

	return true;
}



bool PDESWESphere_BenchmarksCombined::setup_3_benchmarkDetection()
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
		PDESWESphereBenchmarks_BaseInterface *ts = _registered_benchmarks[i];

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


bool PDESWESphere_BenchmarksCombined::setup_4_benchmarkSetup_1_withoutOps()
{
	assert(benchmark != nullptr);

	benchmark->setup_1_shackData();
	ERROR_FORWARD_WITH_RETURN_BOOLEAN(*benchmark);

	return true;
}


bool PDESWESphere_BenchmarksCombined::setup_5_benchmarkSetup_2_withOps(
		sweet::SphereOperators* io_ops
)
{
	ops = io_ops;

	benchmark->setup_2_withOps(ops);
	ERROR_FORWARD_WITH_RETURN_BOOLEAN(*benchmark);
}

void PDESWESphere_BenchmarksCombined::clear()
{
	benchmark = nullptr;
	shackDict = nullptr;
	ops = nullptr;

	_benchmarksFreeAll();
}



void PDESWESphere_BenchmarksCombined::_benchmarksFreeAll(
		PDESWESphereBenchmarks_BaseInterface *skip_this
)
{
	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		PDESWESphereBenchmarks_BaseInterface *ts = _registered_benchmarks[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	_registered_benchmarks.clear();
}


void PDESWESphere_BenchmarksCombined::printAvailableBenchmarks(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	for (std::size_t i = 0; i < _registered_benchmarks.size(); i++)
	{
		std::cout << _registered_benchmarks[i]->printHelp();
		std::cout << std::endl;
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Benchmarks (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

}


PDESWESphere_BenchmarksCombined::~PDESWESphere_BenchmarksCombined()
{
	clear();
}
