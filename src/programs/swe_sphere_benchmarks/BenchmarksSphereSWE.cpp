/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "BenchmarksSphereSWE.hpp"
#include "BenchmarkSphereSWE_zero.hpp"
#include "SWESphereBenchmark_galewsky.hpp"
#include "SWESphereBenchmark_gaussian_bump.hpp"
#include "SWESphereBenchmark_gaussian_bumps_pvd.hpp"
#include "SWESphereBenchmark_gaussian_bumps_test_cases.hpp"
#include "SWESphereBenchmark_three_gaussian_bumps.hpp"
#include "SWESphereBenchmark_williamson_1_advection_cos_bell.hpp"
#include "SWESphereBenchmark_williamson_1_advection_gauss_bump.hpp"
#include "SWESphereBenchmark_williamson_2_geostrophic_balance.hpp"
#include "SWESphereBenchmark_williamson_2_geostrophic_balance_linear.hpp"
#include "SWESphereBenchmark_williamson_5_flow_over_mountain.hpp"
#include "SWESphereBenchmark_williamson_6_rossby_haurwitz_wave.hpp"
#include "SWESphereBenchmark_barotropic_vort_modes.hpp"




void BenchmarksSphereSWE::setup(
		SimulationVariables& io_simVars,
		SphereOperators_SphereData& io_op,
		const std::string& i_benchmark_name
)
{
	simVars = &io_simVars;
	op = &io_op;

	reset();

	std::string benchmark_name;

	if (i_benchmark_name == "")
		benchmark_name = simVars->benchmark.benchmark_name;
	else
		benchmark_name = i_benchmark_name;

	if (benchmark_name == "")
	{
		printAvailableBenchmarks();
		std::cout << std::endl;
		SWEETError("Please choose benchmark");
	}

	benchmarks_register_all();

	/*
	 * Find right one
	 */
	master = nullptr;


	for (std::size_t i = 0; i < registered_benchmarks.size(); i++)
	{
		SWESphereBenchmarks_interface *ts = registered_benchmarks[i];

		if (ts->implements_benchmark(benchmark_name))
		{
			if (master != nullptr)
			{
				std::cout << "Processing " << i+1 << "th element" << std::endl;
				SWEETError(std::string("Duplicate implementation for benchmark ") + benchmark_name);
			}

			std::cout << "Benchmark detection: found match with benchmark id " << i+1 << std::endl;
			ts->setup(&io_simVars, &io_op);
			master = ts;
		}
	}

	// Found integrator, freeing others
	benchmarks_free_all(master);

	if (master == nullptr)
	{
		printAvailableBenchmarks();
		std::cout << std::endl;
		SWEETError(std::string("No valid --benchmark-name '") + benchmark_name + std::string("' provided"));
	}

}


BenchmarksSphereSWE::~BenchmarksSphereSWE()
{
	reset();
}


void BenchmarksSphereSWE::reset()
{
	delete master;
	master = nullptr;

	benchmarks_free_all();
}


void BenchmarksSphereSWE::benchmarks_register_all()
{
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_gaussian_bumps_test_cases));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_gaussian_bump));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_gaussian_bumps_pvd));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_galewsky));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_three_gaussian_bumps));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_williamson_1_advection_cos_bell));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_williamson_1_advection_gauss_bump));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_williamson_2_geostrophic_balance));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_williamson_2_geostrophic_balance_linear));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_williamson_5_flow_over_mountain));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_williamson_6_rossby_haurwitz_wave));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new SWESphereBenchmark_barotropic_vort_modes));
	registered_benchmarks.push_back(static_cast<SWESphereBenchmarks_interface*>(new BenchmarkSphereSWE_zero));
}



void BenchmarksSphereSWE::benchmarks_free_all(SWESphereBenchmarks_interface *skip_this)
{
	for (std::size_t i = 0; i < registered_benchmarks.size(); i++)
	{
		SWESphereBenchmarks_interface *ts = registered_benchmarks[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	registered_benchmarks.clear();
}



BenchmarksSphereSWE::BenchmarksSphereSWE()	:
	simVars(nullptr),
	op(nullptr)
{
}



void BenchmarksSphereSWE::printAvailableBenchmarks()
{
	benchmarks_register_all();
	for (std::size_t i = 0; i < registered_benchmarks.size(); i++)
	{
		std::cout << registered_benchmarks[i]->get_help();
		std::cout << std::endl;
	}
	benchmarks_free_all();
}

