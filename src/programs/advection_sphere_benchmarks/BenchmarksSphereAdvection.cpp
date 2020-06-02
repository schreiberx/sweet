/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "BenchmarksSphereAdvection.hpp"
#include "BenchmarksSphereAdvection_nair_lauritzen_sl.hpp"
#include "BenchmarksSphereAdvection_vector_3d_advection_gauss_bump.hpp"
#include "BenchmarksSphereAdvection_vector_uv_advection_gauss_bump.hpp"
#include "BenchmarksSphereAdvection_williamson_1_advection_cos_bell.hpp"
#include "BenchmarksSphereAdvection_williamson_1_advection_gauss_bump.hpp"
#include "BenchmarksSphereAdvection_zero.hpp"




void BenchmarksSphereAdvection::setup(
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
		BenchmarksSphereAdvection_interface *ts = registered_benchmarks[i];

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


void BenchmarksSphereAdvection::reset()
{
	delete master;
	master = nullptr;

	benchmarks_free_all();
}


void BenchmarksSphereAdvection::benchmarks_register_all()
{
	registered_benchmarks.push_back(static_cast<BenchmarksSphereAdvection_interface*>(new BenchmarksSphereAdvection_zero));
	registered_benchmarks.push_back(static_cast<BenchmarksSphereAdvection_interface*>(new BenchmarksSphereAdvection_nair_lauritzen_sl));
	registered_benchmarks.push_back(static_cast<BenchmarksSphereAdvection_interface*>(new BenchmarksSphereAdvection_vector_3d_advection_gauss_bump));
	registered_benchmarks.push_back(static_cast<BenchmarksSphereAdvection_interface*>(new BenchmarksSphereAdvection_vector_uv_advection_gauss_bump));
	registered_benchmarks.push_back(static_cast<BenchmarksSphereAdvection_interface*>(new BenchmarksSphereAdvection_williamson_1_advection_cos_bell));
	registered_benchmarks.push_back(static_cast<BenchmarksSphereAdvection_interface*>(new BenchmarksSphereAdvection_williamson_1_advection_gauss_bump));
}



void BenchmarksSphereAdvection::benchmarks_free_all(BenchmarksSphereAdvection_interface *skip_this)
{
	for (std::size_t i = 0; i < registered_benchmarks.size(); i++)
	{
		BenchmarksSphereAdvection_interface *ts = registered_benchmarks[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	registered_benchmarks.clear();
}



BenchmarksSphereAdvection::BenchmarksSphereAdvection()	:
	simVars(nullptr),
	op(nullptr)
{
}



void BenchmarksSphereAdvection::printAvailableBenchmarks()
{
	benchmarks_register_all();
	for (std::size_t i = 0; i < registered_benchmarks.size(); i++)
	{
		std::cout << registered_benchmarks[i]->get_help();
		std::cout << std::endl;
	}
	benchmarks_free_all();
}

