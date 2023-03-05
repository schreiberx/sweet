/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_HPP_


#include <iostream>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "benchmarks/ShackPDEAdvectionSphereBenchmarks.hpp"
#include "benchmarks/PDEAdvectionSphere_Benchmark_BaseInterface.hpp"
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include <sweet/core/sphere/Sphere.hpp>



class PDEAdvectionSphereBenchmarksCombined
{
	std::vector<PDEAdvectionSphere_Benchmark_BaseInterface*> _registered_benchmarks;

	sweet::ShackDictionary *shackDict;
	ShackPDEAdvectionSphereBenchmarks *shackBenchmarks;

public:
	sweet::ErrorBase error;

	// Sphere operators
	sweet::SphereOperators *ops;

	PDEAdvectionSphere_Benchmark_BaseInterface *master = nullptr;

	PDEAdvectionSphereBenchmarksCombined();

	bool setup_registerBenchmark(
			std::ostream &o_ostream = std::cout,
			const std::string &i_prefix = ""
	);

	bool setup_benchmark(
			sweet::SphereOperators *io_ops
	);

	void clear();


	/*
	 * Special function to register shacks for benchmarks.
	 *
	 * This is in particular important for the --help output function to include all benchmarks.
	 */
	bool shackRegistration(
		sweet::ShackDictionary &io_shackDict
	);

	void benchmarks_register_all();

	void benchmarks_free_all(
			PDEAdvectionSphere_Benchmark_BaseInterface *skip_this = nullptr
	);

public:
	void printAvailableBenchmarks(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);
};



#endif
