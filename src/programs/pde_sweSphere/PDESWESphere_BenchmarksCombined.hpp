/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARK_SPHERE_SWE_HPP_
#define SRC_BENCHMARK_SPHERE_SWE_HPP_

#include "benchmarks/PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include "benchmarks/ShackPDESWESphereBenchmarks.hpp"

#if !SWEET_USE_SPHERE_SPECTRAL_SPACE && SWEET_PARAREAL
	#error "SWEET_USE_SPHERE_SPECTRAL_SPACE not enabled"
#endif



#include <iostream>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "benchmarks/ShackPDESWESphereBenchmarks.hpp"
#include "benchmarks/PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include <sweet/core/sphere/Sphere.hpp>



class PDESWESphere_BenchmarksCombined
{
	std::vector<PDESWESphereBenchmarks_BaseInterface*> _registered_benchmarks;

	sweet::ShackDictionary *shackDict;
	ShackPDESWESphereBenchmarks *shackBenchmarks;

public:
	sweet::ErrorBase error;

	// Sphere operators
	sweet::SphereOperators *ops;

	PDESWESphereBenchmarks_BaseInterface *benchmark = nullptr;

	PDESWESphere_BenchmarksCombined();

	bool setup_1_registerAllBenchmark();

	bool setup_2_shackRegistration(
		sweet::ShackDictionary *io_shackDict
	);

	bool setup_3_benchmarkDetection(
			const std::string &i_benchmark_name = ""
	);

	bool setup_4_benchmarkSetup_1_withoutOps();

	bool setup_5_benchmarkSetup_2_withOps(
			sweet::SphereOperators *io_ops
	);

	void clear();

private:
	void _benchmarksFreeAll(
			PDESWESphereBenchmarks_BaseInterface *skip_this = nullptr
	);

public:
	void printAvailableBenchmarks(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);

public:
	~PDESWESphere_BenchmarksCombined();
};



#endif
