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
#include "benchmarks/PDEAdvectionSphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include <sweet/core/sphere/Sphere.hpp>



class PDEAdvectionSphereBenchmarksCombined
{
	std::vector<PDEAdvectionSphereBenchmarks_BaseInterface*> _registered_benchmarks;

	sweet::ShackDictionary *shackDict;
	ShackPDEAdvectionSphereBenchmarks *shackBenchmarks;

public:
	sweet::ErrorBase error;

	// Sphere operators
	sweet::SphereOperators *ops;

	PDEAdvectionSphereBenchmarks_BaseInterface *benchmark = nullptr;

	PDEAdvectionSphereBenchmarksCombined();

	bool setup_1_registerAllBenchmark();

	bool setup_2_shackRegistration(
		sweet::ShackDictionary *io_shackDict
	);

	bool setup_3_benchmarkDetection();

	bool setup_4_benchmarkSetup_1_withoutOps();

	bool setup_5_benchmarkSetup_2_withOps(
			sweet::SphereOperators *io_ops
	);

	void clear();

private:
	void _benchmarksFreeAll(
			PDEAdvectionSphereBenchmarks_BaseInterface *skip_this = nullptr
	);

public:
	void printAvailableBenchmarks(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);

public:
	~PDEAdvectionSphereBenchmarksCombined();
};


#endif
