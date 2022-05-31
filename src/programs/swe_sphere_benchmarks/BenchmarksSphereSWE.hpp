/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_BENCHMARK_SPHERE_SWE_HPP_
#define SRC_BENCHMARK_SPHERE_SWE_HPP_

#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>

#if !SWEET_USE_SPHERE_SPECTRAL_SPACE && SWEET_PARAREAL
	#error "SWEET_USE_SPHERE_SPECTRAL_SPACE not enabled"
#endif



class BenchmarksSphereSWE
{
public:
	// Simulation variables
	SimulationVariables *simVars;

	// Sphere operators
	SphereOperators_SphereData *op;


	void setup(
			SimulationVariables& io_simVars,
			SphereOperators_SphereData& io_op,
			const std::string& i_benchmark_name = ""
	);

	SWESphereBenchmarks_interface *master = nullptr;

	std::vector<SWESphereBenchmarks_interface*> registered_benchmarks;


	void reset();

	void benchmarks_register_all();
	void benchmarks_free_all(SWESphereBenchmarks_interface *skip_this = nullptr);
	BenchmarksSphereSWE();
	~BenchmarksSphereSWE();

public:
	void printAvailableBenchmarks();
};



#endif /* SRC_INCLUDE_benchmarks_swe_sphere_SWESPHEREBENCHMARKS_HPP_ */
