/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_HPP_

#include "BenchmarksSphereAdvection_interface.hpp"
#include <sweet/core/SimulationVariables.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>


#if !SWEET_USE_SPHERE_SPECTRAL_SPACE
	#error "SWEET_USE_SPHERE_SPECTRAL_SPACE not enabled"
#endif



class BenchmarksSphereAdvection
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

	BenchmarksSphereAdvection_interface *master = nullptr;

	std::vector<BenchmarksSphereAdvection_interface*> registered_benchmarks;


	void reset();

	void benchmarks_register_all();
	void benchmarks_free_all(BenchmarksSphereAdvection_interface *skip_this = nullptr);
	BenchmarksSphereAdvection();

public:
	void printAvailableBenchmarks();
};



#endif /* SRC_INCLUDE_benchmarks_swe_sphere_SWESPHEREBENCHMARKS_HPP_ */
