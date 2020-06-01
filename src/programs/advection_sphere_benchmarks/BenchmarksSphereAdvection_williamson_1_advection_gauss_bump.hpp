/*
 * Author: Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_WILLIAMSON_1_GAUSS_BUMP_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_WILLIAMSON_1_GAUSS_BUMP_HPP_

#include "../swe_sphere_benchmarks/SWESphereBenchmark_williamson_1_advection_gauss_bump.hpp"

#include "BenchmarksSphereAdvection_interface.hpp"
#include <ostream>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>



class BenchmarksSphereAdvection_williamson_1_advection_gauss_bump	: public BenchmarksSphereAdvection_interface
{
	SWESphereBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	BenchmarksSphereAdvection_williamson_1_advection_gauss_bump()
	{
	}

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		return benchmark.implements_benchmark(i_benchmark_name);
	}


	void setup(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
	)
	{
		benchmark.setup(i_simVars, i_ops);
	}


	std::string get_help()
	{
		return benchmark.get_help();
	}


	void get_initial_state(
		std::vector<SphereData_Spectral*> &o_prognostic_fields,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{
		SWEETDebugAssert(o_prognostic_fields.size() == 1, "Only scalar field supported for this benchmark!");

		get_initial_state(*o_prognostic_fields[0], o_vrt, o_div);
	}



	virtual void get_initial_state(
			SphereData_Spectral &o_phi_pert,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div
	)
	{
		benchmark.get_initial_state(o_phi_pert, o_vrt, o_div);
	}
};

#endif
