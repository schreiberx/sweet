/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_WILLIAMSON_1_COS_BELL_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_WILLIAMSON_1_COS_BELL_HPP_

#include "../swe_sphere_benchmarks/SWESphereBenchmark_williamson_1_advection_cos_bell.hpp"

#include "BenchmarksSphereAdvection_interface.hpp"
#include <ostream>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>



class BenchmarksSphereAdvection_williamson_1_advection_cos_bell	: public BenchmarksSphereAdvection_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


	SWESphereBenchmark_williamson_1_advection_cos_bell benchmark;

public:
	BenchmarksSphereAdvection_williamson_1_advection_cos_bell()
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
		simVars = i_simVars;
		ops = i_ops;

		benchmark.setup(i_simVars, i_ops);
	}


	std::string get_help()
	{
		return benchmark.get_help();
	}


	void get_initial_state(
		std::vector<SphereData_Spectral*> &o_prognostic_fields,
		SphereData_Physical &o_u,
		SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 1, "Only scalar field supported for this benchmark!");

		get_initial_state(*o_prognostic_fields[0], o_u, o_v);
	}



	virtual void get_initial_state(
			SphereData_Spectral &o_phi_pert,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v
	)
	{
		SphereData_Spectral vrt(o_phi_pert.sphereDataConfig);
		SphereData_Spectral div(o_phi_pert.sphereDataConfig);

		benchmark.get_initial_state(o_phi_pert, vrt, div);

		ops->vrtdiv_to_uv(vrt, div, o_u, o_v);
	}
};

#endif
