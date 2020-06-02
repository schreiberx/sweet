/*
 * Author: Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_VECTOR_UV_ADVECTION_GAUSS_BUMP_HPP_
#define SRC_BENCHMARKS_SPHERE_VECTOR_UV_ADVECTION_GAUSS_BUMP_HPP_

#include "BenchmarksSphereAdvection_interface.hpp"
#include <ostream>
#include <sweet/SWEETMath.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>

#include "../swe_sphere_benchmarks/SWESphereBenchmark_williamson_1_advection_gauss_bump.hpp"



class BenchmarksSphereAdvection_vector_uv_advection_gauss_bump	: public BenchmarksSphereAdvection_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;

	SWESphereBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	BenchmarksSphereAdvection_vector_uv_advection_gauss_bump()
	{
	}

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		return (
				i_benchmark_name == "vector_uv_advection_gauss_bump"	||
				false
			);
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
		std::ostringstream stream;

		stream << " * Advection test case with 2d vector in lat-lon space:" << std::endl;
		stream << "    + 'vector_uv_advection_gauss_bump'" << std::endl;

		return stream.str();
	}


	/*
	 * Return number of prognostic fields to be used
	 */
	int get_num_prognostic_fields()
	{
		return 2;
	}



	void get_initial_state(
		std::vector<SphereData_Spectral*> &o_prognostic_fields,
		SphereData_Physical &o_u,
		SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 2, "Only a vectorial field (3 elements) supported for this benchmark!");

		const SphereData_Config *sphereDataConfig = o_prognostic_fields[0]->sphereDataConfig;

		SphereData_Spectral tmp(sphereDataConfig);
		SphereData_Spectral vrt(sphereDataConfig);
		SphereData_Spectral div(sphereDataConfig);
		benchmark.get_initial_state(tmp, vrt, div);

		*o_prognostic_fields[0] = vrt;
		*o_prognostic_fields[1] = div;

		ops->vrtdiv_to_uv(vrt, div, o_u, o_v);
	}
};

#endif
