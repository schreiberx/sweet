/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_VELOCITIES_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_VELOCITIES_HPP_

#include "BenchmarksSphereAdvection_interface.hpp"
#include <ostream>
#include <sweet/core/SimulationVariables.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>
#include <sweet/core/SWEETVectorMath.hpp>
#include "../swe_sphere_benchmarks/SWESphereBenchmark_williamson_1_advection_gauss_bump.hpp"



class BenchmarksSphereAdvection_advection_vector_uv_velocities	: public BenchmarksSphereAdvection_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;

	SWESphereBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	BenchmarksSphereAdvection_advection_vector_uv_velocities()
	{
	}

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		return (
				i_benchmark_name == "advection_vector_uv_velocities"	||
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
		stream << "    + 'advection_vector_uv_velocities'" << std::endl;

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
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 2, "Only a vectorial field (3 elements) supported for this benchmark!");

		const sweet::SphereData_Config *sphereDataConfig = o_prognostic_fields[0]->sphereDataConfig;

		SphereData_Spectral tmp(sphereDataConfig);
		SphereData_Spectral vrt(sphereDataConfig);
		SphereData_Spectral div(sphereDataConfig);
		benchmark.get_initial_state(tmp, vrt, div);

		/*
		 * Setup velocity field
		 */
		ops->vrtdiv_to_uv(vrt, div, o_u, o_v);

		/*
		 * IMPORTANT INFORMATION:
		 * Here, we like to use the velocities as prognostic fields.
		 *
		 * The prognostic fields here are the vrt/div of the velocities!!!
		 */
		*o_prognostic_fields[0] = vrt;
		*o_prognostic_fields[1] = div;
	}
};

#endif
