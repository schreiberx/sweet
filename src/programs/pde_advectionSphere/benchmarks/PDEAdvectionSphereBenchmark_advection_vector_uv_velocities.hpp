/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_VELOCITIES_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_VELOCITIES_HPP_

#include "PDEAdvectionSphere_Benchmark_BaseInterface.hpp"
#include <ostream>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/VectorMath.hpp>
#include "PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump.hpp"



class PDEAdvectionSphereBenchmark_advection_vector_uv_velocities	:
	public PDEAdvectionSphere_Benchmark_BaseInterface
{
	PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	PDEAdvectionSphereBenchmark_advection_vector_uv_velocities()
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
			sweet::ShackDictionary *i_shackDict,
			sweet::SphereOperators *i_ops
	)
	{
		shackDict = i_shackDict;
		ops = i_ops;

		benchmark.setup(i_shackDict, i_ops);
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
	int getNumPrognosticFields()
	{
		return 2;
	}



	void getInitialState(
		std::vector<sweet::SphereData_Spectral> &o_prognostic_fields,
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 2, "Only a vectorial field (3 elements) supported for this benchmark!");

		const sweet::SphereDataConfig *sphereDataConfig = o_prognostic_fields[0].sphereDataConfig;

		sweet::SphereData_Spectral tmp(sphereDataConfig);
		sweet::SphereData_Spectral vrt(sphereDataConfig);
		sweet::SphereData_Spectral div(sphereDataConfig);
		benchmark.getInitialState_Spectral(tmp, vrt, div);

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
		o_prognostic_fields[0] = vrt;
		o_prognostic_fields[1] = div;
	}
};

#endif
