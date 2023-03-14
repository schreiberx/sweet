/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_VELOCITIES_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_VELOCITIES_HPP_

#include "PDEAdvectionSphereBenchmarks_BaseInterface.hpp"
#include <ostream>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/VectorMath.hpp>
#include "PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump.hpp"



class PDEAdvectionSphereBenchmark_advection_vector_uv_velocities	:
	public PDEAdvectionSphereBenchmarks_BaseInterface
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



	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		PDEAdvectionSphereBenchmarks_BaseInterface::shackRegistration(io_shackDict);
		benchmark.shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmark);
		return true;
	}

	void setup_1_shackData()
	{
		benchmark.setup_1_shackData();
	}

	void setup_2_withOps(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;

		benchmark.setup_2_withOps(ops);
	}



	std::string printHelp()
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

		const sweet::SphereData_Config *sphereDataConfig = o_prognostic_fields[0].sphereDataConfig;

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
