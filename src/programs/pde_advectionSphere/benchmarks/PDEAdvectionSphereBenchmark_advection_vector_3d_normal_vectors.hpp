/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_VECTOR_3D_ADVECTION_VECTOR_3D_NORMAL_VECTORS_
#define SRC_BENCHMARKS_SPHERE_VECTOR_3D_ADVECTION_VECTOR_3D_NORMAL_VECTORS_

#include "PDEAdvectionSphereBenchmarks_BaseInterface.hpp"
#include <ostream>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/VectorMath.hpp>
#include "PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump.hpp"



class PDEAdvectionSphereBenchmark_advection_vector_3d_normal_vectors	:
	public PDEAdvectionSphereBenchmarks_BaseInterface
{
	PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	PDEAdvectionSphereBenchmark_advection_vector_3d_normal_vectors()
	{
	}

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		return (
				i_benchmark_name == "vector_3d_normal_vector"	||
				false
			);
	}



	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		PDEAdvectionSphereBenchmarks_BaseInterface::shackRegistration(io_shackDict);
		benchmark.shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(benchmark);
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

		benchmark.setup_2_withOps(io_ops);
	}


	std::string printHelp()
	{
		std::ostringstream stream;

		stream << " * Advection test case with 3d vector:" << std::endl;
		stream << "    + 'vector_3d_normal_vector'" << std::endl;

		return stream.str();
	}


	/*
	 * Return number of prognostic fields to be used
	 */
	int getNumPrognosticFields()
	{
		return 3;
	}



	void getInitialState(
		std::vector<sweet::SphereData_Spectral> &o_prognostic_fields,
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 3, "Only a vectorial field (3 elements) supported for this benchmark!");

		const sweet::SphereDataConfig *sphereDataConfig = o_prognostic_fields[0].sphereDataConfig;

		sweet::SphereData_Spectral tmp(sphereDataConfig);
		sweet::SphereData_Spectral vrt(sphereDataConfig);
		sweet::SphereData_Spectral div(sphereDataConfig);
		benchmark.getInitialState_Spectral(tmp, vrt, div);

		ops->vrtdiv_to_uv(vrt, div, o_u, o_v);

		/*
		 * Setup prognostic fields to k vector (perpendicular to point on sphere)
		 */

		sweet::SphereData_Physical x_phys(sphereDataConfig);
		sweet::SphereData_Physical y_phys(sphereDataConfig);
		sweet::SphereData_Physical z_phys(sphereDataConfig);

		x_phys.physical_update_lambda(
				[&](double lon, double lat, double &o_data)
				{
					double ret[3];
					sweet::VectorMath::point_latlon_to_cartesian__scalar(lon, lat, ret[0], ret[1], ret[2]);
					o_data = ret[0];
				}
		);

		y_phys.physical_update_lambda(
				[&](double lon, double lat, double &o_data)
				{
					double ret[3];
					sweet::VectorMath::point_latlon_to_cartesian__scalar(lon, lat, ret[0], ret[1], ret[2]);
					o_data = ret[1];
				}
		);

		z_phys.physical_update_lambda(
				[&](double lon, double lat, double &o_data)
				{
					double ret[3];
					sweet::VectorMath::point_latlon_to_cartesian__scalar(lon, lat, ret[0], ret[1], ret[2]);
					o_data = ret[2];
				}
		);


		o_prognostic_fields[0] = x_phys;
		o_prognostic_fields[1] = y_phys;
		o_prognostic_fields[2] = z_phys;
	}
};

#endif
