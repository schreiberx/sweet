/*
 * Author: Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_VECTOR_3D_ADVECTION_VECTOR_3D_NORMAL_VECTORS_
#define SRC_BENCHMARKS_SPHERE_VECTOR_3D_ADVECTION_VECTOR_3D_NORMAL_VECTORS_

#include "BenchmarksSphereAdvection_interface.hpp"
#include <ostream>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SWEETVectorMath.hpp>
#include "../swe_sphere_benchmarks/SWESphereBenchmark_williamson_1_advection_gauss_bump.hpp"



class BenchmarksSphereAdvection_advection_vector_3d_normal_vectors	: public BenchmarksSphereAdvection_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;

	SWESphereBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	BenchmarksSphereAdvection_advection_vector_3d_normal_vectors()
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

		stream << " * Advection test case with 3d vector:" << std::endl;
		stream << "    + 'vector_3d_normal_vector'" << std::endl;

		return stream.str();
	}


	/*
	 * Return number of prognostic fields to be used
	 */
	int get_num_prognostic_fields()
	{
		return 3;
	}



	void get_initial_state(
		std::vector<SphereData_Spectral*> &o_prognostic_fields,
		SphereData_Physical &o_u,
		SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 3, "Only a vectorial field (3 elements) supported for this benchmark!");

		const SphereData_Config *sphereDataConfig = o_prognostic_fields[0]->sphereDataConfig;

		SphereData_Spectral tmp(sphereDataConfig);
		SphereData_Spectral vrt(sphereDataConfig);
		SphereData_Spectral div(sphereDataConfig);
		benchmark.get_initial_state(tmp, vrt, div);

		ops->vrtdiv_to_uv(vrt, div, o_u, o_v);

		/*
		 * Setup prognostic fields to k vector (perpendicular to point on sphere)
		 */

		SphereData_Physical x_phys(sphereDataConfig);
		SphereData_Physical y_phys(sphereDataConfig);
		SphereData_Physical z_phys(sphereDataConfig);

		x_phys.physical_update_lambda(
				[&](double lon, double lat, double &o_data)
				{
					double ret[3];
					SWEETVectorMath::point_latlon_to_cartesian__scalar(lon, lat, ret[0], ret[1], ret[2]);
					o_data = ret[0];
				}
		);

		y_phys.physical_update_lambda(
				[&](double lon, double lat, double &o_data)
				{
					double ret[3];
					SWEETVectorMath::point_latlon_to_cartesian__scalar(lon, lat, ret[0], ret[1], ret[2]);
					o_data = ret[1];
				}
		);

		z_phys.physical_update_lambda(
				[&](double lon, double lat, double &o_data)
				{
					double ret[3];
					SWEETVectorMath::point_latlon_to_cartesian__scalar(lon, lat, ret[0], ret[1], ret[2]);
					o_data = ret[2];
				}
		);


		*o_prognostic_fields[0] = x_phys;
		*o_prognostic_fields[1] = y_phys;
		*o_prognostic_fields[2] = z_phys;
	}
};

#endif
