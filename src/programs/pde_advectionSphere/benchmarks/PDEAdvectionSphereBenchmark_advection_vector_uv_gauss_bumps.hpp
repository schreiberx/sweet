/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_GAUSS_BUMPS_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_VECTOR_UV_GAUSS_BUMPS_HPP_

#include "PDEAdvectionSphereBenchmarks_BaseInterface.hpp"
#include <ostream>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/VectorMath.hpp>
#include "PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump.hpp"



class PDEAdvectionSphereBenchmark_advection_vector_uv_gauss_bumps	:
	public PDEAdvectionSphereBenchmarks_BaseInterface
{
	PDEAdvectionSphereBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	PDEAdvectionSphereBenchmark_advection_vector_uv_gauss_bumps()
	{
	}



	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		return (
				i_benchmark_name == "advection_vector_uv_gauss_bumps"	||
				i_benchmark_name == "advection_vector_uv_gaussian_bumps"	||
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

		stream << " * Advection test case with 2d vector in lat-lon space, setup with gaussian bumps:" << std::endl;
		stream << "    + 'advection_vector_uv_gauss_bumps'" << std::endl;

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
		std::vector<sweet::SphereData_Spectral> &o_prog_vec,
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	)
	{
		sweet::SphereData_Spectral tmp(ops->sphereDataConfig);
		sweet::SphereData_Spectral vrt(ops->sphereDataConfig);
		sweet::SphereData_Spectral div(ops->sphereDataConfig);
		benchmark.getInitialState_Spectral(tmp, vrt, div);

		// Convert vrt/div to velocity field
		ops->vrtdiv_to_uv(vrt, div, o_u, o_v);

		/*
		 * IMPORTANT INFORMATION:
		 * The prognostic fields here are the vrt/div of the velocities!
		 * So we first convert the bumps in u/v field to vrt/div fields.
		 * Be careful about pole problems! The bump needs to be numerically 0 at the poles!!!
		 */


		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0.0;
		double i_exp_fac = 20.0;

		sweet::SphereData_Physical u_phys(ops->sphereDataConfig);
		u_phys.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				double d = std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				io_data = std::exp(-d*d*i_exp_fac);
			}
		);

		sweet::SphereData_Physical v_phys(ops->sphereDataConfig);
		v_phys.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				double d = std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				io_data = std::exp(-d*d*i_exp_fac);

			}
		);

		ops->uv_to_vrtdiv(
				u_phys, v_phys,
				o_prog_vec[0], o_prog_vec[1]
			);
	}
};

#endif
