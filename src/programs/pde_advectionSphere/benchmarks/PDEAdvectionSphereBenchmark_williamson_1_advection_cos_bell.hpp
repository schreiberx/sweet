/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_1_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_1_HPP_

#include "PDEAdvectionSphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>



class PDEAdvectionSphereBenchmark_williamson_1_advection_cos_bell	:
	public PDEAdvectionSphereBenchmarks_BaseInterface
{
	double _default_sphere_radius = 6.37122e6;

public:
	PDEAdvectionSphereBenchmark_williamson_1_advection_cos_bell()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson1"		||
				benchmark_name == "adv_cosine_bell"	||
				benchmark_name == "advection_cosine_bell"	||
				false
		;
	}


	void setup_1_shackData()
	{
		std::cout << "!!! WARNING !!!" << std::endl;
		std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
		std::cout << "!!! WARNING !!!" << std::endl;

		shackSphereDataOps->sphere_radius = _default_sphere_radius;
	}

	void setup_2_withOps(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;
	}


	std::string printHelp()
	{
		std::ostringstream stream;
		stream << "  WILLIAMSON #1:" << std::endl;
		stream << "     'williamson1'" << std::endl;
		stream << "     'advection_cosine_bell'" << std::endl;
		stream << "     'adv_cosine_bell': Advection test case of cosine bell" << std::endl;
		stream << "         OPTION:" << std::endl;
		stream << "         --advection-rotation-angle=[angle]" << std::endl;
		return stream.str();
	}

	void getInitialState(
		std::vector<sweet::SphereData_Spectral> &o_prognostic_fields,
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	)
	{
		sweet::SphereData_Spectral vrt(ops->sphereDataConfig);
		sweet::SphereData_Spectral div(ops->sphereDataConfig);

		getInitialState_Spectral(o_prognostic_fields[0], vrt, div);
		ops->vrtdiv_to_uv(vrt, div, o_u, o_v);
	}

	void getInitialState_Spectral(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{
		/*
		 * Advection test case
		 * See Williamson test case, eq. (77), (78), (79)
		 */

		/*
		 * Coefficients from benchmark
		 */
		double gravitation = 9.80616;
		double h0 = 1000.0;
		double gh0 = gravitation * h0;
		
		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0.0;

		//theta_c = M_PI*0.5*0.8;
		double a = shackSphereDataOps->sphere_radius;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		sweet::SphereData_Physical phi_phys(o_phi_pert.sphereDataConfig);

		phi_phys.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				double r = a * std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				if (r < R)
					io_data = gh0/2.0*(1.0+std::cos(M_PI*r/R));
				else
					io_data = 0;
			}
		);

		o_phi_pert.loadSphereDataPhysical(phi_phys);


		/*
		 * Both versions are working
		 */
#if 1
		sweet::SphereData_Physical stream_function(o_phi_pert.sphereDataConfig);

		stream_function.physical_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_theta = i_lat;
				double i_lambda = i_lon;
				double alpha = shackPDEAdvBenchmark->sphere_advection_rotation_angle;

				io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
			}
		);

		o_vrt = ops->laplace(stream_function);

#else

		o_vrt.physical_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_theta = i_lat;
				double i_lambda = i_lon;
				double alpha = shackDict->benchmark.sphere_advection_rotation_angle;

				io_data = 2.0*u0/a*(-std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha) + std::sin(i_theta)*std::cos(alpha));
			}
		);
#endif
		o_div.spectral_set_zero();
	}
};

#endif
