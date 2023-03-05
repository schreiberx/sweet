/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_1_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_1_HPP_

#include "SWESphereBenchmarks_helpers.hpp"
#include "SWESphereBenchmarks_interface.hpp"
#include <ostream>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>



class SWESphereBenchmark_williamson_1_advection_cos_bell	: public SWESphereBenchmarks_interface
{
	sweet::ShackDictionary *shackDict = nullptr;
	sweet::SphereOperators *ops = nullptr;


public:
	SWESphereBenchmark_williamson_1_advection_cos_bell()
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


	void setup(
			sweet::ShackDictionary *i_shackDict,
			sweet::SphereOperators *i_ops
	)
	{
		shackDict = i_shackDict;
		ops = i_ops;
	}


	std::string get_help()
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


	void get_initial_state(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{
		double gh0 = shackDict->sim.gravitation * shackDict->sim.h0;

		/*
		 * Advection test case
		 * See Williamson test case, eq. (77), (78), (79)
		 */

		if (shackDict->benchmark.benchmark_override_simvars)
		{
			if (shackDict->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			shackDict->sim.sphere_rotating_coriolis_omega = 7.292e-5;
			shackDict->sim.gravitation = 9.80616;
			shackDict->sim.sphere_radius = 6.37122e6;
			shackDict->sim.h0 = 1000.0;

			// reset operator
			ops->setup(o_phi_pert.sphereDataConfig, &(shackDict->sim));
		}

		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0.0;
		double a = shackDict->sim.sphere_radius;

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
					io_data = shackDict->sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
				else
					io_data = 0;

				io_data *= shackDict->sim.gravitation;
			}
		);

		o_phi_pert.loadSphereDataPhysical(phi_phys);

		sweet::SphereData_Physical stream_function(o_phi_pert.sphereDataConfig);

		stream_function.physical_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_theta = i_lat;
				double i_lambda = i_lon;
				double alpha = shackDict->benchmark.sphere_advection_rotation_angle;

				io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
			}
		);

		o_vrt = ops->laplace(stream_function);
		o_div.spectral_set_zero();

		o_phi_pert -= gh0;
	}
};

#endif
