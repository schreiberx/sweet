/*
 * Author: Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_1_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_1_HPP_

#include "SWESphereBenchmarks_helpers.hpp"
#include "SWESphereBenchmarks_interface.hpp"
#include <ostream>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereData_Config.hpp>



class SWESphereBenchmark_williamson_1_advection_cos_bell	: public SWESphereBenchmarks_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


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
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
	)
	{
		simVars = i_simVars;
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
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{
		double gh0 = simVars->sim.gravitation * simVars->sim.h0;

		/*
		 * Advection test case
		 * See Williamson test case, eq. (77), (78), (79)
		 */

		if (simVars->benchmark.benchmark_override_simvars)
		{
			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
			simVars->sim.gravitation = 9.80616;
			simVars->sim.sphere_radius = 6.37122e6;
			simVars->sim.h0 = 1000.0;

			// reset operator
			ops->setup(o_phi_pert.sphereDataConfig, &(simVars->sim));
		}

		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0.0;
		double a = simVars->sim.sphere_radius;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		SphereData_Physical phi_phys(o_phi_pert.sphereDataConfig);

		phi_phys.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				double r = a * std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				if (r < R)
					io_data = simVars->sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
				else
					io_data = 0;

				io_data *= simVars->sim.gravitation;
			}
		);

		o_phi_pert.loadSphereDataPhysical(phi_phys);

		SphereData_Physical stream_function(o_phi_pert.sphereDataConfig);

		stream_function.physical_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_theta = i_lat;
				double i_lambda = i_lon;
				double alpha = simVars->benchmark.sphere_advection_rotation_angle;

				io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
			}
		);

		o_vrt = ops->laplace(stream_function);
		o_div.spectral_set_zero();

		o_phi_pert -= gh0;
	}
};

#endif
