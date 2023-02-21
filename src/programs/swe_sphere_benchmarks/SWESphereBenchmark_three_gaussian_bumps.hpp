/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP_


#include "SWESphereBenchmarks_helpers.hpp"
#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>



class SWESphereBenchmark_three_gaussian_bumps	: public SWESphereBenchmarks_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


public:
	SWESphereBenchmark_three_gaussian_bumps()
	{
	}


	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
			i_benchmark_name == "three_gaussian_bumps"	||
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
		stream << "  'three_gaussian_bumps':	Three Gaussian bumps on the geopotential field." << std::endl;
		return stream.str();
	}

public:
	SphereData_Physical get_gaussian_bump(
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{
		SphereData_Physical o_h(ops->sphereDataConfig);

		o_h.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					// https://en.wikipedia.org/wiki/Great-circle_distance
					// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
					// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

					double phi1 = asin(mu);
					double phi2 = i_center_lat;
					double lambda1 = lon;
					double lambda2 = i_center_lon;

					double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

					o_data = std::exp(-d*d*i_exp_fac);
				}
		);

		return o_h;
	}


	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{
		if (simVars->benchmark.benchmark_override_simvars)
		{
			if (simVars->timecontrol.current_simulation_time == 0 && 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;

				simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.sphere_radius = 6.37122e6;
				simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

#if 0
				// Scale geopotential to make NL influencing the stiffness stronger
				simVars->sim.h0 *= 0.2;
				simVars->sim.gravitation *= 0.2;
#endif

				ops->setup(ops->sphereDataConfig, &(simVars->sim));

				simVars->sim.outputConfig();
			}
		}

		o_phi_pert.spectral_set_zero();
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.1, M_PI/3, 20.0)*0.1*simVars->sim.h0;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.6, M_PI/5.0, 80.0)*0.1*simVars->sim.h0;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.8, -M_PI/4, 360.0)*0.1*simVars->sim.h0;

		o_vrt.spectral_set_zero();
		o_div.spectral_set_zero();
	}
};

#endif
