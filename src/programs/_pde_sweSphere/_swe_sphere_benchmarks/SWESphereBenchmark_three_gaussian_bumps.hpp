/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP_


#include "SWESphereBenchmarks_helpers.hpp"
#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>



class SWESphereBenchmark_three_gaussian_bumps	: public SWESphereBenchmarks_interface
{
	sweet::ShackDictionary *shackDict = nullptr;
	sweet::SphereOperators *ops = nullptr;


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
		stream << "  'three_gaussian_bumps':	Three Gaussian bumps on the geopotential field." << std::endl;
		return stream.str();
	}

public:
	sweet::SphereData_Physical get_gaussian_bump(
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{
		sweet::SphereData_Physical o_h(ops->sphereDataConfig);

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
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{
		if (shackDict->benchmark.benchmark_override_simvars)
		{
			if (shackDict->timecontrol.current_simulation_time == 0 && 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;

				shackDict->sim.sphere_rotating_coriolis_omega = 7.292e-5;
				shackDict->sim.gravitation = 9.80616;
				shackDict->sim.sphere_radius = 6.37122e6;
				shackDict->sim.h0 = 29400.0/shackDict->sim.gravitation;

#if 0
				// Scale geopotential to make NL influencing the stiffness stronger
				shackDict->sim.h0 *= 0.2;
				shackDict->sim.gravitation *= 0.2;
#endif

				ops->setup(ops->sphereDataConfig, &(shackDict->sim));
			}
		}

		o_phi_pert.spectral_set_zero();
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.1, M_PI/3, 20.0)*0.1*shackDict->sim.h0;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.6, M_PI/5.0, 80.0)*0.1*shackDict->sim.h0;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.8, -M_PI/4, 360.0)*0.1*shackDict->sim.h0;

		o_vrt.spectral_set_zero();
		o_div.spectral_set_zero();
	}
};

#endif
