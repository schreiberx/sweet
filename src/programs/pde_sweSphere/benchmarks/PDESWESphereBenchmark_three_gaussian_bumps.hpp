/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_THREE_GAUSSIAN_BUMPS_HPP_


#include "PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>



class PDESWESphereBenchmark_three_gaussian_bumps	:
		public PDESWESphereBenchmarks_BaseInterface
{
	sweet::ShackDictionary *shackDict = nullptr;
	sweet::SphereOperators *ops = nullptr;


public:
	PDESWESphereBenchmark_three_gaussian_bumps()
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



	void setup_1_shackData()
	{
		std::cout << "!!! WARNING !!!" << std::endl;
		std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
		std::cout << "!!! WARNING !!!" << std::endl;

		shackPDESWESphere->sphere_rotating_coriolis_omega = 7.292e-5;
		shackPDESWESphere->gravitation = 9.80616;
		shackSphereDataOps->sphere_radius = 6.37122e6;
		shackPDESWESphere->h0 = 29400.0/shackPDESWESphere->gravitation;

#if 0
		// Scale geopotential to make NL influencing the stiffness stronger
		shackPDESWESphere->h0 *= 0.2;
		shackPDESWESphere->gravitation *= 0.2;
#endif
	}

	void setup_2_withOps(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;
	}

	void clear()
	{
	}

	std::string printHelp()
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


	void getInitialState(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{

		o_phi_pert.spectral_set_zero();
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.1, M_PI/3, 20.0)*0.1*shackPDESWESphere->h0;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.6, M_PI/5.0, 80.0)*0.1*shackPDESWESphere->h0;
		o_phi_pert += get_gaussian_bump(2.0*M_PI*0.8, -M_PI/4, 360.0)*0.1*shackPDESWESphere->h0;

		o_vrt.spectral_set_zero();
		o_div.spectral_set_zero();
	}
};

#endif
