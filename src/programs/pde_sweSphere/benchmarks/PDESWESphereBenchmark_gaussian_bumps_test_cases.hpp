/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_GAUSSIAN_BUMPS_TEST_CASES_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_GAUSSIAN_BUMPS_TEST_CASES_HPP_


#include "PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>



class PDESWESphereBenchmark_gaussian_bumps_test_cases	:
		public PDESWESphereBenchmarks_BaseInterface
{
public:
	PDESWESphereBenchmark_gaussian_bumps_test_cases()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
			i_benchmark_name == "gaussian_bumps_test_cases"	||
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
		stream << "  'gaussian_bumps_test_cases':	Gaussian bumps, special test case" << std::endl;
		return stream.str();
	}


public:
	void setup_gaussian_bump(
			sweet::SphereData_Physical &o_data,
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{
		o_data.physical_update_lambda_gaussian_grid(
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
	}


	void getInitialState(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{
		sweet::SphereData_Physical tmp(o_phi_pert.sphereDataConfig);

		setup_gaussian_bump(tmp, 2.0*M_PI*0.1, M_PI/3, 1.0);
		o_phi_pert.loadSphereDataPhysical(tmp);
		o_phi_pert *= 0.1;
		o_phi_pert += shackPDESWESphere->h0*shackPDESWESphere->gravitation;

		setup_gaussian_bump(tmp, 2.0*M_PI*0.1, M_PI/3, 1.0);
		o_vrt.loadSphereDataPhysical(tmp);
		o_vrt *= -1e-8;
		//o_vort *= 0;
		setup_gaussian_bump(tmp, 2.0*M_PI*0.1, M_PI/3, 1.0);
		o_div.loadSphereDataPhysical(tmp);
		o_div *= 1e-8;

		/*
		 * Convert forward/backward to velocity space to apply a certain truncation
		 */
		sweet::SphereData_Physical ug(o_phi_pert.sphereDataConfig);
		sweet::SphereData_Physical vg(o_phi_pert.sphereDataConfig);
		ops->vrtdiv_to_uv(o_vrt, o_div, ug, vg);
		ops->uv_to_vrtdiv(ug, vg, o_vrt, o_div);
	}
};

#endif
