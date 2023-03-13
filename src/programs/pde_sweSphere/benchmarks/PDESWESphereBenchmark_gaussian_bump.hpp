/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_GAUSSIAN_BUMP_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_GAUSSIAN_BUMP_HPP_


#include "PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>



class PDESWESphereBenchmark_gaussian_bump	:
		public PDESWESphereBenchmarks_BaseInterface
{
public:
	PDESWESphereBenchmark_gaussian_bump()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				i_benchmark_name == "gaussian_bump"			||
				i_benchmark_name == "gaussian_bump_phi"			||
				i_benchmark_name == "gaussian_bump_phi_pint"		||
				i_benchmark_name == "sharp_gaussian_bump" ||
				i_benchmark_name == "gaussian_bump_vrt"			||
				i_benchmark_name == "gaussian_bump_div"			||
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

		if (benchmark_name == "gaussian_bump_phi_pint")
		{
			shackPDESWESphere->h0 = 29400.0;
		}
		else
		{
			shackPDESWESphere->h0 = 29400.0/shackPDESWESphere->gravitation;
			// Scale geopotential to make NL influencing the stiffness stronger
			//shackPDESWESphere->h0 *= 0.2;
			//shackPDESWESphere->gravitation *= 0.2;
		}
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
		stream << "  'gaussian_bumps_pvd':	Gaussian bump slightly delocated on pot, vrt and div field" << std::endl;
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


		/*
		 * We scale the variables so that cancellation errors are reduced.
		 *
		 * Div is related to geopotential via the laplace operator which requires a scaling of 1/r^2
		 *
		 * Vrt and div need to be of the same order of magnitude due to the stream formulation
		 */
		double phi_scale;
		if (benchmark_name == "gaussian_bump_phi_pint")
			phi_scale = 6000 * shackPDESWESphere->gravitation;
		else
			phi_scale = 0.1*shackPDESWESphere->h0;
		double vrt_scale = phi_scale/(shackSphereDataOps->sphere_radius*shackSphereDataOps->sphere_radius);
		double div_scale = phi_scale/(shackSphereDataOps->sphere_radius*shackSphereDataOps->sphere_radius);

		if (benchmark_name == "gaussian_bump" || benchmark_name == "gaussian_bump_phi")
			o_phi_pert = get_gaussian_bump(M_PI, M_PI/3, 20.0)*phi_scale;
		else if (benchmark_name == "gaussian_bump_phi_pint")
			o_phi_pert = get_gaussian_bump(M_PI, M_PI/4., 20.)*phi_scale;
		else if (benchmark_name == "sharp_gaussian_bump")
			o_phi_pert = get_gaussian_bump(M_PI, M_PI/4, 40.0)*6000;
		else
			o_phi_pert.spectral_set_zero();

		if (benchmark_name == "gaussian_bump_vrt")
			o_vrt = get_gaussian_bump(M_PI + M_PI*0.1, M_PI/3 + M_PI*0.05, 20.0)*vrt_scale;
		else
			o_vrt.spectral_set_zero();

		if (benchmark_name == "gaussian_bump_div")
			o_div = get_gaussian_bump(M_PI + M_PI*0.05, M_PI/3 + M_PI*0.1, 20.0)*div_scale;
		else
			o_div.spectral_set_zero();
	}
};

#endif
