/*
 * Author: Francois Hamon, Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_5_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_5_HPP_


#include "PDESWESphereBenchmarks_HelperGeostropicBalance.hpp"
#include "PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>


class PDESWESphereBenchmark_williamson_5_flow_over_mountain	:
		public PDESWESphereBenchmarks_BaseInterface
{
	sweet::ShackDictionary *shackDict = nullptr;
	sweet::SphereOperators *ops = nullptr;

	PDESWESphereBenchmarks_HelperGeostropicBalance helperGeostropicBalance;

public:
	PDESWESphereBenchmark_williamson_5_flow_over_mountain()
	{
	}


	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		PDESWESphereBenchmarks_BaseInterface::shackRegistration(io_shackDict);
		helperGeostropicBalance.shackRegistration(io_shackDict);
		return true;
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson5"	||
				benchmark_name == "flow_over_mountain" ||
				false
		;
	}



	void setup_1_shackData()
	{
		std::cout << "!!! WARNING !!!" << std::endl;
		std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
		std::cout << "!!! WARNING !!!" << std::endl;

		/// Setup Williamson's parameters
		shackPDESWESphere->sphere_rotating_coriolis_omega = 7.292e-5;
		shackPDESWESphere->gravitation = 9.80616;
		shackSphereDataOps->sphere_radius = 6.37122e6;
		shackPDESWESphere->h0 = 5600;
	}

	void setup_2_withOps(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;

		helperGeostropicBalance.setup(ops);

		// initialize the topography
		shackPDESWEBenchmark->h_topography.setup(io_ops->sphereDataConfig);
	}


	void clear()
	{
		shackPDESWEBenchmark->h_topography.clear();
	}


	std::string printHelp()
	{
		std::ostringstream stream;
		std::cout << "  WILLIAMSON #5:" << std::endl;
		std::cout << "     'williamson5'/'flow_over_mountain': Flow over mountain benchmark" << std::endl;
		return stream.str();
	}



	void getInitialState(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{

		double gh0 = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

		const double u0 = 20.0;

		/*
		 * Setup V=0
		 */
		sweet::SphereData_Physical vg(o_phi_pert.sphereDataConfig);
		vg.physical_set_zero();

		/*
		 * Setup U=...
		 * initial velocity along longitude
		 */
		sweet::SphereData_Physical ug(o_phi_pert.sphereDataConfig);
		ug.physical_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				o_data = u0 * std::cos(phi);
			}
		);

		ops->uv_to_vrtdiv(ug, vg, o_vrt, o_div);

		sweet::SphereData_Physical hg(o_phi_pert.sphereDataConfig);

		helperGeostropicBalance.computeGeostrophicBalance_nonlinear(
				o_vrt,
				o_div,
				o_phi_pert
		);

		o_phi_pert -= gh0;
	}


public:
	static
	void setup_topography_(
			sweet::SphereData_Physical &o_h_topo,
			double i_R            = M_PI/9.,
			double i_h_topo_0     = 2000.,
			double i_center_lon   = 3.*M_PI/2.,
			double i_center_lat   = M_PI/6.
	)
	{
		const double center_lat = i_center_lat;
		const double center_lon = i_center_lon;

		o_h_topo.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{

				const double phi1    = asin(mu);
				const double phi2    = center_lat;
				const double lambda1 = lon;
				const double lambda2 = center_lon;

				const double r_squared = std::min( i_R*i_R, (phi1-phi2)*(phi1-phi2) + (lambda1-lambda2)*(lambda1-lambda2) );

				o_data = i_h_topo_0 * ( 1. - sqrt(r_squared) / i_R );
			}
		);
	}

public:
	void setup_topography()
	{
		if (shackDict == nullptr)
			SWEETError("Benchmarks are not yet initialized");

		if (benchmark_name == "flow_over_mountain")
		{
			// set the topography flag to true
			shackPDESWEBenchmark->use_topography = true;

			// setup the parameters for the flow-over-mountain test case
			const double R			= M_PI/9.;
			const double h_topo_0	 = 2000.;
			const double i_center_lon = 3.*M_PI/2.;
			const double i_center_lat = M_PI/6.;

			shackPDESWEBenchmark->h_topography.physical_set_zero();

			// setup the topography vector
			setup_topography_(
					shackPDESWEBenchmark->h_topography,
					R,
					h_topo_0,
					i_center_lon,
					i_center_lat
			);
		}
		else
		{
			// set the topography flag to false
			shackPDESWEBenchmark->use_topography = false;
		}
	}



};

#endif
