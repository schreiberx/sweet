/*
 * Author: Francois Hamon, Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_6_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_6_HPP_


#include "PDESWESphereBenchmarks_HelperGeostropicBalance.hpp"
#include "PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>


class PDESWESphereBenchmark_williamson_6_rossby_haurwitz_wave	:
		public PDESWESphereBenchmarks_BaseInterface
{
	PDESWESphereBenchmarks_HelperGeostropicBalance helperGeostropicBalance;

public:
	PDESWESphereBenchmark_williamson_6_rossby_haurwitz_wave()
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
				benchmark_name == "williamson6"	||
				benchmark_name == "rossby_haurwitz_wave" ||
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
		shackPDESWESphere->h0 = 8000;
	}

	void setup_2_withOps(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;

		helperGeostropicBalance.setup(ops);
	}


	void clear()
	{
	}



	std::string printHelp()
	{
		std::ostringstream stream;
		std::cout << "  WILLIAMSON #6:" << std::endl;
		std::cout << "     'williamson6'/'rossby_haurwitz_wave': Rossby Haurwitz wave" << std::endl;
		return stream.str();
	}



	void getInitialState(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{
		double gh0 = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

		const double omega = 7.484e-6;
		const double K = omega;
		const int R	= 4;
		const double a = shackSphereDataOps->sphere_radius;

		/*
		 * Setup U=...
		 */
		sweet::SphereData_Physical ug(o_phi_pert.sphereDataConfig);
		ug.physical_update_lambda(
						[&](double lon, double phi, double &o_data)
							{
								o_data = a * omega * cos(phi)
										+ a * K * pow(cos(phi), R-1) * (R * sin(phi)*sin(phi) - cos(phi)*cos(phi)) * cos(R*lon);
							}
						);


		/*
		 * Setup V=...
		 */
		sweet::SphereData_Physical vg(o_phi_pert.sphereDataConfig);
		vg.physical_update_lambda(
						  [&](double lon, double phi, double &o_data)
								{
									o_data = - a * K * R * pow(cos(phi), R-1) * sin(phi) * sin(R*lon);
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

};

#endif
