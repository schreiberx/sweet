/*
 * Author: Francois Hamon, Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_6_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_6_HPP_


#include <benchmarks_sphere/SWESphereBenchmarks_helpers.hpp>
#include <benchmarks_sphere/SWESphereBenchmarks_interface.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereData_Config.hpp>

class SWESphereBenchmark_williamson_6_rossby_haurwitz_wave	: public SWESphereBenchmarks_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;

public:
	SWESphereBenchmark_williamson_6_rossby_haurwitz_wave()
	{
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
		std::cout << "  WILLIAMSON #6:" << std::endl;
		std::cout << "     'williamson6'/'rossby_haurwitz_wave': Rossby Haurwitz wave" << std::endl;
		return stream.str();
	}



	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{

		if (simVars->timecontrol.current_simulation_time == 0)
		{
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
		}

		/// Setup Williamson's parameters
		simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
		simVars->sim.gravitation = 9.80616;
		simVars->sim.sphere_radius = 6.37122e6;
		simVars->sim.h0 = 8000;

		double gh0 = simVars->sim.gravitation*simVars->sim.h0;

		// update operator because we changed the simulation parameters
		ops->setup(o_phi_pert.sphereDataConfig, &(simVars->sim));

		const double omega = 7.484e-6;
		const double K = omega;
		const int R	= 4;
		const double a = simVars->sim.sphere_radius;

		/*
		 * Setup U=...
		 */
		SphereData_Physical ug(o_phi_pert.sphereDataConfig);
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
		SphereData_Physical vg(o_phi_pert.sphereDataConfig);
		vg.physical_update_lambda(
						  [&](double lon, double phi, double &o_data)
								{
									o_data = - a * K * R * pow(cos(phi), R-1) * sin(phi) * sin(R*lon);
								}
						  );

		ops->uv_to_vortdiv(ug, vg, o_vrt, o_div);

		SphereData_Physical hg(o_phi_pert.sphereDataConfig);
		SWESphereBenchmark_williamson_2_geostrophic_balance::computeGeostrophicBalance_nonlinear(
				ops,
				o_vrt,
				o_div,
				o_phi_pert
		);

		o_phi_pert -= gh0;
	}

};

#endif
