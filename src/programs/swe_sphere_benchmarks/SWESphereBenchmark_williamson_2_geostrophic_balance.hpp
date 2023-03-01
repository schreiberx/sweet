/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_2_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_2_HPP_


#include "SWESphereBenchmarks_helpers.hpp"
#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereData_Config.hpp>



class SWESphereBenchmark_williamson_2_geostrophic_balance	: public SWESphereBenchmarks_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


public:
	SWESphereBenchmark_williamson_2_geostrophic_balance()
	{
	}


	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson2"			||
				benchmark_name == "geostrophic_balance"	||
				benchmark_name == "geostrophic_balance_1"	||
				benchmark_name == "geostrophic_balance_2"	||
				benchmark_name == "geostrophic_balance_4"	||
				benchmark_name == "geostrophic_balance_8"	||
				benchmark_name == "geostrophic_balance_16"	||
				benchmark_name == "geostrophic_balance_32"	||
				benchmark_name == "geostrophic_balance_64"	||
				benchmark_name == "geostrophic_balance_128"	||
				benchmark_name == "geostrophic_balance_256"	||
				benchmark_name == "geostrophic_balance_512"	||
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
		stream << "  WILLIAMSON #2:" << std::endl;
		stream << "     'williamson2'" << std::endl;
		stream << "     'geostrophic_balance': Geostrophic balance, one wave (standard)" << std::endl;
		stream << "     'geostrophic_balance_[N]': Geostrophic balance, with [N] waves" << std::endl;
		return stream.str();
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	static void computeGeostrophicBalance_nonlinear(
			const SphereOperators_SphereData *i_ops,
			const SphereData_Spectral &i_vort,
			const SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi
	)
	{
		const SphereData_Config *sphereDataConfig = i_ops->sphereDataConfig;

		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereData_Physical ug(sphereDataConfig);
		SphereData_Physical vg(sphereDataConfig);

		i_ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

		SphereData_Physical vrtg = i_vort.toPhys();

		SphereData_Physical tmpg1 = ug*(vrtg+i_ops->fg);
		SphereData_Physical tmpg2 = vg*(vrtg+i_ops->fg);

		SphereData_Spectral tmpspec1(sphereDataConfig);
		SphereData_Spectral tmpspec2(sphereDataConfig);

		i_ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		o_phi = i_ops->inv_laplace(tmpspec1) - 0.5*(ug*ug+vg*vg);
	}




	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{

		/*
		 * geostrophic_balance / geostrophic_balance_1:
		 * Williamson test case 2 for geostrophic balance.
		 *
		 * WARNING: This test uses a balanced solution for the full non-linear equations
		 * See Williamson paper for accurate setup
		 *
		 * "geostrophic_balance_N" means that N is the multiplier for the frequency
		 * in the direction of the Latitude
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
			simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

			ops->setup(ops->sphereDataConfig, &(simVars->sim));
		}

		double a = simVars->sim.sphere_radius;
		double omega = simVars->sim.sphere_rotating_coriolis_omega;
		double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		double alpha = 0;


		double freq_multiplier = 1.0;

		if (benchmark_name == "geostrophic_balance_2")
			freq_multiplier = 2.0;
		else if (benchmark_name == "geostrophic_balance_4")
			freq_multiplier = 4.0;
		else if (benchmark_name == "geostrophic_balance_8")
			freq_multiplier = 8.0;
		else if (benchmark_name == "geostrophic_balance_16")
			freq_multiplier = 16.0;
		else if (benchmark_name == "geostrophic_balance_32")
			freq_multiplier = 32.0;
		else if (benchmark_name == "geostrophic_balance_64")
			freq_multiplier = 64.0;
		else if (benchmark_name == "geostrophic_balance_128")
			freq_multiplier = 128.0;
		else if (benchmark_name == "geostrophic_balance_256")
			freq_multiplier = 256.0;
		else if (benchmark_name == "geostrophic_balance_512")
			freq_multiplier = 512.0;


		/*
		 * Setup U
		 */
		SphereData_Physical ug(ops->sphereDataConfig);
		ug.physical_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				// Eq. 90, Williamson TC paper
				o_data = u0*(std::cos(phi*freq_multiplier)*std::cos(alpha) + std::cos(lon)*std::sin(phi*freq_multiplier)*std::sin(alpha));
			}
		);

		/*
		 * Setup V
		 */
		SphereData_Physical vg(ops->sphereDataConfig);
		vg.physical_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				// Eq. 91, Williamson TC paper
				o_data = -u0*std::sin(lon*freq_multiplier)*std::sin(alpha);
			}
		);

		ops->uv_to_vrtdiv(ug, vg, o_vrt, o_div);



		/**
		 * TEST for non-divergent test case
		 */
		double div_zero = o_div.toPhys().physical_reduce_max_abs();
		if (div_zero > 1e-12)
		{

			std::cout << "Divergence: " << div_zero << std::endl;
			SWEETError("Divergence should be close to 0, maybe there are some numerical round-off errors?");
		}

		if (freq_multiplier == 1.0)
		{
			/**
			 * TEST for correct vorticity
			 */

			/*
			 * Setup relative vorticity
			 */
			SphereData_Physical vortg(ops->sphereDataConfig);
			vortg.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 94, Williamson TC paper

					// absolute vorticity, but we like the relative one
					//o_data = (2.0*u0/a + 2.0*omega)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));

					// relative vorticity
					o_data = (2.0*u0/a)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));
				}
			);

			double vort_diff = (o_vrt.toPhys() - vortg).physical_reduce_max_abs();
			if (vort_diff > 1e-12)
			{
				std::cout << "Vorticity difference: " << vort_diff << std::endl;
				SWEETError("Vorticity fields differ (should be close to 0), maybe there are some numerical round-off errors?");
			}
		}


		bool use_analytical_geostrophic_setup = simVars->misc.comma_separated_tags.find("geostrophic_balance_analytical_setup") != std::string::npos;


		if (use_analytical_geostrophic_setup)
		{
			std::cout << "[MULE] geostrophic_balance_analytical_setup: 1" << std::endl;

			computeGeostrophicBalance_nonlinear(
					ops,
					o_vrt,
					o_div,
					o_phi_pert
			);
		}
		else
		{
			std::cout << "[MULE] geostrophic_balance_analytical_setup: 0" << std::endl;

			// Squared term in Eq. 95, Williamson TC paper
			SphereData_Physical r2(ops->sphereDataConfig);
			r2.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					o_data = -std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha);
					o_data = o_data*o_data;
				}
			);

			// Eq. 95, Williamson TC paper
			SphereData_Physical phig = -(a*omega*u0 + u0*u0/2.0)*r2;

			o_phi_pert.loadSphereDataPhysical(phig);
		}
	}
};

#endif
