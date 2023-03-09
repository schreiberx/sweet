/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_2_LINEAR_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_WILLIAMSON_2_LINEAR_HPP_


#include "PDESWESphereBenchmarks_HelperGeostropicBalance.hpp"
#include "PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>



class PDESWESphereBenchmark_williamson_2_geostrophic_balance_linear	:
		public PDESWESphereBenchmarks_BaseInterface
{
	sweet::ShackDictionary *shackDict = nullptr;
	sweet::SphereOperators *ops = nullptr;

	sweet::SphereData_Physical fg;

public:
	PDESWESphereBenchmark_williamson_2_geostrophic_balance_linear()
	{

	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson2_linear"			||
				benchmark_name == "geostrophic_balance_linear"	||
				benchmark_name == "geostrophic_balance_linear_2"	||
				benchmark_name == "geostrophic_balance_linear_4"	||
				benchmark_name == "geostrophic_balance_linear_8"	||
				benchmark_name == "geostrophic_balance_linear_16"	||
				benchmark_name == "geostrophic_balance_linear_32"	||
				benchmark_name == "geostrophic_balance_linear_64"	||
				benchmark_name == "geostrophic_balance_linear_128"	||
				benchmark_name == "geostrophic_balance_linear_256"	||
				benchmark_name == "geostrophic_balance_linear_512"	||
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

		if (shackPDESWESphere->sphere_use_fsphere)
			fg = ops->getFG_fSphere(shackPDESWESphere->sphere_fsphere_f0);
		else
			fg = ops->getFG_rotatingSphere(shackPDESWESphere->sphere_rotating_coriolis_omega);
	}


	void clear()
	{
	}


	std::string printHelp()
	{
		std::ostringstream stream;
		stream << "  WILLIAMSON #2 (variant):" << std::endl;
		stream << "     'williamson2_linear'" << std::endl;
		stream << "     'geostrophic_balance_linear': Geostrophic balance for linear SWE" << std::endl;
		return stream.str();
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_linear(
			const sweet::SphereOperators *i_ops,
			const sweet::SphereData_Spectral &i_vort,
			const sweet::SphereData_Spectral &i_div,
			sweet::SphereData_Spectral &o_phi
	)
	{
		const sweet::SphereData_Config *sphereDataConfig = o_phi.sphereDataConfig;

		/*
		 * Compute vorticity and divergence from velocities
		 */
		sweet::SphereData_Physical u(sphereDataConfig);
		sweet::SphereData_Physical v(sphereDataConfig);

		i_ops->vrtdiv_to_uv(i_vort, i_div, u, v);

		sweet::SphereData_Physical vrtg = i_vort.toPhys();

		using namespace sweet;
		sweet::SphereData_Physical tmpg1 = u*fg;
		sweet::SphereData_Physical tmpg2 = v*fg;

		sweet::SphereData_Spectral tmpspec1(sphereDataConfig);
		sweet::SphereData_Spectral tmpspec2(sphereDataConfig);

		i_ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		o_phi = i_ops->inv_laplace(tmpspec1);
	}



	void getInitialState(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{

		/*
		 * Linear-only!!!
		 *
		 * The original version is for the non-linear case
		 */

		double a = shackSphereDataOps->sphere_radius;
		//double omega = shackPDESWESphere->sphere_rotating_coriolis_omega;
		double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		double alpha = 0;

		double freq_multiplier = 1.0;

		if (benchmark_name == "geostrophic_balance_linear_2")
			freq_multiplier = 2.0;
		else if (benchmark_name == "geostrophic_balance_linear_4")
			freq_multiplier = 4.0;
		else if (benchmark_name == "geostrophic_balance_linear_8")
			freq_multiplier = 8.0;
		else if (benchmark_name == "geostrophic_balance_linear_16")
			freq_multiplier = 16.0;
		else if (benchmark_name == "geostrophic_balance_linear_32")
			freq_multiplier = 32.0;
		else if (benchmark_name == "geostrophic_balance_linear_64")
			freq_multiplier = 64.0;
		else if (benchmark_name == "geostrophic_balance_linear_128")
			freq_multiplier = 128.0;
		else if (benchmark_name == "geostrophic_balance_linear_256")
			freq_multiplier = 256.0;
		else if (benchmark_name == "geostrophic_balance_linear_512")
			freq_multiplier = 512.0;


		/*
		 * Setup U
		 */
		sweet::SphereData_Physical ug(ops->sphereDataConfig);
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
		sweet::SphereData_Physical vg(ops->sphereDataConfig);
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
			sweet::SphereData_Physical vortg(ops->sphereDataConfig);
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


		std::cout << "[MULE] geostrophic_balance_analytical_setup: 1" << std::endl;

		computeGeostrophicBalance_linear(
				ops,
				o_vrt,
				o_div,
				o_phi_pert
		);

		/*
		 * TODO: Maybe we need to add
		 */
	}
};

#endif
