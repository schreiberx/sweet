/*
 * lagrangian_test.cpp
 *
 *  Created on: 5 Dec 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "../include/sweet/plane/PlaneData_Physical.hpp"
#include "../include/sweet/plane/PlaneData_Spectral.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include "../programs/swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/Convert_PlaneDataPhysical_to_ScalarDataArray.hpp>
#include <unistd.h>
#include <stdio.h>

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

// Plane data config
PlaneDataConfig planeDataConfigInstance3;
PlaneDataConfig *planeDataConfig3 = &planeDataConfigInstance3;

SimulationVariables simVars;

int main(
		int i_argc,
		char *const i_argv[]
)
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr))
		return -1;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	int res_x = simVars.disc.space_res_physical[0];
	int res_y = simVars.disc.space_res_physical[1];

	int max_res = 1024;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	double prev_cubic_error_rms = -1;
	double prev_cubic_error_max = -1;

	double prev_linear_error_rms = -1;
	double prev_linear_error_max = -1;

	for (; res_x <= max_res && res_y <= max_res; res_x *= 2, res_y *= 2)
	{
		double resolution_factor = 2;
		int res3[2] = {res_x*(int)resolution_factor, res_y*(int)resolution_factor};

		simVars.disc.space_res_physical[0] = res3[0];
		simVars.disc.space_res_physical[1] = res3[1];
		simVars.disc.space_res_spectral[0] = 0;
		simVars.disc.space_res_spectral[1] = 0;

		//simVars.reset();
		planeDataConfigInstance3.setupAutoSpectralSpace(res3, simVars.misc.reuse_spectral_transformation_plans);
		PlaneData_Spectral prog_h3_pert(planeDataConfig3);

		//PlaneData prog_test3(planeDataConfig3);
		{
			PlaneData_Spectral prog_u3(planeDataConfig3);
			PlaneData_Spectral prog_v3(planeDataConfig3);

			PlaneOperators op3(planeDataConfig3, simVars.sim.plane_domain_size);
			SWEPlaneBenchmarksCombined b;

			b.setupInitialConditions(
					prog_h3_pert,
					//prog_test3,
					prog_u3,
					prog_v3,
					simVars,
					op3
			);
		}

		int res[2] = {res_x, res_y};

		simVars.disc.space_res_physical[0] = res[0];
		simVars.disc.space_res_physical[1] = res[1];
		simVars.disc.space_res_spectral[0] = 0;
		simVars.disc.space_res_spectral[1] = 0;

		//simVars.reset();
		planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);
		PlaneData_Spectral prog_h_pert(planeDataConfig);

		{
			PlaneData_Spectral prog_u(planeDataConfig);
			PlaneData_Spectral prog_v(planeDataConfig);

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size);
			SWEPlaneBenchmarksCombined b;

			b.setupInitialConditions(
					prog_h_pert,
					prog_u,
					prog_v,
					simVars,
					op
			);
		}


		/*
		 * Sampling points with different resolution
		 */
		PlaneData_Physical px(planeDataConfig3);
		PlaneData_Physical py(planeDataConfig3);

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position
		for (int j = 0; j < res3[1]; j++)
		{
			for (int i = 0; i < res3[0]; i++)
			{
				px.physical_set_value(j, i, ((double)i)*(simVars.sim.plane_domain_size[0]/(double)res3[0]));
				py.physical_set_value(j, i, ((double)j)*(simVars.sim.plane_domain_size[1]/(double)res3[1]));
			}
		}


		/*
		 * setup sampler
		 */
		PlaneDataSampler sampler2D;
		sampler2D.setup(simVars.sim.plane_domain_size, planeDataConfig);


		{
			/*
			 * sample with BiLinear interpolation
			 */
			PlaneData_Spectral prog_h3_bilinear(planeDataConfig3);

			sampler2D.bilinear_scalar(
					prog_h_pert,	///< input scalar field
					Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(px),
					Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(py),
					prog_h3_bilinear
			);

			double error_rms = (prog_h3_bilinear-prog_h3_pert).toPhys().physical_reduce_rms();
			double error_max = (prog_h3_bilinear-prog_h3_pert).toPhys().physical_reduce_max_abs();

			double rate_rms = 0;
			double rate_max = 0;

			if (prev_linear_error_rms != -1)
			{
				rate_rms = prev_linear_error_rms/error_rms;
				rate_max = prev_linear_error_max/error_max;

				if (std::isnan(rate_rms))
					SWEETError("NaN detected: error_rms");

				if (std::isnan(rate_max))
					SWEETError("NaN detected: error_max");

				if (
						rate_rms < 3.9 || rate_rms > 4 ||
						rate_max < 3.9 || rate_max > 4
				)
				{
					std::cout << "Bilinear: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\trate_rms: " << rate_rms << "\trate_max: " << rate_max << std::endl;
					std::cerr << "Convergence of 4 expected" << std::endl;
					exit(1);
				}
			}

			std::cout << "Bilinear: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\trate_rms: " << rate_rms << "\trate_max: " << rate_max << std::endl;

			prev_linear_error_rms = error_rms;
			prev_linear_error_max = error_max;
		}

		{
			/*
			 * sample with BiCubic interpolation
			 */
			PlaneData_Spectral prog_h3_bicubic(planeDataConfig3);

			sampler2D.bicubic_scalar(
					prog_h_pert,	///< input scalar field
					Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(px),
					Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(py),
					prog_h3_bicubic	///< output field
			);

	//		double error_norm2 = (prog_h3_bicubic-prog_h3).reduce_norm2()/((double)res3[0] * (double)res3[1]);
			double error_rms = (prog_h3_bicubic-prog_h3_pert).toPhys().physical_reduce_rms();
			double error_max = (prog_h3_bicubic-prog_h3_pert).toPhys().physical_reduce_max_abs();

			double rate_rms = 0;
			double rate_max = 0;
			if (prev_cubic_error_rms != -1)
			{
				rate_rms = prev_cubic_error_rms/error_rms;
				rate_max = prev_cubic_error_max/error_max;

				if (std::isnan(error_rms))
					SWEETError("NaN detected: error_rms");

				if (std::isnan(error_max))
					SWEETError("NaN detected: error_max");

				if (
						//rate_rms < 7.9 || rate_rms > 9 ||
						//rate_max < 7.9 || rate_max > 9
						rate_rms < 15 || rate_rms > 17 ||
						rate_max < 15 || rate_max > 17
				)
				{
					std::cout << "Bicubic: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\trate_rms: " << rate_rms << "\trate_max: " << rate_max << std::endl;
					std::cerr << "ERROR: Convergence of 16 expected for bicubic expected" << std::endl;
					//std::cerr << "Convergence of 8 expected for bicubic expected" << std::endl;
					exit(1);
				}
			}

			std::cout << "Bicubic: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\trate_rms: " << rate_rms << "\trate_max: " << rate_max << std::endl;

			prev_cubic_error_rms = error_rms;
			prev_cubic_error_max = error_max;
		}



	}

	return 0;
}
