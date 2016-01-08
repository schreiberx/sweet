/*
 * lagrangian_test.cpp
 *
 *  Created on: 5 Dec 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */


#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Sampler2D.hpp>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <array>

SimulationVariables simVars;

int main(int i_argc, char *i_argv[])
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr))
		return -1;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.res[0];
	std::size_t res_y = simVars.disc.res[1];

	std::size_t max_res = 2048;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	double prev_cubic_error_rms = -1;
	double prev_cubic_error_max = -1;

	double prev_linear_error_rms = -1;
	double prev_linear_error_max = -1;

	for (; res_x <= max_res && res_y <= max_res; res_x *= 2, res_y *= 2)
	{
		simVars.disc.res[0] = res_x;
		simVars.disc.res[1] = res_y;

		std::size_t res[2] = {res_x, res_y};

		DataArray<2> prog_h(res);

		// setup initial conditions
		for (std::size_t j = 0; j < res[1]; j++)
		{
			for (std::size_t i = 0; i < res[0]; i++)
			{
				double x = (double)i*(simVars.sim.domain_size[0]/(double)res[0]);
				double y = (double)j*(simVars.sim.domain_size[1]/(double)res[1]);

				prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
			}
		}


		double resolution_factor = 2;

		std::size_t res3[2] = {res_x*(int)resolution_factor, res_y*(int)resolution_factor};

		DataArray<2> prog_h3(res3);

		// setup initial conditions
		{
			simVars.disc.res[0] = res3[0];
			simVars.disc.res[1] = res3[1];

			for (std::size_t j = 0; j < res3[1]; j++)
			{
				for (std::size_t i = 0; i < res3[0]; i++)
				{
					double x = ((double)i/resolution_factor)*(simVars.sim.domain_size[0]/(double)res3[0]);
					double y = ((double)j/resolution_factor)*(simVars.sim.domain_size[1]/(double)res3[1]);

					prog_h3.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
				}
			}

			simVars.disc.res[0] = res[0];
			simVars.disc.res[1] = res[1];
		}

		/*
		 * Sampling points with different resolution
		 */
		DataArray<2> px(res3);
		DataArray<2> py(res3);

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position
		for (std::size_t j = 0; j < res3[1]; j++)
		{
			for (std::size_t i = 0; i < res3[0]; i++)
			{
				px.set(j, i, ((double)i/resolution_factor)*(simVars.sim.domain_size[0]/(double)res3[0]));
				py.set(j, i, ((double)j/resolution_factor)*(simVars.sim.domain_size[1]/(double)res3[1]));
			}
		}


		/*
		 * setup sampler
		 */
		Sampler2D sampler2D;
		sampler2D.setup(simVars.sim.domain_size, res);


		{
			/*
			 * sample with BiLinear interpolation
			 */
			DataArray<2> prog_h3_bilinear(res3);

			sampler2D.bilinear_scalar(
					prog_h,	///< input scalar field
					px, py,		///< sampling positions
					prog_h3_bilinear	///< output field
			);

	//		double error_norm2 = (prog_h3_bicubic-prog_h3).reduce_norm2()/((double)res3[0] * (double)res3[1]);
			double error_rms = (prog_h3_bilinear-prog_h3).reduce_rms();
			double error_max = (prog_h3_bilinear-prog_h3).reduce_maxAbs();

			double rate_rms = 0;
			double rate_max = 0;

			if (prev_linear_error_rms != -1)
			{
				rate_rms = prev_linear_error_rms/error_rms;
				rate_max = prev_linear_error_max/error_max;

				if (
						rate_rms < 3.9 || rate_rms > 4 ||
						rate_max < 3.9 || rate_max > 4
				)
				{
					std::cout << "Bilinear: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\t" << rate_rms << "\t" << rate_max << std::endl;
					std::cerr << "Convergence of 4 expected" << std::endl;
//					exit(1);
				}
			}

			std::cout << "Bilinear: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\t" << rate_rms << "\t" << rate_max << std::endl;

			prev_linear_error_rms = error_rms;
			prev_linear_error_max = error_max;
		}

		{
			/*
			 * sample with BiCubic interpolation
			 */
			DataArray<2> prog_h3_bicubic(res3);

			sampler2D.bicubic_scalar(
					prog_h,	///< input scalar field
					px, py,		///< sampling positions
					prog_h3_bicubic	///< output field
			);

	//		double error_norm2 = (prog_h3_bicubic-prog_h3).reduce_norm2()/((double)res3[0] * (double)res3[1]);
			double error_rms = (prog_h3_bicubic-prog_h3).reduce_rms();
			double error_max = (prog_h3_bicubic-prog_h3).reduce_maxAbs();

			double rate_rms = 0;
			double rate_max = 0;
			if (prev_cubic_error_rms != -1)
			{
				rate_rms = prev_cubic_error_rms/error_rms;
				rate_max = prev_cubic_error_max/error_max;

				if (
						rate_rms < 7.9 || rate_rms > 9 ||
						rate_max < 7.9 || rate_max > 9
				)
				{
					std::cout << "Bicubic: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\t" << rate_rms << "\t" << rate_max << std::endl;
					std::cerr << "Convergence of 8 expected for bicubic " << std::endl;
					exit(1);
				}
			}

			std::cout << "Bicubic: " << res_x << "x" << res_y << "\t" << error_rms << "\t" << error_max << "\t" << rate_rms << "\t" << rate_max << std::endl;

			prev_cubic_error_rms = error_rms;
			prev_cubic_error_max = error_max;
		}



	}

	return 0;
}
