/*
 * test_sph_operators_div_freeness.hpp
 *
 *  Created on: 07 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef TEST_SPH_OPERATORS_DIV_FREENESS_HPP_
#define TEST_SPH_OPERATORS_DIV_FREENESS_HPP_

#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/FatalError.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereDataErrorCheck.hpp>
#include <sweet/sphere/SphereDataSpectralComplex.hpp>



SimulationVariables simVars;

SphereDataConfig sphereDataConfigInstance;
SphereDataConfig sphereDataConfigExtInstance;

SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigExtInstance;



void run_tests()
{
	simVars.outputConfig();

	double error_threshold = 1e-10;
	error_threshold *= (sphereDataConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << error_threshold << std::endl;

	SphereOperators op(sphereDataConfig, 1);
	SphereOperatorsComplex opComplex(sphereDataConfig, 1);


	if (true)
	{
		std::cout << "****************************************" << std::endl;
		std::cout << "* DIV FREE TESTS ROBERT *" << std::endl;
		std::cout << "****************************************" << std::endl;


		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0;
		double a = 6.37122e6;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		double alpha[] = {0, M_PI/3, M_PI/2};

		SphereDataSpectral zero(sphereDataConfig);
		zero.physical_set_zero();

		SphereDataSpectralComplex zeroc(sphereDataConfig);
		zeroc.physical_set_zero();

		for (int i = 0; i < 3; i++)
		{
			double advection_rotation_angle = alpha[i];

			std::cout << "****************************************" << std::endl;
			std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;

			SphereDataSpectral phi(sphereDataConfig);
			phi.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					double r = a * std::acos(
							std::sin(theta_c)*std::sin(i_theta) +
							std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
					);

					if (r < R)
						io_data = simVars.sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
					else
						io_data = 0;
				}
			);

			SphereDataSpectral u(sphereDataConfig);
			u.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
					double i_theta = i_lat;
					double i_lambda = i_lon;
					io_data =
							u0*(
								std::cos(i_theta)*std::cos(advection_rotation_angle) +
								std::sin(i_theta)*std::cos(i_lambda)*std::sin(advection_rotation_angle)
						);

					io_data *= cos(i_lat);
				}
			);

			SphereDataSpectral v(sphereDataConfig);
			v.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
//					double i_phi = i_lat;
					double i_lambda = i_lon;
					io_data =
						-u0*(
								std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);

					io_data *= cos(i_lat);
				}
			);



			SphereDataSpectralComplex phic(sphereDataConfig);
			phic.physical_update_lambda(
				[&](double i_lambda, double i_theta, std::complex<double> &io_data)
				{
					double r = a * std::acos(
							std::sin(theta_c)*std::sin(i_theta) +
							std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
					);

					if (r < R)
						io_data = simVars.sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
					else
						io_data = 0;
				}
			);

			SphereDataSpectralComplex uc(sphereDataConfig);
			uc.physical_update_lambda(
				[&](double i_lon, double i_lat, std::complex<double> &io_data)
				{
					double i_theta = i_lat;
					double i_lambda = i_lon;
					io_data =
							u0*(
								std::cos(i_theta)*std::cos(advection_rotation_angle) +
								std::sin(i_theta)*std::cos(i_lambda)*std::sin(advection_rotation_angle)
						);

					io_data *= cos(i_lat);
				}
			);

			SphereDataSpectralComplex vc(sphereDataConfig);
			vc.physical_update_lambda(
				[&](double i_lon, double i_lat, std::complex<double> &io_data)
				{
//					double i_phi = i_lat;
					double i_lambda = i_lon;
					io_data =
						-u0*(
								std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);

					io_data *= cos(i_lat);
				}
			);

		}
	}
}

int main(
		int i_argc,
		char *const i_argv[]
)
{
	/*
	 * Initialize NUMA block allocator
	 */
	MemBlockAlloc numaBlockAlloc;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.disc.space_res_spectral[0] == 0)
		FatalError("Set number of spectral modes to use SPH!");

	if (simVars.disc.space_res_physical[0] <= 0)
	{
		sphereDataConfigInstance.setupAutoPhysicalSpace(
						simVars.disc.space_res_spectral[0],
						simVars.disc.space_res_spectral[1],
						&simVars.disc.space_res_physical[0],
						&simVars.disc.space_res_physical[1],
						simVars.misc.reuse_spectral_transformation_plans

				);
	}
	else
	{
		sphereDataConfigInstance.setup(
						simVars.disc.space_res_spectral[0],
						simVars.disc.space_res_spectral[1],
						simVars.disc.space_res_physical[0],
						simVars.disc.space_res_physical[1],
						simVars.misc.reuse_spectral_transformation_plans

				);
	}

	sphereDataConfigExtInstance.setupAdditionalModes(
			&sphereDataConfigInstance,
			simVars.rexi.use_sphere_extended_modes,
			simVars.rexi.use_sphere_extended_modes,
			simVars.misc.reuse_spectral_transformation_plans

		);

	run_tests();
}


#endif /* SRC_TESTOPERATORS_HPP_ */
