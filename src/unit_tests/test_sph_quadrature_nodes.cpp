/*
 * test_sph_quadrature_nodes.hpp
 *
 *  Created on: 3 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_TESTSPH_QUADRATURE_NODES_HPP_
#define SRC_TESTSPH_QUADRATURE_NODES_HPP_


#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <libmath/BandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/SphereDataErrorCheck.hpp>
#include <sweet/sphere/SphereDataSpectralComplex.hpp>


SimulationVariables simVars;

void run_tests(
		SphereDataConfig *sphereDataConfig
)
{
	double epsilon = 1e-12;
	epsilon *= (sphereDataConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;

	std::cout << std::setprecision(10);

	{
		SphereDataSpectral physical(sphereDataConfig);
		physical.physical_update_lambda_cogaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(2.0)*0.5;
			}
		);

		SphereDataSpectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 0 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereDataErrorCheck::check(physical, spectral, "n=0, m=0", epsilon, false, true);
	}

	{
		SphereDataSpectral physical(sphereDataConfig);
		physical.physical_update_lambda_cogaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(6.0)*mu*0.5;
			}
		);

		SphereDataSpectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 1 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereDataErrorCheck::check(physical, spectral, "n=1, m=0", epsilon, false, true);
	}

	{
		SphereDataSpectral physical(sphereDataConfig);
		physical.physical_update_lambda_cogaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(10.0)/4.0 * (3.0*mu*mu - 1.0);
			}
		);

		SphereDataSpectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 2 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereDataErrorCheck::check(physical, spectral, "n=2, m=0", epsilon, false, true);
	}


	{
		SphereDataSpectral physical(sphereDataConfig);
		physical.physical_update_lambda_cogaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(14.0)/4.0*mu * (5.0*mu*mu - 3.0);
			}
		);

		SphereDataSpectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 3 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereDataErrorCheck::check(physical, spectral, "n=3, m=0", epsilon, false, true);
	}
}




int main(
		int i_argc,
		char *const i_argv[]
)
{
	MemBlockAlloc numaBlockAlloc;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.disc.space_res_spectral[0] == 0)
		FatalError("Set number of spectral modes to use SPH!");

	SphereDataConfig sphereDataConfig;
	sphereDataConfig.setupAutoPhysicalSpace(
					simVars.disc.space_res_spectral[0],
					simVars.disc.space_res_spectral[1],
					&simVars.disc.space_res_physical[0],
					&simVars.disc.space_res_physical[1]
			);

	run_tests(&sphereDataConfig);
}



#endif
