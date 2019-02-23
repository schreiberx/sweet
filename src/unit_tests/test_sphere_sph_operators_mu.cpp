/*
 * test_sphere_sph_operators_mu.cpp
 *
 *  Created on: 05 Nov 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/FatalError.hpp>
#include <sweet/MemBlockAlloc.hpp>



SimulationVariables simVars;

SphereDataConfig sphereDataConfigInstance;
SphereDataConfig sphereDataConfigExtInstance;

SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigExtInstance;


void test_header(const std::string &i_str)
{
	std::cout << "**********************************************" << std::endl;
	std::cout << i_str << std::endl;
	//std::cout << "**********************************************" << std::endl;
}

void run_tests()
{
	simVars.outputConfig();

	//double eps = 1e-10;
	double eps = 1e-8;
	eps *= std::sqrt(sphereDataConfig->spectral_modes_n_max)*std::sqrt(sphereDataConfig->spectral_modes_m_max);
	std::cout << "Using max allowed error of " << eps << std::endl;

	// Use earth radius of 1
	SphereOperators op(sphereDataConfig, 1);


	{
		test_header("Testing Coriolis operators");

		SphereDataSpectral one(sphereDataConfig);
		one.physical_set_all_value(1.0);

		/*
		 * Coriolis frequency
		 */
		SphereDataSpectral f(sphereDataConfig);
		f.physical_update_lambda_gaussian_grid(
			[&](double i_lon, double i_lat_gauss, double &io_data)
			{
				io_data = i_lat_gauss*2.0*simVars.sim.sphere_rotating_coriolis_omega;
			}
		);
		double max_f = f.physical_reduce_max_abs();	// 2*\Omega

		f.spectral_print();


		double eps = 1e-12;

		eps *= std::sqrt(sphereDataConfig->spectral_modes_n_max)*std::sqrt(sphereDataConfig->spectral_modes_m_max);
		std::cout << "Using max allowed error of " << eps << std::endl;

		for (int m_lon = 0; m_lon <= sphereDataConfig->spectral_modes_m_max; m_lon++)
		{
			for (int m_lat = m_lon; m_lat <= sphereDataConfig->spectral_modes_n_max; m_lat++)
			{
				std::cout << std::endl;
				std::cout << "Mode lon=" << m_lon << ", lat=" << m_lat << std::endl;

				SphereDataSpectral test_fun(sphereDataConfig);
				test_fun.spectral_set_zero();
				test_fun.spectral_set(m_lat, m_lon, 1);

				SphereDataSpectral physical_test_mul_f = (test_fun*f);
				SphereDataSpectral operator_test_mul_f = op.mu(test_fun)*2.0*simVars.sim.sphere_rotating_coriolis_omega;

				double forward_diff_error_linf = (physical_test_mul_f - operator_test_mul_f).physical_reduce_max_abs();
				forward_diff_error_linf /= max_f;

				std::cout << "forward error_linf: " << forward_diff_error_linf << std::endl;

				if (forward_diff_error_linf > eps)
				{
					std::cout << "test_fun:" << std::endl;
					test_fun.spectral_print();

					std::cout << "physical_test_mul_f:" << std::endl;
					physical_test_mul_f.spectral_print();

					std::cout << "operator_test_mul_f:" << std::endl;
					operator_test_mul_f.spectral_print();

					FatalError(" + ERROR! max error exceeds threshold");
				}

				/*
				 * Is this frequency already on the limit of the represantable spectrum?
				 */
				bool freq_limit = m_lat == sphereDataConfig->spectral_modes_n_max;

				SphereDataSpectral invert_physical_test_mul_f = physical_test_mul_f / f;

				double backward_error_linf = (test_fun-invert_physical_test_mul_f).physical_reduce_max_abs();
				std::cout << "backward error_linf: " << backward_error_linf << std::endl;

				if (backward_error_linf > eps)
				{
					if (freq_limit)
					{
						std::cout << "ERROR, but ignoring this, since we're on the limit of the representable spectrum" << std::endl;
					}
					else
					{
						std::cout << "physical_test_mul_f:" << std::endl;
						physical_test_mul_f.spectral_print();

						std::cout << "test_fun:" << std::endl;
						test_fun.spectral_print();

						std::cout << "invert_physical_test_mul_f:" << std::endl;
						invert_physical_test_mul_f.spectral_print();

						FatalError(" + ERROR! max error exceeds threshold");
					}
				}
			}
		}

/*
		(phi*one).spectral_print();

		phi = phi*phi;
		phi.spectral_print();

		phi = phi*phi;
		phi.spectral_print();

		phi = phi*phi;
		phi.spectral_print();
*/

	}


};



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

	std::cout << "All test successful" << std::endl;
}

