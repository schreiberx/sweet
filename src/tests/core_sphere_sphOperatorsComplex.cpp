/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/sphere/SphereComplex.hpp>
#include <sweet/core/SWEETError.hpp>

#include "core_sphere_sphOperators/SphereTestSolutions_Gaussian.hpp"


void test_header(const std::string &i_str)
{
	std::cout << "**********************************************" << std::endl;
	std::cout << i_str << std::endl;
	//std::cout << "**********************************************" << std::endl;
}

void run_tests(
		sweet::SphereData_Config *i_sphereDataConfig,
		sweet::ShackSphereDataOps *i_shackSphereDataOps
)
{
	//double eps = 1e-10;
	double eps = 1e-8;
	eps *= std::sqrt(i_sphereDataConfig->spectral_modes_n_max)*std::sqrt(i_sphereDataConfig->spectral_modes_m_max);
	std::cout << "Using max allowed error of " << eps << std::endl;

	// Use earth radius of 1
	i_shackSphereDataOps->sphere_radius = 1.0;
	sweet::SphereOperatorsComplex ops(i_sphereDataConfig, i_shackSphereDataOps);


	if (true)
	{
		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0;
		double a = 6.37122e6;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		double alpha[] = {0, M_PI/3, M_PI/2};

		{
			test_header("Testing divergence freeness (non-Robert)");

			for (int i = 0; i < sizeof(alpha)/sizeof(double); i++)
			{
				double advection_rotation_angle = alpha[i];

				std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;

				sweet::SphereData_PhysicalComplex u(i_sphereDataConfig);
				u.physical_update_lambda(
					[&](double i_lon, double i_lat, std::complex<double> &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
								u0*(
									std::cos(i_theta)*std::cos(advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(advection_rotation_angle)
							);
					}
				);

				sweet::SphereData_PhysicalComplex v(i_sphereDataConfig);
				v.physical_update_lambda(
					[&](double i_lon, double i_lat, std::complex<double> &io_data)
					{
						//double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
									std::sin(i_lambda)*std::sin(advection_rotation_angle)
							);
					}
				);

				sweet::SphereData_SpectralComplex vort(i_sphereDataConfig);
				sweet::SphereData_SpectralComplex div(i_sphereDataConfig);
				ops.uv_to_vrtdiv(u, v, vort, div);

				double div_max_error = div.toPhys().physical_reduce_max_abs();
				std::cout << " + div_max_error: " << div_max_error << std::endl;

				if (div_max_error > eps)
					SWEETError(" + ERROR! max error exceeds threshold");
			}
		}
	}


	{
		SphereTestSolutions_Gaussian testSolutions;

		if (true)
		{
			test_header("Testing Multiplication (a*b) with b=123.0");

			sweet::SphereData_SpectralComplex data(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex data_phys(i_sphereDataConfig);

			data_phys.physical_update_lambda(
				[&](double x, double y, std::complex<double> &io_data)
				{
					io_data = y;
				}
			);

			data.loadSphereDataPhysical(data_phys);
			data = data*123.0;

			sweet::SphereData_SpectralComplex data2(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex data2_phys(i_sphereDataConfig);
			data2_phys.physical_update_lambda(
				[&](double x, double y, std::complex<double> &io_data)
				{
					io_data = y*123.0;
				}
			);
			data2.loadSphereDataPhysical(data2_phys);

			double div_max_error = (data-data2).toPhys().physical_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETError(" + ERROR! max error exceeds threshold");
		}



		{
			SphereTestSolutions_Gaussian testSolutions;

			if (true)
			{
				test_header("Testing Multiplication (a *= b) with b=123.0");

				sweet::SphereData_SpectralComplex data(i_sphereDataConfig);
				sweet::SphereData_PhysicalComplex data_phys(i_sphereDataConfig);
				data_phys.physical_update_lambda(
						[&](double x, double y, std::complex<double> &io_data)
						{
							io_data = y;
						}
				);
				data.loadSphereDataPhysical(data_phys);

				data *= 123.0;

				sweet::SphereData_SpectralComplex data2(i_sphereDataConfig);
				sweet::SphereData_PhysicalComplex data2_phys(i_sphereDataConfig);
				data2_phys.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y*123.0;
					}
				);
				data2.loadSphereDataPhysical(data2_phys);

				double div_max_error = (data-data2).toPhys().physical_reduce_max_abs();
				std::cout << " + div_max_error: " << div_max_error << std::endl;

				if (div_max_error > eps)
					SWEETError(" + ERROR! max error exceeds threshold");
			}
		}

		if (true)
		{
			test_header("Testing add (a+b) operation with 123.0");

			sweet::SphereData_SpectralComplex data(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex data_phys(i_sphereDataConfig);
			data_phys.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y;
					}
			);

			data.loadSphereDataPhysical(data_phys);
			data = data + 123.0;

			sweet::SphereData_SpectralComplex data2(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex data2_phys(i_sphereDataConfig);
			data2_phys.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y+123.0;
					}
			);
			data2.loadSphereDataPhysical(data2_phys);

			double div_max_error = (data-data2).toPhys().physical_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETError(" + ERROR! max error exceeds threshold");
		}


		if (true)
		{
			test_header("Testing add (a+=b) operation with 123.0");

			sweet::SphereData_SpectralComplex data(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex data_phys(i_sphereDataConfig);
			data_phys.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y;
					}
			);
			data.loadSphereDataPhysical(data_phys);

			data += 123.0;

			sweet::SphereData_SpectralComplex data2(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex data2_phys(i_sphereDataConfig);
			data2_phys.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y+123.0;
					}
			);
			data2.loadSphereDataPhysical(data2_phys);

			double div_max_error = (data-data2).toPhys().physical_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETError(" + ERROR! max error exceeds threshold");
		}

		if (true)
		{
			test_header("Testing Gaussian latitude coordinates");

			sweet::SphereData_SpectralComplex h(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex h_phys(i_sphereDataConfig);
			h_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function__grid_gaussian(a,b,val);
						c = val;
					}
			);
			h.loadSphereDataPhysical(h_phys);

			sweet::SphereData_SpectralComplex hphi(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex hphi_phys(i_sphereDataConfig);
			hphi_phys.physical_update_lambda(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function_phi__grid_phi(a,b,val);
						c = val;
					}
			);
			hphi.loadSphereDataPhysical(hphi_phys);

			double div_max_error = (h-hphi).toPhys().physical_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETError(" + ERROR! max error exceeds threshold");
		}


		if (true)
		{
			test_header("Testing multiplication with Gaussian latitude");

			// mu*F(\lambda,\mu)
			sweet::SphereData_SpectralComplex h(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex h_phys(i_sphereDataConfig);
			h_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function__grid_gaussian(a,b,val);
						c = val;
					}
			);
			h.loadSphereDataPhysical(h_phys);
			h = ops.mu(h);

			sweet::SphereData_SpectralComplex result(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex result_phys(i_sphereDataConfig);
			result_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.correct_result_mu__grid_gaussian(a,b,val);
						c = val;
					}
			);
			result.loadSphereDataPhysical(result_phys);

			double div_max_error = (h-result).toPhys().physical_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETError(" + ERROR! max error exceeds threshold");
		}

		if (true)
		{
			test_header("Testing multiplication with pow2 of Gaussian latitude");

			// mu*mu*F(\lambda,\mu)
			sweet::SphereData_SpectralComplex h(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex h_phys(i_sphereDataConfig);
			h_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double val;
						testSolutions.test_function__grid_gaussian(a,b,val);
						c = val;
					}
			);
			h.loadSphereDataPhysical(h_phys);
			h = ops.mu2(h);

			sweet::SphereData_SpectralComplex result(i_sphereDataConfig);
			sweet::SphereData_PhysicalComplex result_phys(i_sphereDataConfig);
			result_phys.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &i_data)
					{
						double val;
						testSolutions.test_function__grid_gaussian(lat, mu, val);
						i_data = val;
						i_data *= mu*mu;
					}
			);
			result.loadSphereDataPhysical(result_phys);


			double div_max_error = (h-result).toPhys().physical_reduce_max_abs();
			std::cout << " + div_max_error: " << div_max_error << std::endl;

			if (div_max_error > eps)
				SWEETError(" + ERROR! max error exceeds threshold");
		}
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{

	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::ShackSphereDataOps *shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	sweet::SphereData_Config sphereData_Config;
	sphereData_Config.setupAuto(shackSphereDataOps);
	run_tests(&sphereData_Config, shackSphereDataOps);

	std::cout << "All test successful" << std::endl;
}

