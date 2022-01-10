/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "test_sphere_sph_operators/SphereTestSolutions_Gaussian.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SWEETError.hpp>



SimulationVariables simVars;

SphereData_Config sphereDataConfigInstance;
SphereData_Config *sphereDataConfig = &sphereDataConfigInstance;


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
	simVars.sim.sphere_radius = 1.0;
	SphereOperators_SphereData op(sphereDataConfig, &(simVars.sim));


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

				SphereData_Physical u(sphereDataConfig);
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
					}
				);

				SphereData_Physical v(sphereDataConfig);
				v.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						//double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
									std::sin(i_lambda)*std::sin(advection_rotation_angle)
							);
					}
				);

				SphereData_Spectral vort(sphereDataConfig);
				SphereData_Spectral div(sphereDataConfig);
				op.uv_to_vrtdiv(u, v, vort, div);

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

			SphereData_Spectral data(sphereDataConfig);
			SphereData_Physical data_phys(sphereDataConfig);

			data_phys.physical_update_lambda(
				[&](double x, double y, double &io_data)
				{
					io_data = y;
				}
			);

			data.loadSphereDataPhysical(data_phys);
			data = data*123.0;

			SphereData_Spectral data2(sphereDataConfig);
			SphereData_Physical data2_phys(sphereDataConfig);
			data2_phys.physical_update_lambda(
				[&](double x, double y, double &io_data)
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

				SphereData_Spectral data(sphereDataConfig);
				SphereData_Physical data_phys(sphereDataConfig);
				data_phys.physical_update_lambda(
						[&](double x, double y, double &io_data)
						{
							io_data = y;
						}
				);
				data.loadSphereDataPhysical(data_phys);

				data *= 123.0;

				SphereData_Spectral data2(sphereDataConfig);
				SphereData_Physical data2_phys(sphereDataConfig);
				data2_phys.physical_update_lambda(
					[&](double x, double y, double &io_data)
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

			SphereData_Spectral data(sphereDataConfig);
			SphereData_Physical data_phys(sphereDataConfig);
			data_phys.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y;
					}
			);

			data.loadSphereDataPhysical(data_phys);
			data = data + 123.0;

			SphereData_Spectral data2(sphereDataConfig);
			SphereData_Physical data2_phys(sphereDataConfig);
			data2_phys.physical_update_lambda(
					[&](double x, double y, double &io_data)
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

			SphereData_Spectral data(sphereDataConfig);
			SphereData_Physical data_phys(sphereDataConfig);
			data_phys.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y;
					}
			);
			data.loadSphereDataPhysical(data_phys);

			data += 123.0;

			SphereData_Spectral data2(sphereDataConfig);
			SphereData_Physical data2_phys(sphereDataConfig);
			data2_phys.physical_update_lambda(
					[&](double x, double y, double &io_data)
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

			SphereData_Spectral h(sphereDataConfig);
			SphereData_Physical h_phys(sphereDataConfig);
			h_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c)
					{
						testSolutions.test_function__grid_gaussian(a,b,c);
					}
			);
			h.loadSphereDataPhysical(h_phys);

			SphereData_Spectral hphi(sphereDataConfig);
			SphereData_Physical hphi_phys(sphereDataConfig);
			hphi_phys.physical_update_lambda(
					[&](double a, double b, double &c){testSolutions.test_function_phi__grid_phi(a,b,c);}
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
			SphereData_Spectral h(sphereDataConfig);
			SphereData_Physical h_phys(sphereDataConfig);
			h_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h.loadSphereDataPhysical(h_phys);
			h = op.mu(h);

			SphereData_Spectral result(sphereDataConfig);
			SphereData_Physical result_phys(sphereDataConfig);
			result_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_mu__grid_gaussian(a,b,c);}
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
			SphereData_Spectral h(sphereDataConfig);
			SphereData_Physical h_phys(sphereDataConfig);
			h_phys.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h.loadSphereDataPhysical(h_phys);
			h = op.mu2(h);

			SphereData_Spectral result(sphereDataConfig);
			SphereData_Physical result_phys(sphereDataConfig);
			result_phys.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &i_data){
						testSolutions.test_function__grid_gaussian(lat, mu, i_data);
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

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.disc.space_res_spectral[0] == 0)
		SWEETError("Set number of spectral modes to use SPH!");

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

	run_tests();

	std::cout << "All test successful" << std::endl;
}

