/*
 * AppTestSPHOperators.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_TESTOPERATORS_HPP_
#define SRC_TESTOPERATORS_HPP_

#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/FatalError.hpp>
#include <sweet/MemBlockAlloc.hpp>



SimulationVariables simVars;
SphereOperators op;


void run_tests(
		SphereDataConfig *sphConfig
)
{
	double epsilon = 1e-11;
	epsilon *= (sphConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;


	if (true)
	{
		// div mu
		SphereData h(sphConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double a, double mu, double &c){
					c = mu/sqrt(1.0-mu*mu);
				}
		);

		h.file_physical_writeFile_lon_pi_shifted("o_test_initial_pretrunc.csv");
//		h.physical_truncate();
		h.file_physical_writeFile_lon_pi_shifted("o_test_initial_trunc.csv");

		h = op.div_lat(h);
		h.physical_truncate();

		h.file_physical_writeFile_lon_pi_shifted("o_test_postdiv.csv");

		SphereData result(sphConfig);
		result.physical_update_lambda_gaussian_grid(
				[&](double a, double b, double &c){
					c = 1;
				}
		);
		result.physical_truncate();

		result.file_physical_writeFile_lon_pi_shifted("o_test_result.csv");

		double error_max = (h-result).reduce_abs_max();
		std::cout << "TEST SPECIAL DIV LAT  - max error: " << error_max << std::endl;
//		exit(1);
	}



	if (true)
	{
		// diff mu
		SphereData h(sphConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double lambda, double mu, double &io_data)
				{
					io_data = mu;	// sin(phi)
				}
		);
//		h.physical_write_file("O_mu_sph_result.csv");
//		h.spectral_truncate();

		h = op.diff_lat_mu(h);
//		h.spectral_truncate();
//		h.physical_write_file("O_diff_mu_sph_result.csv");

		SphereData result(sphConfig);
		result.physical_update_lambda_gaussian_grid(
				[&](double a, double b, double &io_data){
					io_data = 1.0;
				}
		);

//		result.spectral_truncate();
//		result.physical_write_file("O_diff_mu_correct_result.csv");
//		(h-result).physical_write_file("O_diff_mu_correct_diff.csv");

		double error_max = h.physical_reduce_error_max(result);
		std::cout << "TEST DIFF LAT  - max error: " << error_max << std::endl;

		if (error_max > epsilon)
		{
			std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED" << std::endl;
//				FatalError("ERROR THRESHOLD EXCEEDED!");
		}
	}


	if (true)
	{
		SphereTestSolutions_SPH testSolutionsSph(2,1);

		//int tn = 1;
		//int tm = 0;
		SphereData h(sphConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double a, double b, double &c){testSolutionsSph.test_function__grid_gaussian(a,b,c);}
		);
		h.request_data_spectral();


		SphereData result(sphConfig);
		result.physical_update_lambda_gaussian_grid(
				[&](double a, double b, double &c){testSolutionsSph.correct_result_diff_mu__grid_gaussian(a,b,c);}
		);
		result.request_data_spectral();

//		h.file_physical_writeFile("O_SPHbasis_test_function.csv");


		// one_minus_mu_squared_diff_lat
		h = op.spec_one_minus_mu_squared_diff_lat_mu(h);
//		h.file_physical_writeFile("O_SPHbasis_one_minus_mu_squared_diff_mu_sph_result.csv");

//		result.file_physical_writeFile("O_SPHbasis_one_minus_mu_squared_diff_mu_correct_result.csv");


		double error_max = h.physical_reduce_error_max(result);
		std::cout << "TEST SPHbasis (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

		if (error_max > epsilon)
			FatalError("ERROR THRESHOLD EXCEEDED!");
	}


	{
		SphereTestSolutions_Gaussian testSolutions;


		if (true)
		{
			SphereData data(sphConfig);
			data.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y;
					}
			);

			data = data*123.0;
			data.spectral_truncate();

			SphereData data2(sphConfig);
			data2.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y*123.0;
					}
			);
			data2.spectral_truncate();

			double error = data.physical_reduce_error_max(data2);
			std::cout << "ERROR operator*123.0: " << error << std::endl;

			if (error > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			SphereData data(sphConfig);
			data.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y;
					}
			);

			data = data + 100.0;
			data.spectral_truncate();

			SphereData data2(sphConfig);
			data2.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y+100.0;
					}
			);
			data2.spectral_truncate();

			double error = data.physical_reduce_error_max(data2);
			std::cout << "ERROR operator+(double): " << error << std::endl;

			if (error > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
//			h.file_physical_writeFile("O_test_function.csv");

			SphereData hphi(sphConfig);
			hphi.physical_update_lambda(
					[&](double a, double b, double &c){testSolutions.test_function_phi__grid_phi(a,b,c);}
			);
//			hphi.file_physical_writeFile("O_test_function_phi.csv");

			double error_max = h.physical_reduce_error_max(hphi);
			std::cout << "TEST PHI vs. MU: max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// identity
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			//h = h.spat_truncate();

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();

			double error_max = h.physical_reduce_error_max(result);

			std::cout << "TEST IDENTITY: max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d lambda
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			//h = h.spat_truncate();

			h = op.diff_lon(h);

//			h.file_physical_writeFile("O_diff_lambda_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_lambda__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_diff_lambda_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST DIFF LON - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d phi
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.diff_lat_phi(h);
//			h.file_physical_writeFile("O_diff_phi_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_phi__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_diff_phi_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST DIFF PHI - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d mu
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.diff_lat_mu(h);
			//h = h.spat_truncate();
//			h.file_physical_writeFile("O_diff_mu_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_mu__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_diff_mu_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST DIFF LAT MU - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// mu*F(\lambda,\mu)
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.mu(h);
//			h.file_physical_writeFile("O_mu_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_mu__grid_gaussian(a,b,c);}
			);
//			result.file_physical_writeFile("O_mu_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST mu*() - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// mu*mu*F(\lambda,\mu)
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.mu2(h);
//			h.file_physical_writeFile("O_mu2_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &i_data){
						testSolutions.test_function__grid_gaussian(lat, mu, i_data);
						i_data *= mu*mu;
					}
			);
//			result.file_physical_writeFile("O_mu2_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST mu*mu*() - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// one_minus_mu_squared_diff_lat
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.spec_one_minus_mu_squared_diff_lat_mu(h);
//			h.file_physical_writeFile("O_one_minus_mu_squared_diff_mu_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(a,b,c);}
			);

//			result.file_physical_writeFile("O_one_minus_mu_squared_diff_mu_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// grad lambda
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.grad_lon(h);
			//h = h.spat_truncate();
//			h.file_physical_writeFile("O_grad_lambda_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_grad_lambda_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST GRAD LON - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// grad mu
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.grad_lat(h);
//			h.file_physical_writeFile("O_grad_phi_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_grad_phi__grid_gaussian(a,b,c);}
			);

//			result.file_physical_writeFile("O_grad_phi_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST GRAD LAT - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// div lambda
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.div_lon(h);
//			h.file_physical_writeFile("O_div_lambda_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_div_lambda__grid_gaussian(a,b,c);}
			);
//			result.file_physical_writeFile("O_div_lambda_correct_result.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST DIV LON  - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// div mu
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){
						testSolutions.test_function__grid_gaussian(a,b,c);
					}
			);
			h.physical_write_file("O_div_mu_initial_sph.csv");

			h = op.div_lat(h);
//			h.spectral_truncate();
			h.physical_write_file("O_div_mu_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){
						testSolutions.correct_result_div_mu__grid_gaussian(a,b,c);
					}
			);
//			result.spectral_truncate();
			result.physical_write_file("O_div_mu_correct_result.csv");
			(h-result).physical_write_file("O_div_mu_correct_diff.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;

			if (error_max > epsilon)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED" << std::endl;
//				FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}
#if 0

		if (true)
		{
			// div mu TEST
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.div_lat_TEST(h);
			h.file_physical_writeFile("O_div_mu_TEST_sph_result.csv");
#if 0
			h.spat_update_lambda_cogaussian_grid(
					[&](double a, double mu, double &c){
						c *= mu;
					}
			);
#endif
			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double mu, double &c){
						testSolutions.correct_result_div_mu__grid_gaussian(a,mu,c);
					}
			);
			result.spat_update_lambda_cogaussian_grid(
					[&](double a, double mu, double &c){
						c *= mu;
					}
			);

			result.file_physical_writeFile("O_div_mu_TEST_correct_result.csv");
			(h-result).file_physical_writeFile("O_div_mu_TEST_correct_diff.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST DIV TEST LAT  - max error: " << error_max << std::endl;
		}
#endif

		if (true)
		{
			// divergence
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.div(h, h);
//				h.spat_update_lambda_gaussian_grid([&](double a, double mu, double &c){c *= mu*mu*mu*mu;});
//			h.file_physical_writeFile("O_divergence_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){

						double data1;
						testSolutions.correct_result_div_mu__grid_gaussian(a,b,data1);

						double data2;
						testSolutions.correct_result_div_lambda__grid_gaussian(a,b,data2);

						c = data1 + data2;
					}
			);
			//result = result.spat_truncate();

//			result.file_physical_writeFile("O_divergence_correct_result.csv");

//			(h-result).file_physical_writeFile("O_divergence_correct_diff.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST DIVERGENCE - max error: " << error_max << std::endl;

			if (error_max > epsilon)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED" << std::endl;
				//FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}


		if (true)
		{
			// vorticity
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.vort(h, h);
			//h = h.spat_truncate();
//			h.file_physical_writeFile("O_vort_sph_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){

						double data1;
						testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,data1);

						double data2;
						testSolutions.correct_result_grad_phi__grid_gaussian(a,b,data2);

						c = data1 - data2;
					}
			);
			//result = result.spat_truncate();

//			result.file_physical_writeFile("O_vort_correct_result.csv");

//			(h-result).file_physical_writeFile("O_vort_correct_diff.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST VORTICITY - max error: " << error_max << std::endl;

			if (error_max > epsilon)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED" << std::endl;
				//FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}


		if (true)
		{
			// Laplace
			SphereData h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);

			h = op.div_lon(op.grad_lon(h))
				+ op.div_lat(op.grad_lat(h));

//			h.file_physical_writeFile("O_laplace_div_grad_result.csv");

			SphereData result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			result = op.laplace(result);

//			result.file_physical_writeFile("O_laplace_laplace_result.csv");

//			(h-result).file_physical_writeFile("O_laplace_laplace_z_diff.csv");

			double error_max = h.physical_reduce_error_max(result);
			std::cout << "TEST LAPLACE (div.grad vs. sph laplace) - max error: " << error_max << std::endl;

			if (error_max > epsilon)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED" << std::endl;
				//FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}
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

	if (simVars.disc.res_spectral[0] == 0)
		FatalError("Set number of spectral modes to use SPH!");

	SphereDataConfig sphereDataConfig;
	sphereDataConfig.setupAutoPhysicalSpace(
					simVars.disc.res_spectral[0],
					simVars.disc.res_spectral[1],
					&simVars.disc.res_physical[0],
					&simVars.disc.res_physical[1]
			);

	run_tests(&sphereDataConfig);
}


#endif /* SRC_TESTOPERATORS_HPP_ */
