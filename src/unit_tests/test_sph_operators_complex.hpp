/*
 * AppTestSPHOperatorsComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTOPERATORS_COMPLEX_HPP_
#define SRC_TESTOPERATORS_COMPLEX_HPP_

#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperators.hpp>



class AppTestSPHOperatorsComplex
{
public:
	SimulationVariables simVars;
	SphereOperatorsComplex opComplex;


	void run(
			SphereDataConfig *i_sphConfig
	)
	{
//		opComplex.setup(i_sphConfig);

		double epsilon = 1e-10;

		if (true)
		{
			SphereTestSolutions_SPH testSolutionsSph(2,1);

			int tn = 1;
			int tm = 0;
			SphereDataComplex h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutionsSph.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h.request_data_spectral();


			SphereDataComplex result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutionsSph.correct_result_diff_mu__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			result.request_data_spectral();

			h.spat_write_file("O_SPHbasis_test_function.csv");


			// one_minus_mu_squared_diff_lat
			h = opComplex.spec_one_minus_mu_squared_diff_lat_mu(h);
			h.spat_write_file("O_SPHbasis_one_minus_mu_squared_diff_mu_sph_result.csv");

			result.spat_write_file("O_SPHbasis_one_minus_mu_squared_diff_mu_correct_result.csv");


			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST SPHbasis (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;
		}


		{
			SphereTestSolutions_Gaussian testSolutions;



			if (true)
			{
				std::complex<double> scalar = {12.0, 34.0};
				SphereDataComplex data(i_sphConfig);
				data.spat_update_lambda(
						[&](double x, double y, cplx &io_data)
						{
							io_data = y;
						}
				);

				data = data*scalar;
				data.spec_truncate();

				SphereDataComplex data2(i_sphConfig);
				data2.spat_update_lambda(
						[&](double x, double y, cplx &io_data)
						{
							io_data = y*scalar;
						}
				);
				data2.spec_truncate();

				double error = data.spat_reduce_error_max(data2);
				std::cout << "ERROR operator*" << scalar << ": " << error << std::endl;
			}

			if (true)
			{
				SphereDataComplex data(i_sphConfig);
				data.spat_update_lambda(
						[&](double x, double y, std::complex<double> &io_data)
						{
							io_data = y;
						}
				);

				data = data + 100.0;
				data.spec_truncate();

				SphereDataComplex data2(i_sphConfig);
				data2.spat_update_lambda(
						[&](double x, double y, std::complex<double> &io_data)
						{
							io_data = y+100.0;
						}
				);
				data2.spec_truncate();

				double error = data.spat_reduce_error_max(data2);
				std::cout << "ERROR operator+(double): " << error << std::endl;
			}

			if (true)
			{
				std::complex<double> offset(200.0, 100.0);

				SphereDataComplex data(i_sphConfig);
				data.spat_update_lambda(
						[&](double x, double y, std::complex<double> &io_data)
						{
							io_data.real(x*y);
							io_data.imag(x+y);
						}
				);

				data = data+offset;
				data.spec_truncate();


				SphereDataComplex data2(i_sphConfig);
				data2.spat_update_lambda(
						[&](double x, double y, std::complex<double> &io_data)
						{
							io_data.real(x*y+offset.real());
							io_data.imag(x+y+offset.imag());
						}
				);
				data2.spec_truncate();

				double error = data.spat_reduce_error_max(data2);
				std::cout << "ERROR operator+(cplx): " << error << std::endl;
			}



			if (true)
			{
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h.spat_write_file("O_test_function.csv");

				SphereDataComplex hphi(i_sphConfig);
				hphi.spat_update_lambda(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function_phi__grid_phi(a,b,tmp);
							c = tmp;
						}
				);
				hphi.spat_write_file("O_test_function_phi.csv");

				double error_max = h.spat_reduce_error_max(hphi);
				std::cout << "TEST PHI vs. MU: max error: " << error_max << std::endl;
			}

			if (true)
			{
				// identity
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//h = h.spat_truncate();

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//result = result.spat_truncate();

				double error_max = h.spat_reduce_error_max(result);

				std::cout << "TEST IDENTITY: max error: " << error_max << std::endl;
			}


			if (true)
			{
				// d/d lambda
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//h = h.spat_truncate();

				h = opComplex.diff_lon(h);

				h.spat_write_file("O_diff_lambda_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_diff_lambda__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_diff_lambda_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIFF LON - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// d/d phi
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.diff_lat_phi(h);
				h.spat_write_file("O_diff_phi_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_diff_phi__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_diff_phi_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIFF PHI - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// d/d mu
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.diff_lat_mu(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_diff_mu_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_diff_mu__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_diff_mu_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIFF LAT MU - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// grad lambda
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.grad_lon(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_grad_lambda_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_grad_lambda_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST GRAD LON - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// mu*F(\lambda,\mu)
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.mu(h);
				h.spat_write_file("O_mu_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_mu__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				result.spat_write_file("O_mu_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST mu*() - max error: " << error_max << std::endl;
			}

			if (true)
			{
				// mu*mu*F(\lambda,\mu)
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.mu2(h);
				h.spat_write_file("O_mu2_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &i_data){
							double tmp;
							testSolutions.test_function__grid_gaussian(lat, mu, tmp);
							i_data = tmp;

							i_data *= mu*mu;
						}
				);
				result.spat_write_file("O_mu2_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST mu*mu*() - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// one_minus_mu_squared_diff_lat
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.spec_one_minus_mu_squared_diff_lat_mu(h);
				h.spat_write_file("O_one_minus_mu_squared_diff_mu_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);

				result.spat_write_file("O_one_minus_mu_squared_diff_mu_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// grad mu
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.grad_lat(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_grad_phi_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp = 0;
							testSolutions.correct_result_grad_phi__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);

				result.spat_write_file("O_grad_phi_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST GRAD LAT - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// div lambda
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.div_lon(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_div_lambda_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_div_lambda__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_div_lambda_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIV LON  - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// div mu
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				h = opComplex.div_lat(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_div_mu_sph_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.correct_result_div_mu__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				//result = result.spat_truncate();

				result.spat_write_file("O_div_mu_correct_result.csv");

				(h-result).spat_write_file("O_div_mu_correct_diff.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// Laplace
				SphereDataComplex h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);

				h = opComplex.div_lon(opComplex.grad_lon(h))
					+ opComplex.div_lat(opComplex.grad_lat(h));

				//h = opComplex.mu2(h);
				//h = opComplex.mu2(h);

				h.spat_write_file("O_laplace_div_grad_result.csv");

				SphereDataComplex result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, std::complex<double> &c)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(a,b,tmp);
							c = tmp;
						}
				);
				result = opComplex.laplace(result);

				//result = opComplex.mu2(result);
				//result = opComplex.mu2(result);

				result.spat_write_file("O_laplace_laplace_result.csv");

				(h-result).spat_write_file("O_laplace_laplace_z_diff.csv");

				result.spec_truncate();
				h.spec_truncate();

				double error_max = (h-result).spat_reduce_error_max();
				std::cout << "TEST LAPLACE  - max error: " << error_max << std::endl;
			}
		}
	}
};


#endif /* SRC_TESTOPERATORS_HPP_ */
