/*
 * AppTestSPHOperatorsComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_TESTOPERATORS_COMPLEX_HPP_
#define SRC_TESTOPERATORS_COMPLEX_HPP_

#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>


SimulationVariables simVars;
SphereOperatorsComplex opComplex;

void run_tests(
		SphereDataConfig *sphConfig
)
{
	double epsilon = 1e-11;
	epsilon *= (sphConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;


	if (true)
	{
		SphereTestSolutions_SPH testSolutionsSph(2,1);

		SphereDataComplex h(sphConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double a, double b, std::complex<double> &c)
				{
					double tmp;
					testSolutionsSph.test_function__grid_gaussian(a,b,tmp);
					c = tmp;
				}
		);

//		h.physical_write_file("O_SPHbasis_test_function.csv");


		SphereDataComplex result(sphConfig);
		result.physical_update_lambda_gaussian_grid(
				[&](double a, double b, std::complex<double> &c)
				{
					double tmp;
					testSolutionsSph.correct_result_diff_mu__grid_gaussian(a,b,tmp);
					c = tmp;
				}
		);


		// one_minus_mu_squared_diff_lat
		h = opComplex.spectral_one_minus_mu_squared_diff_lat_mu(h);

		double error_max = h.physical_reduce_max(result);
		std::cout << "TEST SPHbasis (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

		if (error_max > epsilon)
			FatalError("ERROR THRESHOLD EXCEEDED!");
	}


	{
		SphereTestSolutions_Gaussian testSolutions;



		if (true)
		{
			std::complex<double> scalar = {12.0, 34.0};
			SphereDataComplex data(sphConfig);
			data.physical_update_lambda(
					[&](double x, double y, cplx &io_data)
					{
						io_data = y;
					}
			);
			data.request_data_spectral();
			data.request_data_physical();

			data = data*scalar;

			SphereDataComplex data2(sphConfig);
			data2.physical_update_lambda(
					[&](double x, double y, cplx &io_data)
					{
						io_data = y*scalar;
					}
			);
			data2.request_data_spectral();
			data2.request_data_physical();

			double error = data.physical_reduce_max(data2);
			std::cout << "ERROR operator*" << scalar << ": " << error << std::endl;

			if (error > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			SphereDataComplex data(sphConfig);
			data.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y;
					}
			);

			data = data + 100.0;
			data.request_data_spectral();
			data.request_data_physical();

			SphereDataComplex data2(sphConfig);
			data2.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y+100.0;
					}
			);
			data2.request_data_spectral();
			data2.request_data_physical();

			double error = data.physical_reduce_max(data2);
			std::cout << "ERROR operator+(double): " << error << std::endl;

			if (error > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			std::complex<double> offset(200.0, 100.0);

			SphereDataComplex data(sphConfig);
			data.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data.real(x*y);
						io_data.imag(x+y);
					}
			);

			data = data+offset;
			data.request_data_spectral();
			data.request_data_physical();


			SphereDataComplex data2(sphConfig);
			data2.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data.real(x*y+offset.real());
						io_data.imag(x+y+offset.imag());
					}
			);
			data2.request_data_spectral();
			data2.request_data_physical();

			double error = data.physical_reduce_max(data2);
			std::cout << "ERROR operator+(cplx): " << error << std::endl;

			if (error > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}



		if (true)
		{
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
//			h.physical_write_file("O_test_function.csv");

			SphereDataComplex hphi(sphConfig);
			hphi.physical_update_lambda(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function_phi__grid_phi(a,b,tmp);
						c = {tmp, 0};
					}
			);
//			hphi.physical_write_file("O_test_function_phi.csv");

			double error_max = h.physical_reduce_max(hphi);
			std::cout << "TEST PHI vs. MU: max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}



		if (true)
		{
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			h = opComplex.robert_div_lon(h);

			SphereDataComplex hphi(sphConfig);
			hphi.physical_update_lambda_gaussian_grid(
			[&](double lambda, double mu, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lambda, mu, tmp);
						double phi = asin(mu);
						tmp /= cos(phi)*cos(phi);
						c = {tmp, 0};
					}
			);

			double error_max = h.physical_reduce_max(hphi);
			std::cout << "ROBERT DIV LONGITUDE: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			h = opComplex.robert_div_lat(h);

			SphereDataComplex hphi(sphConfig);
			hphi.physical_update_lambda_gaussian_grid(
			[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_mu__grid_gaussian(a, b, tmp);
						c = {tmp, 0};
					}
			);

			double error_max = h.physical_reduce_max(hphi);
			std::cout << "ROBERT DIV LATITUDE: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}



		if (true)
		{
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			h = opComplex.robert_grad_lon(h);

			SphereDataComplex hphi(sphConfig);
			hphi.physical_update_lambda_gaussian_grid(
			[&](double lambda, double mu, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lambda, mu, tmp);
						c = {tmp, 0};
					}
			);

			double error_max = h.physical_reduce_max(hphi);
			std::cout << "ROBERT GRAD LONGITUDE: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			h = opComplex.robert_grad_lat(h);

			SphereDataComplex hphi(sphConfig);
			hphi.physical_update_lambda_gaussian_grid(
			[&](double lambda, double mu, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_mu__grid_gaussian(lambda, mu, tmp);
						double phi = asin(mu);
						tmp *= cos(phi)*cos(phi);
						c = {tmp, 0};
					}
			);

			double error_max = h.physical_reduce_max(hphi);
			std::cout << "ROBERT GRAD LATITUDE: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// identity
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = h.physical_truncate();

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();

			double error_max = h.physical_reduce_max(result);

			std::cout << "TEST IDENTITY: max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d lambda
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = h.physical_truncate();

			h = opComplex.diff_lon(h);

//			h.physical_write_file("O_diff_lambda_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_lambda__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();
//			result.physical_write_file("O_diff_lambda_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIFF LON - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d phi
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.diff_lat_phi(h);
			h = h.physical_truncate();
//			h.physical_write_file("O_diff_phi_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_phi__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();
//			result.physical_write_file("O_diff_phi_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIFF PHI - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d mu
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.diff_lat_mu(h);
			h = h.physical_truncate();
//			h.physical_write_file("O_diff_mu_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_mu__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();
//			result.physical_write_file("O_diff_mu_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIFF LAT MU - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// mu*F(\lambda,\mu)
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.mu(h);
//			h.physical_write_file("O_mu_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_mu__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
//			result.physical_write_file("O_mu_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST mu*() - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// mu*mu*F(\lambda,\mu)
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.mu2(h);
//			h.physical_write_file("O_mu2_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &i_data){
						double tmp;
						testSolutions.test_function__grid_gaussian(lat, mu, tmp);
						i_data = tmp;

						i_data *= mu*mu;
					}
			);
//			result.physical_write_file("O_mu2_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST mu*mu*() - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// one_minus_mu_squared_diff_lat
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.spectral_one_minus_mu_squared_diff_lat_mu(h);
//			h.physical_write_file("O_one_minus_mu_squared_diff_mu_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

//			result.physical_write_file("O_one_minus_mu_squared_diff_mu_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}



		if (true)
		{
			// grad lambda
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.grad_lon(h);
			h = h.physical_truncate();
//			h.physical_write_file("O_grad_lambda_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();
//			result.physical_write_file("O_grad_lambda_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST GRAD LON - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// grad mu
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.grad_lat(h);
			h = h.physical_truncate();
//			h.physical_write_file("O_grad_phi_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp = 0;
						testSolutions.correct_result_grad_phi__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();

//			result.physical_write_file("O_grad_phi_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST GRAD LAT - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// div lambda
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.div_lon(h);
			h = h.physical_truncate();

//			h.physical_write_file("O_div_lambda_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_div_lambda__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();
//			result.physical_write_file("O_div_lambda_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIV LON  - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// div mu
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.div_lat(h);
			h = h.physical_truncate();

//			h.physical_write_file("O_div_mu_sph_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_div_mu__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			result = result.physical_truncate();

//			result.physical_write_file("O_div_mu_correct_result.csv");

//			(h-result).physical_write_file("O_div_mu_correct_diff.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;

			if (error_max > epsilon)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED" << std::endl;
//				FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}


		if (true)
		{
			// Laplace
			SphereDataComplex h(sphConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			h = opComplex.div_lon(opComplex.grad_lon(h))
				+ opComplex.div_lat(opComplex.grad_lat(h));

			h.physical_truncate();

			//h = opComplex.mu2(h);
			//h = opComplex.mu2(h);

//			h.physical_write_file("O_laplace_div_grad_result.csv");

			SphereDataComplex result(sphConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			result = opComplex.laplace(result);

			result.physical_truncate();

			//result = opComplex.mu2(result);
			//result = opComplex.mu2(result);

//			result.physical_write_file("O_laplace_laplace_result.csv");

//			(h-result).physical_write_file("O_laplace_laplace_z_diff.csv");

			//result.physical_truncate();
			//h.physical_truncate();

			double error_max = (h-result).physical_reduce_max_abs();
			std::cout << "TEST LAPLACE  - max error: " << error_max << std::endl;

			if (error_max > epsilon)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED" << std::endl;
//				FatalError("ERROR THRESHOLD EXCEEDED!");
			}
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
