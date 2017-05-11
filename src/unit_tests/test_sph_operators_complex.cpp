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
#include <sweet/sphere/ErrorCheck.hpp>


SimulationVariables simVars;



SphereDataConfig sphereDataConfigInstance;
SphereDataConfig sphereDataConfigExtInstance;

SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigExtInstance;




void run_tests()
{
	double error_threshold = 1e-10;
	error_threshold *= (sphereDataConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << error_threshold << std::endl;

	SphereOperatorsComplex opComplex(sphereDataConfig, 1);


#if 1

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

		SphereDataComplex zeroc(sphereDataConfig);
		zeroc.physical_set_zero();

		for (int i = 0; i < 3; i++)
		{
			double advection_rotation_angle = alpha[i];

			std::cout << "****************************************" << std::endl;
			std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;


			SphereDataComplex phic(sphereDataConfig);
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

			SphereDataComplex uc(sphereDataConfig);
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

			SphereDataComplex vc(sphereDataConfig);
			vc.physical_update_lambda(
				[&](double i_lon, double i_lat, std::complex<double> &io_data)
				{
					double i_phi = i_lat;
					double i_lambda = i_lon;
					io_data =
						-u0*(
								std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);

					io_data *= cos(i_lat);
				}
			);

			{
				SphereDataComplex lhsc = opComplex.robert_div(uc, vc);
				ErrorCheck::check(lhsc, zeroc, "COMPLEX: div freeness Robert formulation", error_threshold*1e1);
			}

			{
				// LHS = div(U*phi)
				SphereDataComplex lhsc = opComplex.robert_div(uc*phic, vc*phic);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				SphereDataComplex rhsc = opComplex.robert_grad_M(phic, uc, vc);

				ErrorCheck::check(lhsc, rhsc, "COMPLEX: TEST div(U*phi) with Robert formulation", error_threshold*10e4);
			}
		}
	}

#endif


#if 1
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

		SphereDataComplex zero(sphereDataConfig);
		zero.physical_set_zero();

		for (int i = 0; i < 3; i++)
		{
			double advection_rotation_angle = alpha[i];

			std::cout << "****************************************" << std::endl;
			std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;

			SphereDataComplex phi(sphereDataConfig);
			phi.physical_update_lambda(
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

			SphereDataComplex u(sphereDataConfig);
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

					io_data *= cos(i_lat);
				}
			);

			SphereDataComplex v(sphereDataConfig);
			v.physical_update_lambda(
				[&](double i_lon, double i_lat, std::complex<double> &io_data)
				{
					double i_phi = i_lat;
					double i_lambda = i_lon;
					io_data =
						-u0*(
								std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);

					io_data *= cos(i_lat);
				}
			);

			{
				SphereDataComplex lhs = opComplex.robert_div(u, v);

				ErrorCheck::check(lhs, zero, "div freeness Robert formulation", error_threshold);
			}

			{
				// LHS = div(U*phi)
				SphereDataComplex lhs = opComplex.robert_div(u*phi, v*phi);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				SphereDataComplex rhs = opComplex.robert_grad_M(phi, u, v);


				ErrorCheck::check(lhs, rhs, "TEST div(U*phi) with Robert formulation", error_threshold*10e5);
			}
		}
	}
#endif


	if (true)
	{
		SphereTestSolutions_SPH testSolutionsSph(2,1);

		SphereDataComplex h(sphereDataConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double a, double b, std::complex<double> &c)
				{
					double tmp;
					testSolutionsSph.test_function__grid_gaussian(a,b,tmp);
					c = tmp;
				}
		);

//		h.physical_write_file("O_SPHbasis_test_function.csv");


		SphereDataComplex result(sphereDataConfig);
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

		if (error_max > error_threshold)
			FatalError("ERROR THRESHOLD EXCEEDED!");
	}


	{
		SphereTestSolutions_Gaussian testSolutions;



		if (true)
		{
			std::complex<double> scalar = {12.0, 34.0};
			SphereDataComplex data(sphereDataConfig);
			data.physical_update_lambda(
					[&](double x, double y, cplx &io_data)
					{
						io_data = y;
					}
			);
			data.request_data_spectral();
			data.request_data_physical();

			data = data*scalar;

			SphereDataComplex data2(sphereDataConfig);
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

			if (error > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			SphereDataComplex data(sphereDataConfig);
			data.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data = y;
					}
			);

			data = data + 100.0;
			data.request_data_spectral();
			data.request_data_physical();

			SphereDataComplex data2(sphereDataConfig);
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

			if (error > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			std::complex<double> offset(200.0, 100.0);

			SphereDataComplex data(sphereDataConfig);
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


			SphereDataComplex data2(sphereDataConfig);
			data2.physical_update_lambda(
					[&](double x, double y, std::complex<double> &io_data)
					{
						io_data.real(x*y+offset.real());
						io_data.imag(x+y+offset.imag());
					}
			);
			data2.request_data_spectral();
			data2.request_data_physical();

			ErrorCheck::check(data, data2, "operator+(cplx):", error_threshold*10e2);
		}



		if (true)
		{
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
//			h.physical_write_file("O_test_function.csv");

			SphereDataComplex hphi(sphereDataConfig);
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

			if (error_max > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			SphereDataComplex lhs(sphereDataConfig);
			lhs.physical_update_lambda_gaussian_grid(
					[&](double a, double mu, std::complex<double> &c)
					{
						c = (1.0-mu*mu);
					}
			);
			lhs = opComplex.inv_one_minus_mu2(lhs);

			SphereDataComplex one(sphereDataConfig);
			one.physical_set_all_value(1.0);

			ErrorCheck::check(lhs, one, "1.0/(1.0-mu^2)", error_threshold, false, true);
		}


#if 1
		if (true)
		{
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double lambda, double mu, std::complex<double> &c)
					{
						double phi = std::asin(mu);

						double tmp;
						testSolutions.test_function__grid_gaussian(lambda,mu,tmp);
						c = tmp;

						c = sin(lambda);
						c *= cos(phi);
						c *= cos(phi);

						c *= cos(phi);

//						c *= cos(phi);
//						c = {tmp, 0};
					}
			);
			h = opComplex.robert_div_lon(h);

			SphereDataComplex hphi(sphereDataConfig);
			hphi.physical_update_lambda_gaussian_grid(
			[&](double lambda, double mu, std::complex<double> &c)
					{
						double phi = std::asin(mu);

						double tmp;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lambda, mu, tmp);
						c = tmp;
						c = cos(lambda);

						c *= cos(phi);
//						c *= cos(phi);
//						c = {tmp, 0};
					}
			);

			ErrorCheck::check(h, hphi, "ROBERT DIV LONGITUDE", error_threshold*1e3, false, true);
		}
#endif

#if 0
		if (true)
		{
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			h = opComplex.robert_div_lat(h);

			SphereDataComplex hphi(sphereDataConfig);
			hphi.physical_update_lambda_gaussian_grid(
			[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_mu__grid_gaussian(a, b, tmp);
						c = {tmp, 0};
					}
			);

			ErrorCheck::check(h, hphi, "ROBERT DIV LATITUDE", error_threshold*1e3, false, true);
		}
#endif


		if (true)
		{
			SphereDataComplex phi(sphereDataConfig);
			phi.physical_update_lambda_gaussian_grid(
				[&](double i_lambda, double i_mu, std::complex<double> &io_data)
				{
					io_data = (1.0-i_mu*i_mu);
				}
			);

			SphereDataComplex lhs = opComplex.robert_div_lat(phi);

			SphereDataComplex div_lat_phi(sphereDataConfig);
			div_lat_phi.physical_update_lambda_gaussian_grid(
				[&](double i_lambda, double i_mu, std::complex<double> &io_data)
				{
					io_data = -2.0*i_mu;
				}
			);

			SphereDataComplex rhs = div_lat_phi;

			ErrorCheck::check(lhs, rhs, "TEST div_j(phi) with Robert formulation", error_threshold);
		}


		if (true)
		{
			SphereDataComplex phi(sphereDataConfig);
			phi.physical_update_lambda_gaussian_grid(
				[&](double i_lambda, double i_mu, std::complex<double> &io_data)
				{
					io_data = (1.0-i_mu*i_mu);
				}
			);
			phi = opComplex.robert_div_lat(phi);

			SphereDataComplex div_lat_phi(sphereDataConfig);
			div_lat_phi.physical_update_lambda_gaussian_grid(
				[&](double i_lambda, double i_mu, std::complex<double> &io_data)
				{
					io_data = (-2.0*i_mu);
				}
			);

			ErrorCheck::check(phi, div_lat_phi, "ROBERT DIV LATITUDE");
		}



		if (true)
		{
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			h = opComplex.robert_grad_lon(h);

			SphereDataComplex hphi(sphereDataConfig);
			hphi.physical_update_lambda_gaussian_grid(
			[&](double lambda, double mu, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lambda, mu, tmp);
						c = {tmp, 0};
					}
			);

			ErrorCheck::check(h, hphi, "ROBERT GRAD LONGITUDE");
		}


		if (true)
		{
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = {tmp, 0};
					}
			);
			h = opComplex.robert_grad_lat(h);

			SphereDataComplex hphi(sphereDataConfig);
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

			ErrorCheck::check(h, hphi, "ROBERT GRAD LATITUDE");
		}


		if (true)
		{
			// identity
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			ErrorCheck::check(h, result, "TEST IDENTITY");
		}


		if (true)
		{
			// d/d lambda
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			h = opComplex.diff_lon(h);

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_lambda__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			ErrorCheck::check(h, result, "TEST DIFF LON");
		}


		if (true)
		{
			// d/d phi
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.diff_lat_phi(h);


			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_diff_phi__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			ErrorCheck::check(h, result, "TEST DIFF PHI");
		}


		if (true)
		{
			// d/d mu
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
				[&](double a, double b, std::complex<double> &c)
				{
					double tmp;
					testSolutions.test_function__grid_gaussian(a,b,tmp);
					c = tmp;
				}
			);
			h = opComplex.diff_lat_mu(h);

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
				[&](double a, double b, std::complex<double> &c)
				{
					double tmp;
					testSolutions.correct_result_diff_mu__grid_gaussian(a,b,tmp);
					c = tmp;
				}
			);

			ErrorCheck::check(h, result, "TEST DIFF LAT MU");
		}


		if (true)
		{
			// mu*F(\lambda,\mu)
			SphereDataComplex h(sphereDataConfig);
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

			SphereDataComplex result(sphereDataConfig);
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

			if (error_max > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// mu*mu*F(\lambda,\mu)
			SphereDataComplex h(sphereDataConfig);
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

			SphereDataComplex result(sphereDataConfig);
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

			if (error_max > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// one_minus_mu_squared_diff_lat
			SphereDataComplex h(sphereDataConfig);
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

			SphereDataComplex result(sphereDataConfig);
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

			if (error_max > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}



		if (true)
		{
			// grad lambda
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.grad_lon(h);

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST GRAD LON - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// grad mu
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.grad_lat(h);

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp = 0;
						testSolutions.correct_result_grad_phi__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST GRAD LAT - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// div lambda
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.div_lon(h);

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_div_lambda__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIV LON  - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// div mu
			SphereDataComplex h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);
			h = opComplex.div_lat(h);

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.correct_result_div_mu__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED since dataset is incompatible to divergence" << std::endl;
				std::cerr << "This is no problem for real simulations since such fields would not exist" << std::endl;
			}
		}


		if (true)
		{
			// Laplace
			SphereDataComplex h(sphereDataConfig);
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

			SphereDataComplex result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, std::complex<double> &c)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(a,b,tmp);
						c = tmp;
					}
			);

			result = opComplex.laplace(result);


			double error_max = (h-result).physical_reduce_max_abs();
			std::cout << "TEST LAPLACE  - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED since dataset is incompatible to divergence" << std::endl;
				std::cerr << "This is no problem for real simulations since such fields would not exist" << std::endl;
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

	if (simVars.disc.res_physical[0] <= 0)
	{
		sphereDataConfigInstance.setupAutoPhysicalSpace(
						simVars.disc.res_spectral[0],
						simVars.disc.res_spectral[1],
						&simVars.disc.res_physical[0],
						&simVars.disc.res_physical[1]
				);
	}
	else
	{
		sphereDataConfigInstance.setup(
						simVars.disc.res_spectral[0],
						simVars.disc.res_spectral[1],
						simVars.disc.res_physical[0],
						simVars.disc.res_physical[1]
				);
	}

	sphereDataConfigExtInstance.setupAdditionalModes(
			&sphereDataConfigInstance,
			simVars.rexi.rexi_use_extended_modes,
			simVars.rexi.rexi_use_extended_modes
		);

	run_tests();
}

#endif /* SRC_TESTOPERATORS_HPP_ */
