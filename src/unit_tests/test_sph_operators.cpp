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
#include <sweet/sphere/ErrorCheck.hpp>



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

		SphereData zero(sphereDataConfig);
		zero.physical_set_zero();

		for (int i = 0; i < 3; i++)
		{
			double advection_rotation_angle = alpha[i];

			std::cout << "****************************************" << std::endl;
			std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;

			SphereData phi(sphereDataConfig);
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

			SphereData u(sphereDataConfig);
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

			SphereData v(sphereDataConfig);
			v.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
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
				SphereData lhs = op.robert_div(u, v);

				ErrorCheck::check(lhs, zero, "div freeness Robert formulation", error_threshold);
			}

			{
				// LHS = div(U*phi)
				SphereData lhs = op.robert_div(u*phi, v*phi);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				SphereData rhs = op.robert_grad_M(phi, u, v);

				ErrorCheck::check(lhs, rhs, "TEST div(U*phi) with Robert formulation", error_threshold*10e4);
			}
		}
	}
#endif



#if 1
	if (true)
	{
		std::cout << "****************************************" << std::endl;
		std::cout << "* DIV FREE TESTS NON-ROBERT *" << std::endl;
		std::cout << "****************************************" << std::endl;


		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0;
		double a = 6.37122e6;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		double alpha[] = {0, M_PI/3, M_PI/2};

		SphereData zero(sphereDataConfig);
		zero.physical_set_zero();

		for (int i = 0; i < 3; i++)
		{
			double advection_rotation_angle = alpha[i];

			std::cout << "****************************************" << std::endl;
			std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;

			SphereData phi(sphereDataConfig);
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

			SphereData u(sphereDataConfig);
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

			SphereData v(sphereDataConfig);
			v.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
					double i_phi = i_lat;
					double i_lambda = i_lon;
					io_data =
						-u0*(
								std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);
				}
			);

			{
				SphereData lhs = op.div(u, v);

				ErrorCheck::check(lhs, zero, "div freeness NON-Robert formulation", error_threshold, true /* Ignored, because NON-Robert */);
			}

			{
				// LHS = div(U*phi)
				SphereData lhs = op.div(u*phi, v*phi);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				SphereData rhs = op.grad_lon(phi)*u + op.grad_lat(phi)*v;

				ErrorCheck::check(lhs, rhs, "DIV(U*phi) NON-Robert", error_threshold, true /* Ignored, because NON-Robert */);
			}
		}
	}
#endif


#if 1
	if (true)
	{
		// 1.0/sqrt(1.0-mu^2)
		SphereData h(sphereDataConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &ret)
				{
					ret = (1.0-mu*mu);;
				}
		);

		h = op.inv_one_minus_mu2(h);

		h.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &ret)
				{
					//ret *= (1.0-mu*mu);
				}
		);

		SphereData one(sphereDataConfig);
		one.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &ret)
				{
					ret = 1.0;
				}
		);

		ErrorCheck::check(h, one, "mu^2 * 1.0/(1-mu^2)", error_threshold, false, true);
	}
#endif


#if 1
	if (true)
	{
		// 1.0/sqrt(1.0-mu^2)
		SphereData h(sphereDataConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &ret)
				{
					ret = (1.0-mu*mu);
				}
		);
		h = op.inv_one_minus_mu2(h);

		SphereData one(sphereDataConfig);
		one.physical_set_all_value(1.0);


		ErrorCheck::check(h, one, "1.0/(1-mu^2) ver2", error_threshold, false, true);
	}

#endif


#if 1
	if (true)
	{
		double u0 = 123;
		double advection_rotation_angle = M_PI*0.3;

		SphereData u(sphereDataConfig);
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

				io_data *= std::cos(i_lat);
			}
		);

		SphereData v(sphereDataConfig);
		v.physical_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_phi = i_lat;
				double i_lambda = i_lon;
				io_data =
					-u0*(
							std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);

				io_data *= std::cos(i_lat);
			}
		);

		SphereData lhs = op.robert_div(u, v);

		SphereData zero(sphereDataConfig);
		zero.spectral_set_zero();

		ErrorCheck::check(lhs, zero, "TEST div(u,v)=0 with Robert formulation", error_threshold*1e3);
	}
#endif

#if 0
	if (true)
	{
		double u0 = 123;
		double advection_rotation_angle = M_PI*0.3;

		SphereData u(sphereDataConfig);
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

		SphereData v(sphereDataConfig);
		v.physical_update_lambda(
			[&](double i_lon, double i_lat, double &io_data)
			{
				double i_phi = i_lat;
				double i_lambda = i_lon;
				io_data =
					-u0*(
							std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);
			}
		);

		SphereData lhs = op.div(u, v);

		SphereData zero(sphereDataConfig);
		zero.spectral_set_zero();

		ErrorCheck::check(lhs, zero, "TEST div(u,v)=0 with non-Robert formulation", error_threshold);
	}
#endif

#if 1
	if (true)
	{
		SphereData lhs(sphereDataConfig);
		lhs.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	  std::cos(i_theta)*std::cos(i_theta)*std::cos(i_theta)
							* std::cos(i_lambda)*std::sin(i_lambda);
			}
		);

		lhs = op.diff_lon(lhs);


		SphereData rhs(sphereDataConfig);
		rhs.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	  std::cos(i_theta)*std::cos(i_theta)*std::cos(i_theta)
							* std::cos(2*i_lambda);
			}
		);

		ErrorCheck::check(lhs, rhs, "TEST diff_lon(phi)", error_threshold*1e1);
	}
#endif


#if 0
	if (true)
	{
		SphereData phi(sphereDataConfig);
		phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data = std::cos(i_theta)*std::cos(i_theta)*std::cos(i_theta) * std::cos(i_lambda);
			}
		);

		SphereData div_lon_phi(sphereDataConfig);
		div_lon_phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data = -std::cos(i_theta)*std::cos(i_theta)*std::sin(i_lambda);
			}
		);

		SphereData lhs = op.div_lon(phi);
		SphereData rhs = div_lon_phi;

		ErrorCheck::check(lhs, rhs, "TEST div_i(phi) with NON Robert formulation", error_threshold);
	}
#endif


#if 1

	if (true)
	{
		SphereData lhs(sphereDataConfig);
		lhs.physical_update_lambda_cogaussian_grid(
			[&](double i_lambda, double i_comu, double &io_data)
			{
				io_data = std::sin(i_lambda)*i_comu*i_comu*i_comu;
			}
		);
		lhs = op.robert_div_lon(lhs);


		SphereData div_lon_phi(sphereDataConfig);
		div_lon_phi.physical_update_lambda_cogaussian_grid(
			[&](double i_lambda, double i_comu, double &io_data)
			{
				io_data = std::cos(i_lambda)*i_comu;
			}
		);

		SphereData rhs = div_lon_phi;

		ErrorCheck::check(lhs, rhs, "TEST div_i(phi) with Robert formulation", error_threshold);
	}

	if (true)
	{
		SphereData phi(sphereDataConfig);
		phi.physical_update_lambda_gaussian_grid(
			[&](double i_lambda, double i_mu, double &io_data)
			{
				io_data = (1.0-i_mu*i_mu);
			}
		);

		SphereData div_lat_phi(sphereDataConfig);
		div_lat_phi.physical_update_lambda_gaussian_grid(
			[&](double i_lambda, double i_mu, double &io_data)
			{
				io_data = (-2.0*i_mu);
			}
		);

		SphereData lhs = op.robert_div_lat(phi);
		SphereData rhs = div_lat_phi;

		ErrorCheck::check(lhs, rhs, "TEST div_j(phi) with Robert formulation", error_threshold);
	}
#endif


#if 1
	if (true)
	{
		double lambda_c = 3.0*M_PI/2.0;
//		double theta_c = 0;
		double theta_c = M_PI/3.0;
		double a = 6.37122e6;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		double alpha[] = {0, 0.0000001, M_PI/3, M_PI/2};

		for (int i = 0; i < 4; i++)
		{
			double advection_rotation_angle = alpha[i];

			std::cout << std::endl;
			std::cout << "Using rotation angle " << advection_rotation_angle << std::endl;

			SphereData phi(sphereDataConfig);
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


			SphereData u(sphereDataConfig);
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

					io_data *= std::cos(i_lat);
				}
			);

			SphereData v(sphereDataConfig);
			v.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
					double i_theta = i_lat;
					double i_lambda = i_lon;
					io_data =
						-u0*(
								std::sin(i_lambda)*std::sin(advection_rotation_angle)
						);

					io_data *= std::cos(i_lat);
				}
			);

			SphereData zero(sphereDataConfig);
			zero.physical_set_zero();

			{
				SphereData ret = op.robert_div(u, v);

				ErrorCheck::check(ret, zero, "TEST div freeness", error_threshold*10e2);
			}

#if 0
			{
	#if 1
				// LHS = div(U*phi)
				SphereData lhs = op.robert_div(u*phi, v*phi)*(1.0/a);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				//SphereData rhs = (op.robert_grad_lon_M(phi)*u + op.robert_grad_lat_M(phi)*v)*(1.0/a);
				SphereData rhs = op.robert_grad_M(phi, u ,v)*(1.0/a);

	#else
				u = zero;
				v = zero;
				SphereData lhs = op.robert_div(u*phi, v*phi);
				SphereData rhs = op.robert_grad_lon_M(phi)*u + op.robert_grad_lat_M(phi)*v;
	#endif

				ErrorCheck::check(lhs, rhs, "TEST div(U*phi) - grad(phi)*U with Robert formulation", error_threshold);
			}
#endif
		}
	}
#endif



	if (true)
	{
		SphereTestSolutions_SPH testSolutionsSph(2,1);

		//int tn = 1;
		//int tm = 0;
		SphereData h(sphereDataConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double a, double b, double &c){testSolutionsSph.test_function__grid_gaussian(a,b,c);}
		);
		h.request_data_spectral();


		SphereData result(sphereDataConfig);
		result.physical_update_lambda_gaussian_grid(
				[&](double a, double b, double &c){testSolutionsSph.correct_result_diff_mu__grid_gaussian(a,b,c);}
		);
		result.request_data_spectral();

		// one_minus_mu_squared_diff_lat
		h = op.spectral_one_minus_mu_squared_diff_lat_mu(h);

		double error_max = h.physical_reduce_max_abs(result);
		std::cout << "TEST SPHbasis (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

		if (error_max > error_threshold)
			FatalError("ERROR THRESHOLD EXCEEDED!");
	}

#if 1
	if (true)
	{
		SphereTestSolutions_Gaussian testSolutions;

		SphereData h(sphereDataConfig);
		h.physical_update_lambda_gaussian_grid(
				[&](double lambda, double mu, double &c)
				{
					double phi = std::asin(mu);

					testSolutions.test_function__grid_gaussian(lambda,mu,c);
					c *= cos(phi)*cos(phi);

					c *= cos(phi);
				}
		);
		h = op.robert_div_lon(h);

		SphereData result(sphereDataConfig);
		result.physical_update_lambda_gaussian_grid(
		[&](double lambda, double mu, double &c)
				{
					double phi = std::asin(mu);

					testSolutions.correct_result_diff_lambda__grid_gaussian(lambda, mu, c);

					c *= cos(phi);
				}
		);

		ErrorCheck::check(h, result, "ROBERT DIV LONGITUDE", error_threshold*1e2);
	}
#endif


	{
		SphereTestSolutions_Gaussian testSolutions;


		if (true)
		{
			SphereData data(sphereDataConfig);
			data.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y;
					}
			);

			data = data*123.0;

			SphereData data2(sphereDataConfig);
			data2.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y*123.0;
					}
			);

			ErrorCheck::check(data, data2, "operator*123.0");
		}

		if (true)
		{
			SphereData data(sphereDataConfig);
			data.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y;
					}
			);

			data = data + 100.0;

			SphereData data2(sphereDataConfig);
			data2.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y+100.0;
					}
			);

			ErrorCheck::check(data, data2, "operator+(double)");
		}

		if (true)
		{
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);

			SphereData hphi(sphereDataConfig);
			hphi.physical_update_lambda(
					[&](double a, double b, double &c){testSolutions.test_function_phi__grid_phi(a,b,c);}
			);

			ErrorCheck::check(h, hphi, " vs. MU");
		}

		if (true)
		{
			// identity
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "TEST IDENTITY");
		}


		if (true)
		{
			// d/d lambda
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);

			h = op.diff_lon(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_lambda__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "TEST DIFF LON");
		}


		if (true)
		{
			// d/d phi
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.diff_lat_phi(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_phi__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "TEST DIFF PHI");
		}


		if (true)
		{
			// d/d mu
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double mu, double &c){c = mu;}
			);
			h = op.diff_lat_mu(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double mu, double &c){c = 1.0;}
			);


			ErrorCheck::check(h, result, "TEST DIFF LAT MU", error_threshold*1e2, false, true);
		}

		if (true)
		{
			// mu*F(\lambda,\mu)
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.mu(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_mu__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "TEST mu*()");
		}

		if (true)
		{
			// mu*mu*F(\lambda,\mu)
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.mu2(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &i_data){
						testSolutions.test_function__grid_gaussian(lat, mu, i_data);
						i_data *= mu*mu;
					}
			);

			ErrorCheck::check(h, result, "TEST mu*mu*()");
		}


		if (true)
		{
			// one_minus_mu_squared_diff_lat
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.spectral_one_minus_mu_squared_diff_lat_mu(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "TEST (1-mu*mu)*d/dmu");
		}


		if (true)
		{
			// grad lambda
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.grad_lon(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "TEST GRAD LON");
		}


		if (true)
		{
			// grad mu
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.grad_lat(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_grad_phi__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "LAT");
		}


		if (true)
		{
			// div lambda
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.div_lon(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_div_lambda__grid_gaussian(a,b,c);}
			);

			ErrorCheck::check(h, result, "LAT");
		}


		if (true)
		{
			// div mu
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){
						testSolutions.test_function__grid_gaussian(a,b,c);
					}
			);

			h = op.div_lat(h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){
						testSolutions.correct_result_div_mu__grid_gaussian(a,b,c);
					}
			);

			double error_max = h.physical_reduce_max_abs(result);
			std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED since DIV should be only applied to a vector field, however we do this on a scalar field" << std::endl;
//				FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}


#if 0
		if (true)
		{
			// div mu TEST
			SphereData h(sphereDataConfig);
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
			SphereData result(sphereDataConfig);
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

			double error_max = h.physical_reduce_max_abs(result);
			std::cout << "TEST DIV TEST LAT  - max error: " << error_max << std::endl;
		}
#endif


		if (true)
		{
			// divergence
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.div(h, h);
//				h.spat_update_lambda_gaussian_grid([&](double a, double mu, double &c){c *= mu*mu*mu*mu;});
//			h.file_physical_writeFile("O_divergence_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){

						double data1;
						testSolutions.correct_result_div_mu__grid_gaussian(a,b,data1);

						double data2;
						testSolutions.correct_result_div_lambda__grid_gaussian(a,b,data2);

						c = data1 + data2;
					}
			);

			double error_max = h.physical_reduce_max_abs(result);
			std::cout << "TEST DIVERGENCE - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED since DIV should be only applied to a vector field, however we do this on a scalar field" << std::endl;
				//FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}


		if (true)
		{
			// vorticity
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.vort(h, h);

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double lambda, double mu, double &c){

						double data1;
						testSolutions.correct_result_div_lambda__grid_gaussian(lambda,mu,data1);

						double data2;
						testSolutions.correct_result_div_mu__grid_gaussian(lambda,mu,data2);

						c = data1 - data2;
					}
			);

			double error_max = h.physical_reduce_max_abs(result);
			std::cout << "TEST VORTICITY - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED since DIV should be only applied to a vector field, however we do this on a scalar field" << std::endl;
				//FatalError("ERROR THRESHOLD EXCEEDED!");
			}
		}


		if (true)
		{
			// Laplace
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);

			h = op.div_lon(op.grad_lon(h))
				+ op.div_lat(op.grad_lat(h));

//			h.file_physical_writeFile("O_laplace_div_grad_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			result = op.laplace(result);

//			result.file_physical_writeFile("O_laplace_laplace_result.csv");

//			(h-result).file_physical_writeFile("O_laplace_laplace_z_diff.csv");

			double error_max = h.physical_reduce_max_abs(result);
			std::cout << "TEST LAPLACE (div.grad vs. sph laplace) - max error: " << error_max << std::endl;

			if (error_max > error_threshold)
			{
				std::cerr << "EXCEEDED ERROR THRESHOLD IGNORED since DIV should be only applied to a vector field, however we do this on a scalar field" << std::endl;
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
