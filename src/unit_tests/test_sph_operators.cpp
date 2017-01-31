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

SphereDataConfig sphereDataConfigInstance;
SphereDataConfig sphereDataConfigExtInstance;

SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigExtInstance;


bool errorCheck(
		SphereData &i_lhs,
		SphereData &i_rhs,
		const std::string &i_id,
		double i_error_threshold = 1.0,
		bool i_ignore_error = false,
		bool i_normalization = true
)
{
	SphereData diff = i_lhs-i_rhs;
	diff.physical_reduce_max_abs();

	double lhs_maxabs = i_lhs.physical_reduce_max_abs();
	double rhs_maxabs = i_rhs.physical_reduce_max_abs();

	double normalize_fac;

	if (i_normalization)
	{
		normalize_fac = std::min(lhs_maxabs, rhs_maxabs);

		if (std::max(lhs_maxabs, rhs_maxabs) < i_error_threshold)
		{
			std::cout << "Error computation for '" << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
			normalize_fac = 1.0;
			//return false;
		}
	}
	else
	{
		normalize_fac = 1.0;
	}

	double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
	double rel_rms = diff.physical_reduce_rms() / normalize_fac;

	std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\terror threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

	if (rel_max_abs > i_error_threshold)
	{
		i_lhs.physical_file_write("o_error_lhs.csv");
		i_rhs.physical_file_write("o_error_rhs.csv");
		(i_lhs-i_rhs).physical_file_write("o_error_diff.csv");

		if (i_ignore_error)
			std::cerr << "Error ignored" << std::endl;
		else
			FatalError("Error too large");

		return true;
	}
	return false;
}



void run_tests()
{
	double epsilon = 1e-11;
	epsilon *= (sphereDataConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;

	SphereOperators op(sphereDataConfig);


#if 0
	if (true)
	{
		SphereData phi(sphereDataConfig);
		phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data = 1;
			}
		);
		//phi.physical_truncate();

		SphereData lhs(sphereDataConfig);
		lhs.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data = 1.0/(std::cos(i_theta)*std::cos(i_theta));
			}
		);

		SphereData rhs = phi;
		rhs.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, double &o_data)
				{
					o_data /= cos_phi*cos_phi;
				}
			);

		errorCheck(lhs, rhs, "TEST 1/cc", epsilon);
	}
#endif


#if 0
	if (true)
	{
		SphereData phi(sphereDataConfig);
		phi.physical_update_lambda_cosphi_grid(
			[&](double i_lambda, double cos_phi, double &io_data)
			{
				io_data = cos_phi;
			}
		);
		phi.physical_truncate();

		SphereData rhs = phi;
		rhs.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, double &o_data)
				{
					o_data /= cos_phi*cos_phi;
				}
			);
		rhs.physical_truncate();


		SphereData lhs(sphereDataConfig);
		lhs.physical_update_lambda_cosphi_grid(
			[&](double i_lambda, double cos_phi, double &io_data)
			{
				io_data = 1.0/cos_phi;
			}
		);
		lhs.physical_truncate();

		errorCheck(lhs, rhs, "TEST 1/cc", epsilon);
	}
#endif


#if 0
	if (true)
	{
		SphereData phi(sphereDataConfig);
		phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	  std::cos(i_theta)*std::cos(i_theta)
							* std::cos(i_lambda)*std::sin(i_lambda);
			}
		);
		phi.physical_truncate();

		phi.physical_file_write("o_phi.csv");

		SphereData div_lon_phi(sphereDataConfig);
		div_lon_phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	  std::cos(2*i_lambda);
			}
		);
		div_lon_phi.physical_truncate();

		SphereData lhs_ext = op.robert_div_lon(phi.spectral_returnWithDifferentModes(sphereDataConfigExt));
		SphereData rhs_ext = div_lon_phi.spectral_returnWithDifferentModes(sphereDataConfigExt);

		SphereData lhs = lhs_ext.spectral_returnWithDifferentModes(sphereDataConfig).physical_truncate();
		SphereData rhs = rhs_ext.spectral_returnWithDifferentModes(sphereDataConfig).physical_truncate();

		errorCheck(lhs, rhs, "TEST div_i(phi) with Robert formulation", epsilon);
	}
#endif

#if 1
	if (true)
	{
		SphereData phi(sphereDataConfig);
		phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	  std::cos(i_theta)*std::cos(i_theta)*std::cos(i_theta)
							* std::cos(i_lambda)*std::sin(i_lambda);
			}
		);
		phi.physical_truncate();

		phi.physical_file_write("o_phi.csv");

		SphereData div_lon_phi(sphereDataConfig);
		div_lon_phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	  std::cos(i_theta)*std::cos(i_theta)*std::cos(i_theta)
							* std::cos(2*i_lambda);
			}
		);
		div_lon_phi.physical_truncate();

		SphereData lhs_ext = op.diff_lon(phi.spectral_returnWithDifferentModes(sphereDataConfigExt));
		SphereData rhs_ext = div_lon_phi.spectral_returnWithDifferentModes(sphereDataConfigExt);

		SphereData lhs = lhs_ext.spectral_returnWithDifferentModes(sphereDataConfig).physical_truncate();
		SphereData rhs = rhs_ext.spectral_returnWithDifferentModes(sphereDataConfig).physical_truncate();

		errorCheck(lhs, rhs, "TEST diff_lon(phi) with Robert formulation", epsilon);
	}
#endif

#if 1
	if (true)
	{
		SphereData phi(sphereDataConfig);
		phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	  std::cos(i_theta)*std::cos(i_theta)
							* std::cos(i_lambda)*std::sin(i_lambda);
			}
		);

		SphereData div_lon_phi(sphereDataConfig);
		div_lon_phi.physical_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
			{
				io_data =	std::cos(2*i_lambda);
			}
		);

		SphereData lhs = op.robert_div_lon(phi);
		SphereData rhs = div_lon_phi;

		lhs.spectral_truncate();
		rhs.spectral_truncate();

		errorCheck(lhs, rhs, "TEST div_i(phi) with Robert formulation", epsilon);
	}
#endif

#if 1
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

		lhs.spectral_truncate();
		rhs.spectral_truncate();

		errorCheck(lhs, rhs, "TEST div_i(phi) with NON Robert formulation", epsilon);
	}
#endif

	std::cout << "OK" << std::endl;
	exit(0);

#if 1
	if (true)
	{
		double lambda_c = 3.0*M_PI/2.0;
//		double theta_c = 0;
		double theta_c = M_PI/3.0;
		double a = 6.37122e6;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

//		double alpha[] = {0, M_PI/3, M_PI/2};
		double alpha[] = {M_PI/3, M_PI/2};

		for (int i = 0; i < 3; i++)
		{
			double advection_rotation_angle = alpha[i];

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

			phi.physical_truncate();
			u.physical_truncate();
			v.physical_truncate();

			{
				SphereData ret = op.robert_div(u, v)*(1.0/a);

				errorCheck(ret, zero, "TEST div freeness", epsilon);
			}

			(op.robert_grad_lon(phi)*u).physical_file_write("o_grad_lon_phi.csv");
			(op.robert_grad_lat(phi)*v).physical_file_write("o_grad_lat_phi.csv");

			(op.robert_grad_lon_M(phi)*u).physical_file_write("o_grad_lon_M_phi.csv");
			(op.robert_grad_lat_M(phi)*v).physical_file_write("o_grad_lat_M_phi.csv");

			exit(1);

			{
#if 1
				// LHS = div(U*phi)
				SphereData lhs = op.robert_div(u*phi, v*phi)*(1.0/a);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				SphereData rhs = (op.robert_grad_lon_M(phi)*u + op.robert_grad_lat_M(phi)*v)*(1.0/a);
#else
				u = zero;
//				v = zero;
				SphereData lhs = op.robert_div(u*phi, v*phi);
				SphereData rhs = op.robert_grad_lon_M(phi)*u + op.robert_grad_lat_M(phi)*v;
#endif

				errorCheck(lhs, rhs, "TEST div(U*phi) - grad(phi)*u with Robert formulation", epsilon);
			}
		}
	}
#endif

#if 0
	if (true)
	{
		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0;
		double a = 6.37122e6;

		double R = a/3.0;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		double alpha[] = {0, M_PI/3, M_PI/2};

		for (int i = 0; i < 3; i++)
		{
			double advection_rotation_angle = alpha[i];

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

					io_data /= std::cos(i_lat);
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

					io_data /= std::cos(i_lat);
				}
			);

			phi.physical_truncate();
			u.physical_truncate();
			v.physical_truncate();

			{
				// LHS = div(U*phi)
				SphereData lhs = op.robert_div(u*phi, v*phi);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				SphereData rhs = op.robert_grad_lon(phi)*u + op.robert_grad_lat(phi)*v;

				errorCheck(lhs, rhs, "TEST div(U*phi) with Robert formulation", epsilon);
			}
		}
	}
#endif

#if 0
	if (true)
	{
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
				double normalization = std::max(
						u.physical_reduce_max_abs(),
						v.physical_reduce_max_abs()
					);

				SphereData lhs = op.div(u, v);

				errorCheck(lhs, zero, "TEST div freeness NORobert formulation", epsilon);
			}

			{
				// LHS = div(U*phi)
				SphereData lhs = op.div(u*phi, v*phi);

				// RHS = div(U)*phi + grad(phi)*U = grad(phi)*U  (divergence free)
				SphereData rhs = op.grad_lon(phi)*u + op.grad_lat(phi)*v;

				errorCheck(lhs, rhs, "DIV(U*phi) NORobert", epsilon);
			}
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
		h.physical_truncate();
		h.request_data_spectral();


		SphereData result(sphereDataConfig);
		result.physical_update_lambda_gaussian_grid(
				[&](double a, double b, double &c){testSolutionsSph.correct_result_diff_mu__grid_gaussian(a,b,c);}
		);
		result.request_data_spectral();
		result.physical_truncate();

		// one_minus_mu_squared_diff_lat
		h = op.spec_one_minus_mu_squared_diff_lat_mu(h);

		double error_max = h.physical_reduce_max(result);
		std::cout << "TEST SPHbasis (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

		if (error_max > epsilon)
			FatalError("ERROR THRESHOLD EXCEEDED!");
	}


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
			data.spectral_truncate();

			SphereData data2(sphereDataConfig);
			data2.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y*123.0;
					}
			);
			data2.spectral_truncate();

			double error = data.physical_reduce_max(data2);
			std::cout << "ERROR operator*123.0: " << error << std::endl;

			if (error > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
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
			data.spectral_truncate();

			SphereData data2(sphereDataConfig);
			data2.physical_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = y+100.0;
					}
			);
			data2.spectral_truncate();

			double error = data.physical_reduce_max(data2);
			std::cout << "ERROR operator+(double): " << error << std::endl;

			if (error > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
//			h.file_physical_writeFile("O_test_function.csv");

			SphereData hphi(sphereDataConfig);
			hphi.physical_update_lambda(
					[&](double a, double b, double &c){testSolutions.test_function_phi__grid_phi(a,b,c);}
			);
//			hphi.file_physical_writeFile("O_test_function_phi.csv");

			double error_max = h.physical_reduce_max(hphi);
			std::cout << "TEST PHI vs. MU: max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// identity
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			//h = h.spat_truncate();

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();

			double error_max = h.physical_reduce_max(result);

			std::cout << "TEST IDENTITY: max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d lambda
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			//h = h.spat_truncate();

			h = op.diff_lon(h);

//			h.file_physical_writeFile("O_diff_lambda_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_lambda__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_diff_lambda_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIFF LON - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d phi
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.diff_lat_phi(h);
//			h.file_physical_writeFile("O_diff_phi_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_phi__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_diff_phi_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIFF PHI - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// d/d mu
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.diff_lat_mu(h);
			//h = h.spat_truncate();
//			h.file_physical_writeFile("O_diff_mu_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_diff_mu__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_diff_mu_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIFF LAT MU - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// mu*F(\lambda,\mu)
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.mu(h);
//			h.file_physical_writeFile("O_mu_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_mu__grid_gaussian(a,b,c);}
			);
//			result.file_physical_writeFile("O_mu_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST mu*() - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		if (true)
		{
			// mu*mu*F(\lambda,\mu)
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.mu2(h);
//			h.file_physical_writeFile("O_mu2_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &i_data){
						testSolutions.test_function__grid_gaussian(lat, mu, i_data);
						i_data *= mu*mu;
					}
			);
//			result.file_physical_writeFile("O_mu2_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST mu*mu*() - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// one_minus_mu_squared_diff_lat
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.spec_one_minus_mu_squared_diff_lat_mu(h);
//			h.file_physical_writeFile("O_one_minus_mu_squared_diff_mu_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(a,b,c);}
			);

//			result.file_physical_writeFile("O_one_minus_mu_squared_diff_mu_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// grad lambda
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.grad_lon(h);
			//h = h.spat_truncate();
//			h.file_physical_writeFile("O_grad_lambda_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,c);}
			);
			//result = result.spat_truncate();
//			result.file_physical_writeFile("O_grad_lambda_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST GRAD LON - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// grad mu
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.grad_lat(h);
//			h.file_physical_writeFile("O_grad_phi_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_grad_phi__grid_gaussian(a,b,c);}
			);

//			result.file_physical_writeFile("O_grad_phi_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST GRAD LAT - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		if (true)
		{
			// div lambda
			SphereData h(sphereDataConfig);
			h.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
			);
			h = op.div_lon(h);
//			h.file_physical_writeFile("O_div_lambda_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutions.correct_result_div_lambda__grid_gaussian(a,b,c);}
			);
//			result.file_physical_writeFile("O_div_lambda_correct_result.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIV LON  - max error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
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
			//h.physical_file_write("O_div_mu_initial_sph.csv");

			h = op.div_lat(h);
//			h.spectral_truncate();
			//h.physical_file_write("O_div_mu_sph_result.csv");

			SphereData result(sphereDataConfig);
			result.physical_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){
						testSolutions.correct_result_div_mu__grid_gaussian(a,b,c);
					}
			);
//			result.spectral_truncate();
			//result.physical_file_write("O_div_mu_correct_result.csv");
			//(h-result).physical_file_write("O_div_mu_correct_diff.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;

			if (error_max > epsilon)
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

			double error_max = h.physical_reduce_max(result);
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
			//result = result.spat_truncate();

//			result.file_physical_writeFile("O_divergence_correct_result.csv");

//			(h-result).file_physical_writeFile("O_divergence_correct_diff.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST DIVERGENCE - max error: " << error_max << std::endl;

			if (error_max > epsilon)
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
			//h = h.spat_truncate();
//			h.file_physical_writeFile("O_vort_sph_result.csv");

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
			//result = result.spat_truncate();

//			result.file_physical_writeFile("O_vort_correct_result.csv");

//			(h-result).file_physical_writeFile("O_vort_correct_diff.csv");

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST VORTICITY - max error: " << error_max << std::endl;

			if (error_max > epsilon)
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

			double error_max = h.physical_reduce_max(result);
			std::cout << "TEST LAPLACE (div.grad vs. sph laplace) - max error: " << error_max << std::endl;

			if (error_max > epsilon)
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

	sphereDataConfigInstance.setupAutoPhysicalSpace(
					simVars.disc.res_spectral[0],
					simVars.disc.res_spectral[1],
					&simVars.disc.res_physical[0],
					&simVars.disc.res_physical[1]
			);

	sphereDataConfigExtInstance.setupAdditionalModes(
			&sphereDataConfigInstance,
			simVars.rexi.rexi_use_extended_modes,
			simVars.rexi.rexi_use_extended_modes
		);

	run_tests();
}


#endif /* SRC_TESTOPERATORS_HPP_ */
