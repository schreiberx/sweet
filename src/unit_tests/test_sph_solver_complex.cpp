/*
 * AppTestSPHSolverComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_TESTSPHSOLVERS_COMPLEX_HPP_
#define SRC_TESTSPHSOLVERS_COMPLEX_HPP_


#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <libmath/BandedMatrixPhysicalComplex.hpp>
#include <rexi/REXI_Terry.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalComplex.hpp>
#include <sweet/sphere/ErrorCheck.hpp>

#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>



SphereDataConfig *sphereDataConfig;

SimulationVariables simVars;


/**
 * Run with
 */
void run_tests()
{
	double max_scalar = std::max(std::max(simVars.sim.gravitation, simVars.sim.h0), simVars.sim.f0);

	double epsilon = 1e-11*max_scalar;
	epsilon *= (sphereDataConfig->spectral_modes_n_max);

	std::cout << "Using max allowed error of " << epsilon << std::endl;

	std::cout << std::setprecision(10);

	SphereTestSolutions_Gaussian testSolutions;

	if (sphereDataConfig->spectral_modes_n_max < 32)
	{
		std::cerr << "WARNING: AT LEAST 32 MODES REQUIRED for proper accuracy!!!" << std::endl;
	}

	/**
	 * Use test function as expected result
	 */
	SphereDataComplex x_result(sphereDataConfig);
	x_result.physical_update_lambda_gaussian_grid(
			[&](double lat, double mu, std::complex<double> &io_data){
				double tmp;
				testSolutions.test_function__grid_gaussian(lat,mu,tmp);
				io_data = tmp;
			}
	);

	double r = simVars.sim.earth_radius;
	double two_omega = 2.0*simVars.sim.coriolis_omega;

	REXI_Terry<__float128> rexi("phi0", 0.2, 64);

	SphereOperatorsComplex opComplex(sphereDataConfig, 1);


	for (std::size_t i = 0; i < rexi.alpha.size(); i+=10)
	{
		std::complex<double> &alpha = rexi.alpha[i];
		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "Testing for alpha " << alpha << std::endl;
		std::cout << "***************************************************" << std::endl;


		/*
		 * Test Zx = c*Phi(mu)
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			sphSolver.solver_component_scalar_phi(alpha);

			/*
			 * Setup RHS = scalar_a * phi(lambda,mu)
			 */
			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = tmp*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Zx = c*Phi(mu)", epsilon);
		}

		/*
		 * Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			sphSolver.solver_component_scalar_phi(alpha);
			sphSolver.solver_component_mu_phi();

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat,mu,fun);

						io_data = mu*fun+alpha*fun;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)", epsilon);
		}

		/*
		 * Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)
		 */
		if (true)
		{

			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> a(alpha);
			sphSolver.solver_component_scalar_phi(a);
			sphSolver.solver_component_one_minus_mu_mu_diff_mu_phi();

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double one_minus_m2_dfun;
						testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(lat, mu, one_minus_m2_dfun);

						io_data = one_minus_m2_dfun + alpha*fun;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)", epsilon);
		}

		std::cout << "************************************************************" << std::endl;
		std::cout << "*** BASIC TESTS FINISHED ***********************************" << std::endl;
		std::cout << "************************************************************" << std::endl;
		std::cout << "*** TESTING REXI COMPONENTS ********************************" << std::endl;
		std::cout << "************************************************************" << std::endl;

		/*
		 * Test F1 = alpha^4 * Phi(mu)
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = (alpha*alpha)*(alpha*alpha);
			sphSolver.solver_component_rexi_z1(scalar, r);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = tmp*scalar;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z1 = alpha^4*Phi(mu)", epsilon);
		}

		/*
		 * Test Z2 = mu^2 * Phi(mu)
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = alpha*alpha*two_omega*two_omega;
			//sphSolver.solver_component_scalar_phi(scalar_a);
			sphSolver.solver_component_rexi_z2(scalar, r);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data){
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = scalar*(mu*mu)*tmp;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z2 = mu^2*Phi(lam,mu)", epsilon);
		}


		/*
		 * Test Z3 = mu^4 * Phi(mu)
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 4);

			std::complex<double> scalar = two_omega*two_omega*two_omega*two_omega;
			//sphSolver.solver_component_scalar_phi(scalar_a);
			sphSolver.solver_component_rexi_z3(scalar, r);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = scalar*(mu*mu)*(mu*mu)*tmp;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z3 = mu^4*Phi(lam,mu)", epsilon*1e+1);
		}


		/*
		 * Test Z4 = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = -alpha*alpha*two_omega;
			sphSolver.solver_component_rexi_z4(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dlambda_fun;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lat, mu, dlambda_fun);

						io_data = scalar/(r*r)*dlambda_fun + fun*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z4 = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)", epsilon);
		}
		/*
		 * Test Z4robert = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = -alpha*alpha*two_omega;
			sphSolver.solver_component_rexi_z4robert(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dlambda_fun;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lat, mu, dlambda_fun);

						io_data = scalar/(r*r)*dlambda_fun + fun*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z4robert = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)", epsilon);
		}


		/*
		 * Test Z5 = grad_j(mu) mu^2 grad_i(Phi(lam,mu))
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = two_omega*two_omega*two_omega;
			sphSolver.solver_component_rexi_z5(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dgrad_lambda_fun;
						testSolutions.correct_result_grad_lambda__grid_gaussian(lat, mu, dgrad_lambda_fun);

						io_data = (scalar/(r*r))*std::sqrt(1.0-mu*mu)*mu*mu*dgrad_lambda_fun + fun*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			std::cout << "WARNING: This component is very sensitive to small alpha values vs. large two_omega^3 !" << std::endl;
			ErrorCheck::check(x_numerical, x_result, "Test Z5 = grad_j(mu) mu^2 grad_i(Phi(lam,mu))", epsilon, two_omega>=100000);
		}


		/*
		 * Test Z5robert = grad_j(mu) mu^2 grad_i(Phi(lam,mu))
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = two_omega*two_omega*two_omega;
			scalar = two_omega*two_omega;//*two_omega;
			sphSolver.solver_component_rexi_z5robert(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dgrad_lambda_fun;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lat, mu, dgrad_lambda_fun);

						io_data = (scalar/(r*r))*mu*mu*dgrad_lambda_fun + fun*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			std::cout << "WARNING: This component is very sensitive to small alpha values vs. large two_omega^3 !" << std::endl;
			ErrorCheck::check(x_numerical, x_result, "Test Z5robert = grad_j(mu) mu^2 grad_i(Phi(lam,mu))", epsilon);
		}


		/*
		 * Test Z6 = grad_j(mu) mu grad_j(Phi(lam,mu))
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = 2.0*alpha*two_omega*two_omega;
			sphSolver.solver_component_rexi_z6(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dgrad_mu_fun;
						testSolutions.correct_result_grad_phi__grid_gaussian(lat, mu, dgrad_mu_fun);

						io_data = scalar/(r*r)*std::sqrt(1.0-mu*mu)*mu*dgrad_mu_fun + fun*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z6 = grad_j(mu) mu grad_j(Phi(lam,mu))", epsilon);
		}


		/*
		 * Test Z6robert = grad_j(mu) mu grad_j(Phi(lam,mu))
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = 2.0*alpha*two_omega*two_omega;
			sphSolver.solver_component_rexi_z6robert(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dgrad_mu_fun;
						testSolutions.correct_result_diff_mu__grid_gaussian(lat, mu, dgrad_mu_fun);

						io_data = scalar/(r*r)*mu*(1-mu*mu)*dgrad_mu_fun + fun*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z6robert = grad_j(mu) mu grad_j(Phi(lam,mu))", epsilon);
		}


		/*
		 * Test Z7 = laplace(Phi(lam,mu))
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = 1.0;
			sphSolver.solver_component_rexi_z7(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphereDataConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat, mu, tmp);
						io_data = tmp;
					}
			);

			b = opComplex.laplace(b,1)*scalar/(r*r) + b*alpha;

			SphereDataComplex x_numerical = sphSolver.solve(b);
			ErrorCheck::check(x_numerical, x_result, "Test Z7 = laplace(Phi(lam,mu))", epsilon);
		}


		/*
		 * Test Z8 = mu*mu*laplace(Phi(lam,mu))
		 */
		if (true)
		{
			SphBandedMatrixPhysicalComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphereDataConfig, 2);

			std::complex<double> scalar = 1.0;
			sphSolver.solver_component_rexi_z8(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b = opComplex.mu2(opComplex.laplace(x_result,1))*scalar/(r*r) + x_result*alpha;
			SphereDataComplex x_numerical = sphSolver.solve(b);

			ErrorCheck::check(x_numerical, x_result, "Test Z8 = mu*mu*laplace(Phi(lam,mu))", epsilon);
		}



		/*
		 * Test direct Helmholtz problem solver
		 * (a + b D^2) x = rhs
		 */
		if (true)
		{
			std::complex<double> a = alpha;
			std::complex<double> b = alpha*2.0-1.0;
			double r = simVars.sim.earth_radius;

			SphereDataComplex testc = a*x_result + (b/(r*r))*opComplex.laplace(x_result,1);

			SphereDataComplex x_numerical = testc.spectral_solve_helmholtz(a, b, r);

			ErrorCheck::check(x_numerical, x_result, "Test Zx = a + b*laplace", epsilon);
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

	SphereDataConfig p_sphereDataConfig;
	p_sphereDataConfig.setupAutoPhysicalSpace(
					simVars.disc.res_spectral[0],
					simVars.disc.res_spectral[1],
					&simVars.disc.res_physical[0],
					&simVars.disc.res_physical[1],
					simVars.misc.shtns_use_plans
			);

	sphereDataConfig = &p_sphereDataConfig;
	run_tests();
}



#endif /* SRC_TESTSPHSOLVERS_HPP_ */
