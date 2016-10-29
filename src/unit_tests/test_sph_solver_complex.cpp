/*
 * AppTestSPHSolverComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_TESTSPHSOLVERS_COMPLEX_HPP_
#define SRC_TESTSPHSOLVERS_COMPLEX_HPP_


#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphBandedMatrixComplex.hpp>
#include <rexi/REXI.hpp>
#include <libmath/BandedMatrixComplex.hpp>



SimulationVariables simVars;
SphereOperatorsComplex opComplex;



/**
 * Run with
 *
 * 	$ ./build/sh_example T32 P2
 */
void run_tests(
		SphereDataConfig *sphConfig
)
{
	double epsilon = 1e-12;
	epsilon *= (sphConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;

	std::cout << std::setprecision(10);

	SphereTestSolutions_Gaussian testSolutions;

	if (sphConfig->spectral_modes_n_max < 32)
	{
		std::cerr << "WARNING: AT LEAST 32 MODES REQUIRED for proper accuracy!!!" << std::endl;
	}

	/**
	 * Use test function as expected result
	 */
	SphereDataComplex x_result(sphConfig);
	x_result.physical_update_lambda_gaussian_grid(
			[&](double lat, double mu, std::complex<double> &io_data){
				double tmp;
				testSolutions.test_function__grid_gaussian(lat,mu,tmp);
				io_data = tmp;
			}
	);

	double r = simVars.sim.earth_radius;
	double two_omega = 2.0*simVars.sim.coriolis_omega;


	REXI rexi;
	rexi.setup(0.2, 64);


	for (std::size_t i = 0; i < rexi.alpha.size(); i++)
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
			std::cout << "Test Zx = c*Phi(mu)";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			sphSolver.solver_component_scalar_phi(alpha);

			/*
			 * Setup RHS = scalar_a * phi(lambda,mu)
			 */
			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = tmp*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		/*
		 * Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)
		 */
		if (true)
		{
			std::cout << "Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			sphSolver.solver_component_scalar_phi(alpha);
			sphSolver.solver_component_mu_phi();
//				sphSolver.lhs.print();

			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat,mu,fun);

						io_data = mu*fun+alpha*fun;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);

//				x_numerical.spat_write_file("O_numerical.csv");
//				x_result.spat_write_file("O_result.csv");
//				(x_result-x_numerical).spat_write_file("O_diff.csv");

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		/*
		 * Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)
		 */
		if (true)
		{
			std::cout << "Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> a(alpha);
			sphSolver.solver_component_scalar_phi(a);
			sphSolver.solver_component_one_minus_mu_mu_diff_mu_phi();

			SphereDataComplex b(sphConfig);
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

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
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
			std::cout << "Test Z1 = alpha^4*Phi(mu)";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = (alpha*alpha)*(alpha*alpha);
			sphSolver.solver_component_rexi_z1(scalar, r);

			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = tmp*scalar;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		/*
		 * Test Z2 = mu^2 * Phi(mu)
		 */
		if (true)
		{
			std::cout << "Test Z2 = mu^2*Phi(lam,mu)";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = alpha*alpha*two_omega*two_omega;
			//sphSolver.solver_component_scalar_phi(scalar_a);
			sphSolver.solver_component_rexi_z2(scalar, r);

			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data){
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = scalar*(mu*mu)*tmp;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		/*
		 * Test Z3 = mu^4 * Phi(mu)
		 */
		if (true)
		{
			std::cout << "Test Z3 = mu^4*Phi(lam,mu)";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 4);

			std::complex<double> scalar = two_omega*two_omega*two_omega*two_omega;
			//sphSolver.solver_component_scalar_phi(scalar_a);
			sphSolver.solver_component_rexi_z3(scalar, r);

			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = scalar*(mu*mu)*(mu*mu)*tmp;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		/*
		 * Test Z4 = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)
		 */
		if (true)
		{
			std::cout << "Test Z4 = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = -alpha*alpha*two_omega;
			sphSolver.solver_component_rexi_z4(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphConfig);
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

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		/*
		 * Test Z5 = grad_j(mu) mu^2 grad_i(Phi(lam,mu))
		 */
		if (true)
		{
			std::cout << "Test Z5 = grad_j(mu) mu^2 grad_i(Phi(lam,mu))";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = two_omega*two_omega*two_omega;
			sphSolver.solver_component_rexi_z5(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dgrad_lambda_fun;
						testSolutions.correct_result_grad_lambda__grid_gaussian(lat, mu, dgrad_lambda_fun);

						io_data = scalar/(r*r)*std::sqrt(1.0-mu*mu)*mu*mu*dgrad_lambda_fun + fun*alpha;
					}
			);

			SphereDataComplex x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		/*
		 * Test Z6 = grad_j(mu) mu grad_j(Phi(lam,mu))
		 */
		if (true)
		{
			std::cout << "Test Z6 = grad_j(mu) mu grad_j(Phi(lam,mu))";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = 2.0*alpha*two_omega*two_omega;
			sphSolver.solver_component_rexi_z6(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphConfig);
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

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		/*
		 * Test Z7 = laplace(Phi(lam,mu))
		 */
		if (true)
		{
			std::cout << "Test Z7 = laplace(Phi(lam,mu))";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = 1.0;
			sphSolver.solver_component_rexi_z7(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat, mu, tmp);
						io_data = tmp;
					}
			);

			b = opComplex.laplace(b)*scalar/(r*r) + b*alpha;

			SphereDataComplex x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}



		/*
		 * Test Z8 = mu*mu*laplace(Phi(lam,mu))
		 */
		if (true)
		{
			std::cout << "Test Z8 = mu*mu*laplace(Phi(lam,mu))";

			SphBandedMatrixComplex<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = 1.0;
			sphSolver.solver_component_rexi_z8(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereDataComplex b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat, mu, tmp);
						io_data = tmp;
					}
			);

			SphereDataComplex c = opComplex.mu2(opComplex.laplace(b))*scalar/(r*r) + b*alpha;

			SphereDataComplex x_numerical = sphSolver.solve(c);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}


		/*
		 * Test direct Helmholtz problem solver
		 * (a + b D^2) x = rhs
		 */
		if (true)
		{
			std::cout << "Test Zx = a + b*laplace";

			SphereDataComplex testb(sphConfig);
			testb.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat, mu, tmp);
						io_data = tmp;
					}
			);

			std::complex<double> a = alpha;
			std::complex<double> b = alpha*2.0-1.0;
			double r = simVars.sim.earth_radius;

			SphereDataComplex testc = a*testb + (b/(r*r))*opComplex.laplace(testb);

			SphereDataComplex x_numerical = testc.spectral_solve_helmholtz(a, b, r);

			double error_max = x_numerical.physical_reduce_error_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
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



#endif /* SRC_TESTSPHSOLVERS_HPP_ */
