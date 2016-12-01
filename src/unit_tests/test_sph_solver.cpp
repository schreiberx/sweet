/*
 * AppTestSPHSolver.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_TESTSPHSOLVERS_HPP_
#define SRC_TESTSPHSOLVERS_HPP_


#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <libmath/BandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SWESphBandedMatrixPhysicalReal.hpp>


SimulationVariables simVars;
SphereOperators op;

void run_tests(
		SphereDataConfig *sphConfig
)
{
	double epsilon = 1e-12;
	epsilon *= (sphConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;

	std::cout << std::setprecision(10);

	{
		SphereTestSolutions_Gaussian testSolutions;

		if (sphConfig->spectral_modes_n_max < 32)
		{
			std::cerr << "WARNING: AT LEAST 32 MODES REQUIRED for proper accuracy!!!" << std::endl;
		}

		/**
		 * Use test function as expected result
		 */
		SphereData x_result(sphConfig);
		x_result.physical_update_lambda_gaussian_grid(
				[&](double lat, double mu, double &io_data){
					testSolutions.test_function__grid_gaussian(lat,mu,io_data);
				}
		);

		std::complex<double> alpha(3.0, 0.0);
		double r = simVars.sim.earth_radius;
		double two_omega = 2.0*simVars.sim.coriolis_omega;



		/*
		 * Test Zx = c*Phi(mu)
		 */
		if (true)
		{
			std::cout << "Test Zx = c*Phi(mu)";

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			sphSolver.solver_component_scalar_phi(alpha);

			/*
			 * Setup RHS = scalar_a * phi(lambda,mu)
			 */
			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						double tmp;
						testSolutions.test_function__grid_gaussian(lat,mu,tmp);
						io_data = tmp*alpha.real();
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			sphSolver.solver_component_scalar_phi(alpha);
			sphSolver.solver_component_mu_phi();
//				sphSolver.lhs.print();

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat,mu,fun);

						io_data = mu*fun+alpha.real()*fun;
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

//				x_numerical.spat_write_file("O_numerical.csv");
//				x_result.spat_write_file("O_result.csv");
//				(x_result-x_numerical).spat_write_file("O_diff.csv");

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> a(alpha.real());
			sphSolver.solver_component_scalar_phi(a);
			sphSolver.solver_component_one_minus_mu_mu_diff_mu_phi();

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double one_minus_m2_dfun;
						testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(lat, mu, one_minus_m2_dfun);

						io_data = one_minus_m2_dfun + alpha.real()*fun;
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = (alpha*alpha)*(alpha*alpha);
			sphSolver.solver_component_rexi_z1(scalar, r);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						testSolutions.test_function__grid_gaussian(lat,mu,io_data);
						io_data *= scalar.real();
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = alpha*alpha*two_omega*two_omega;
			//sphSolver.solver_component_scalar_phi(scalar_a);
			sphSolver.solver_component_rexi_z2(scalar, r);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data){
						testSolutions.test_function__grid_gaussian(lat,mu,io_data);
						io_data = scalar.real()*(mu*mu)*io_data;
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

		/*
		 * Test Z2 = mu^2 * Phi(mu)
		 */
#if 0
		if (true)
		{
			std::cout << "Test Z2 = mu^2*Phi(lam,mu) (DOUBLE MODES)";

			SphereDataConfig sphConfigDouble;
			sphConfigDouble.setupDouble(sphConfig);

			SphBandedMatrixPhysicalReal< std::complex<double> > sphSolver;
			sphSolver.setup(&sphConfigDouble, 2);

			std::complex<double> scalar = alpha*alpha*two_omega*two_omega;
			sphSolver.solver_component_rexi_z2(scalar, r);

			SphereData b(&sphConfigDouble);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data){
						testSolutions.test_function__grid_gaussian(lat,mu,io_data);
						io_data = scalar.real()*(mu*mu)*io_data;
					}
			);

			SphereData x_numerical_double = sphSolver.solve(b);

			SphereData x_numerical(sphConfig);

			x_numerical.spec_set_zero();
			x_numerical.spec_update_lambda(
					[&](int n, int m, std::complex<double> &io_data)
					{
						io_data = x_numerical_double.spec_get(n, m);
					}
			);

			double error_max = x_numerical.physical_reduce_max(x_result);
			std::cout << " ||| Error: " << error_max << std::endl;

			if (error_max > epsilon)
				FatalError("ERROR THRESHOLD EXCEEDED!");
		}

#endif

		/*
		 * Test Z3 = mu^4 * Phi(mu)
		 */
		if (true)
		{
			std::cout << "Test Z3 = mu^4*Phi(lam,mu)";

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 4);

			std::complex<double> scalar = two_omega*two_omega*two_omega*two_omega;
			//sphSolver.solver_component_scalar_phi(scalar_a);
			sphSolver.solver_component_rexi_z3(scalar, r);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data){
						testSolutions.test_function__grid_gaussian(lat,mu,io_data);
						io_data = scalar.real()*(mu*mu)*(mu*mu)*io_data;
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = -alpha*alpha*two_omega;
			sphSolver.solver_component_rexi_z4(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dlambda_fun;
						testSolutions.correct_result_diff_lambda__grid_gaussian(lat, mu, dlambda_fun);

						io_data = scalar.real()/(r*r)*dlambda_fun + fun*alpha.real();
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = two_omega*two_omega*two_omega;
			sphSolver.solver_component_rexi_z5(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dgrad_lambda_fun;
						testSolutions.correct_result_grad_lambda__grid_gaussian(lat, mu, dgrad_lambda_fun);

						io_data = scalar.real()/(r*r)*std::sqrt(1.0-mu*mu)*mu*mu*dgrad_lambda_fun + fun*alpha.real();
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = 2.0*alpha*two_omega*two_omega;
			sphSolver.solver_component_rexi_z6(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						double fun;
						testSolutions.test_function__grid_gaussian(lat, mu, fun);

						double dgrad_mu_fun;
						testSolutions.correct_result_grad_phi__grid_gaussian(lat, mu, dgrad_mu_fun);

						io_data = scalar.real()/(r*r)*std::sqrt(1.0-mu*mu)*mu*dgrad_mu_fun + fun*alpha.real();
					}
			);

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = 1.0;
			sphSolver.solver_component_rexi_z7(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						testSolutions.test_function__grid_gaussian(lat, mu, io_data);
					}
			);

			b = op.laplace(b)*scalar.real()/(r*r) + b*alpha.real();

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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

			SphBandedMatrixPhysicalReal<std::complex<double>> sphSolver;
			sphSolver.setup(sphConfig, 2);

			std::complex<double> scalar = 1.0;
			sphSolver.solver_component_rexi_z8(scalar, r);

			// ADD OFFSET FOR NON-SINGULAR SOLUTION
			sphSolver.solver_component_scalar_phi(alpha);

			SphereData b(sphConfig);
			b.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data)
					{
						testSolutions.test_function__grid_gaussian(lat, mu, io_data);
					}
			);

			b = op.mu2(op.laplace(b))*scalar.real()/(r*r) + b*alpha.real();

			SphereData x_numerical = sphSolver.solve(b);

			double error_max = x_numerical.physical_reduce_max(x_result);
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
