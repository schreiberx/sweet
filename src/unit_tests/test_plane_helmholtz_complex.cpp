
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <sweet/Stopwatch.hpp>
#include "../include/sweet/plane/PlaneData.hpp"
#include "../include/sweet/plane/PlaneDataComplex.hpp"
#include <sweet/SimulationVariables.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>
#include <complex>

#include "../programs/swe_plane/SWE_Plane_REXI.hpp"
#include "../programs/rexiswe/RexiSWE_HelmholtzSolver.hpp"

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

SimulationVariables simVars;

#if SWEET_DEBUG
#include <fenv.h>
static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked.  */

  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif

typedef std::complex<double> complex;


int main(int i_argc, char *i_argv[])
{
#if SWEET_DEBUG
	trapfpe();
#endif

	SimulationVariables simVars;
	simVars.disc.use_spectral_basis_diffs = 1;

	const char *bogus_var_names[] =
	{
			"use-specdiff-for-complex-array",	/// use finite differences for complex array
			"helmholtz-solver-id",		/// Which Helmholtz solver to use
			"helmholtz-smoother-sor",		/// Helmholtz solver overrelaxation
			"helmholtz-smoother-eps",		/// Helmholtz solver error threshold
			nullptr
	};

	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "      --use-specdiff-for-complex-array=[0/1]	Use finite-differences for derivatives in spectral space" << std::endl;
		std::cout << "      --helmholtz-solver-id=[0/1/2]			Which Helmholtz solver should we use" << std::endl;
		std::cout << "                          0: Spectral solver (default)" << std::endl;
		std::cout << "                          1: Iterative solver (Jacobi)" << std::endl;
		std::cout << "                          2: Iterative solver (CG)" << std::endl;
		std::cout << "                          3: MG solver (Jacobi smoother)" << std::endl;
		std::cout << "                          4: MG solver (CG smoother)" << std::endl;
		std::cout << "      --helmholtz-smoother-sor=[float]			Over relaxation parameter, default:1 (no overrelaxation)" << std::endl;
		std::cout << "      --helmholtz-smoother-eps=[float]			Error threshold" << std::endl;
		std::cout << std::endl;
		return -1;
	}


	/*
	 * use finite differences for differential operators in complex array
	 */
	bool use_spectral_differences_for_complex_array = 0;
	int helmholtz_solver_id = 0;
	double helmholtz_solver_sor = 1;
	double helmholtz_solver_eps = 1e-7;

	if (simVars.bogus.var[0] != "")
		use_spectral_differences_for_complex_array = atof(simVars.bogus.var[0].c_str());

	if (simVars.bogus.var[1] != "")
		helmholtz_solver_id = atof(simVars.bogus.var[1].c_str());

	if (simVars.bogus.var[2] != "")
		helmholtz_solver_sor = atof(simVars.bogus.var[2].c_str());

	if (simVars.bogus.var[3] != "")
		helmholtz_solver_eps = atof(simVars.bogus.var[3].c_str());


	if (!use_spectral_differences_for_complex_array)
	{
		std::cout << "********************************************************" << std::endl;
		std::cout << "*** Using finite-differences for complex array" << std::endl;
		std::cout << "********************************************************" << std::endl;
	}

	if (simVars.disc.use_spectral_basis_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	double freq_x = 2.0;
	double freq_y = 2.0;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.res_physical[0];
	std::size_t res_y = simVars.disc.res_physical[1];

	std::size_t max_res = 2048;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	for (; res_x <= max_res && res_y <= max_res; res_x *= 2, res_y *= 2)
	{
		/*
		 * error tolerance for machine accuracy
		 *
		 * We assume 1e-12 for double precision
		 */
//		double eps = 1e-9*tolerance_increase;
//		double eps_conv = 1e-3*tolerance_increase;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};


		simVars.disc.res_physical[0] = res[0];
		simVars.disc.res_physical[1] = res[1];
		simVars.reset();

		planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical, simVars.misc.reuse_spectral_transformation_plans);

		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		PlaneDataComplex h_cart(res);



		/**
		 * Test iterative solver for Helmholtz problem
		 */
		{

			PlaneDataComplex h_diff2_x(res);
			PlaneDataComplex h_diff2_y(res);

//			Operators2D op(parameters.discretization.res, parameters.sim.domain_size, parameters.disc.use_spectral_diffs);

			for (int j = 0; j < simVars.disc.res_physical[1]; j++)
			{
				for (int i = 0; i < simVars.disc.res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.res_physical[1];

					// H to reconstruct
					h_cart.p_physical_set(
						j, i,
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y),
						sin(2.0*freq_y*M_PIl*y)*cos(2.0*freq_x*M_PIl*x)
					);

					h_diff2_x.p_physical_set(
						j, i,
						freq_x*freq_x*M_PIl*M_PIl*(-1.0)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(simVars.sim.domain_size[0]*simVars.sim.domain_size[0]),
						-sin(2.0*freq_y*M_PIl*y)*2.0*freq_x*M_PIl*2.0*freq_x*M_PIl*cos(2.0*freq_x*M_PIl*x)/(simVars.sim.domain_size[1]*simVars.sim.domain_size[1])
					);

					h_diff2_y.p_physical_set(
						j, i,
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*freq_y*M_PIl*cos(freq_y*M_PIl*y)/(simVars.sim.domain_size[1]*simVars.sim.domain_size[1]),
						2.0*freq_y*2.0*freq_y*M_PIl*M_PIl*(-1.0)*sin(2.0*freq_y*M_PIl*y)*cos(2.0*freq_x*M_PIl*x)/(simVars.sim.domain_size[0]*simVars.sim.domain_size[0])
					);
				}
			}



			SWE_Plane_REXI rexiSWE;
			rexiSWE.setup(
					1,			// time step size
					0.2,		// h
					16,			// REXI M
					0,			// REXI N (0 for auto detection)
					simVars.sim.f0,
					res,
					simVars.sim.domain_size,
					true,		// use only half of REXI
					use_spectral_differences_for_complex_array,	// use finite differences
					helmholtz_solver_id,						// iterative solver
					helmholtz_solver_eps
				);


			PlaneDataComplex op_diff2_c_x(res);
			PlaneDataComplex op_diff2_c_y(res);
			op_diff2_c_x.op_setup_diff2_x(simVars.sim.domain_size, use_spectral_differences_for_complex_array);
			op_diff2_c_y.op_setup_diff2_y(simVars.sim.domain_size, use_spectral_differences_for_complex_array);


			double inv_helm_h[2];
			inv_helm_h[0] = (double)res[0]/(double)simVars.sim.domain_size[0];
			inv_helm_h[1] = (double)res[1]/(double)simVars.sim.domain_size[1];

			double scalar_Dx = inv_helm_h[0]*inv_helm_h[0];
			double scalar_Dy = inv_helm_h[1]*inv_helm_h[1];
			double scalar_C = -(2.0*(inv_helm_h[0]*inv_helm_h[0]) + 2.0*(inv_helm_h[1]*inv_helm_h[1]));

			PlaneDataComplex h(res);
			h.setAll(0, 0);

			double tau = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 1);

			PlaneOperators op(simVars.disc.res_physical, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

			for (std::size_t i = 0; i < rexiSWE.rexi.alpha.size(); i++)
			{
				std::complex<double> alpha = -rexiSWE.rexi.alpha[i]/tau;
				std::complex<double> kappa = alpha*alpha + simVars.sim.f0*simVars.sim.f0;

				double gh = simVars.sim.gravitation * simVars.setup.h0;

				std::cout << "TAU: " << tau << std::endl;
				std::cout << "KAPPA: " << kappa << std::endl;
				std::cout << "gh: " << gh << std::endl;

				// compute RHS
				PlaneDataComplex rhs(res);

				// analytical rhs
				//rhs = kappa*h_cart - gh*(h_diff2_x + h_diff2_y);

				// numerical rhs
//				rhs  = kappa*h_cart - gh*(op_diff2_c_x(h_cart.toSpec()).toCart() + op_diff2_c_y(h_cart.toSpec()).toCart());
				rhs  = kappa*h_cart - gh*(h_cart.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C));

				rhs = rhs*(1.0/tau);

//				h.setAll(0, 0);

				int iter_max = 200000;

				Stopwatch watch;
				watch.start();

					bool retval = true;
					if (helmholtz_solver_id == 0)
					{
						if (simVars.misc.verbosity > 1)
							std::cout << "REXI: Using helmholtz_spectral_solver" << std::endl;

						rexiSWE.helmholtz_spectral_solver_cart(	// DIRECT SPECTRAL SOLVER
								kappa,
								simVars.sim.gravitation * simVars.setup.h0,
								rhs,
								h
							);
					}
					else if (helmholtz_solver_id == 1)
					{
						if (simVars.misc.verbosity > 1)
							std::cout << "REXI: Using helmholtz_iterative_smoother_jacobi" << std::endl;

						retval = RexiSWE_HelmholtzSolver::smoother_jacobi(	// ITERATIVE JACOBI SOLVER
								kappa,
								simVars.sim.gravitation * simVars.setup.h0,
								rhs,
								h,
								simVars.sim.domain_size,
								helmholtz_solver_eps,
								iter_max,
								helmholtz_solver_sor,	// SOR omega
								simVars.misc.verbosity
							);
					}
					else if (helmholtz_solver_id == 2)
					{
						if (simVars.misc.verbosity > 1)
							std::cout << "REXI: Using helmholtz_iterative_smoother_conjugate_gradient" << std::endl;

						retval = RexiSWE_HelmholtzSolver::smoother_conjugate_gradient(	// CONJUGATE GRADIENT ITERATIVE SOLVER
							kappa,
							simVars.sim.gravitation * simVars.setup.h0,
							rhs,
							h,
							simVars.sim.domain_size,
							helmholtz_solver_eps,
							iter_max,
							-123,
							simVars.misc.verbosity
						);
					}
					else if (helmholtz_solver_id == 3)
					{
						if (simVars.misc.verbosity > 1)
							std::cout << "REXI: Using helmholtz_iterative_solve_mg(RexiSWE::helmholtz_iterative_smoother_jacobi)" << std::endl;

						retval = RexiSWE_HelmholtzSolver::multigrid(	// MG with jacobi
							kappa,
							simVars.sim.gravitation * simVars.setup.h0,
							rhs,
							h,
							RexiSWE_HelmholtzSolver::smoother_jacobi,
							simVars.sim.domain_size,
							helmholtz_solver_eps,
							iter_max,
							helmholtz_solver_sor,
							2,
							simVars.misc.verbosity
						);
					}
					else if (helmholtz_solver_id == 4)
					{
						if (simVars.misc.verbosity > 1)
							std::cout << "REXI: Using helmholtz_iterative_solve_mg(RexiSWE::helmholtz_iterative_smoother_conjugate_gradient)" << std::endl;

						retval = RexiSWE_HelmholtzSolver::multigrid(	// MG with CG
							kappa,
							simVars.sim.gravitation * simVars.setup.h0,
							rhs,
							h,
							RexiSWE_HelmholtzSolver::smoother_conjugate_gradient,
							simVars.sim.domain_size,
							helmholtz_solver_eps,
							iter_max,
							-999,	/// sor obsolete for CG
							-1,
							simVars.misc.verbosity
						);
					}
#if 0
					else if (helmholtz_solver_id == 5)
					{
						if (simVars.misc.verbosity > 1)
							std::cout << "REXI: Using helmholtz_iterative_smoother_conjugate_gradient_REAL" << std::endl;

						retval = RexiSWE_HelmholtzSolver::smoother_conjugate_gradient_real(	// CONJUGATE GRADIENT ITERATIVE SOLVER
							kappa,
							simVars.sim.gravitation * simVars.setup.h0,
							rhs,
							h,
							op,
							simVars.sim.domain_size,
							helmholtz_solver_eps,
							iter_max,
							-123,
							simVars.misc.verbosity
						);
					}
#endif
					else
					{
						std::cout << "INVALID SOLVER "<< std::endl;
						exit(-1);
					}

				watch.stop();

				h = h*tau;

				double residual = RexiSWE_HelmholtzSolver::helmholtz_iterative_get_residual_rms(
						kappa,
						simVars.sim.gravitation * simVars.setup.h0,
						rhs,
						h,
						simVars.sim.domain_size
					);

				// compare with analytical solution (e.g. for convergence test)
				double error_analytical = (h_cart-h).reduce_rms();

				// compare with numerical solution (has to be close to requested numerical accuracy)
				double error_numerical = ((h * kappa - simVars.sim.gravitation * simVars.setup.h0 * (h.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C))) - rhs*tau).reduce_rms();

				std::cout << "    Computed solution with residual " << residual << " in " << watch() << " seconds" << std::endl;
				std::cout << "                           RMS error on h (analytical): " << error_analytical << std::endl;
				std::cout << "                           RMS error on residual (numeric): " << error_numerical << std::endl;

				if (!retval)
				{
					std::cout << "CONVERGENCE NOT REACHED after " << iter_max << " iterations!!! STOP" << std::endl;
					exit(-1);
				}
			}
		}
		std::cout << "TEST REXI: DONE" << std::endl;
	}


	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
