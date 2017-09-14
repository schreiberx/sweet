/*
 * AppTestSPHSolverComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_TESTSPH_REXI_SOLVER_COMPLEX_HPP_
#define SRC_TESTSPH_REXI_SOLVER_COMPLEX_HPP_


#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks_sphere/SphereTestSolutions_SPH.hpp>
#include <rexi/REXI_Terry.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>

#include "../programs/swe_sphere_rexi/SWE_Sphere_TS_l_rexi.hpp"

#include <sweet/sphere/GenerateConsistentGradDivSphereData.hpp>
#include <sweet/sphere/ErrorCheck.hpp>


SimulationVariables simVars;

SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

SphereDataConfig sphereDataConfigRexiAddedModes;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigRexiAddedModes;



/**
 * Run with
 *
 * 	$ ./build/sh_example T32 P2
 */
void run_tests()
{
	std::cerr << "WARNING: THESE TESTS ARE FORMULATED FOR VELOCITY, BUT VORTICITY/DIVERGENCE is used" << std::endl;
	std::cerr << "WARNING: THESE TESTS ARE FORMULATED FOR VELOCITY, BUT VORTICITY/DIVERGENCE is used" << std::endl;
	std::cerr << "WARNING: THESE TESTS ARE FORMULATED FOR VELOCITY, BUT VORTICITY/DIVERGENCE is used" << std::endl;
	std::cerr << "Therefore, it's deactivated" << std::endl;
	return;

	simVars.outputConfig();

	double max_scalar = std::max(std::max(simVars.sim.gravitation, simVars.sim.h0), simVars.sim.f0);

	double epsilon = 1e-10*max_scalar;
	epsilon *= std::sqrt(sphereDataConfig->spectral_modes_n_max);

	std::cout << std::setprecision(20);

	if (sphereDataConfig->spectral_modes_n_max < 32)
	{
		std::cerr << "WARNING: AT LEAST 32 MODES REQUIRED for proper accuracy!!!" << std::endl;
	}


	double r = simVars.sim.earth_radius;
	double ir = 1.0/r;
	double two_omega = 2.0*simVars.sim.coriolis_omega;


	if (simVars.rexi.use_sphere_extended_modes == 0)
	{
		sphereDataConfigExt = sphereDataConfig;
	}
	else
	{
		// Add modes only along latitude since these are the "problematic" modes
		sphereDataConfigRexiAddedModes.setupAdditionalModes(
				sphereDataConfig,
				simVars.rexi.use_sphere_extended_modes,	// TODO: Extend SPH wrapper to also support m != n to set this guy to 0
				simVars.rexi.use_sphere_extended_modes
		);
		sphereDataConfigExt = &sphereDataConfigRexiAddedModes;
	}

	SphereOperatorsComplex opComplex(sphereDataConfig, 1);
	SphereOperatorsComplex opComplexExt(sphereDataConfigExt, 1);

	REXI_Terry<> rexi("phi0", simVars.rexi.h, simVars.rexi.M);

	if (!simVars.misc.sphere_use_robert_functions)
		FatalError("Only Robert formulation allowed");

	double timestep_size = 1.0;

	std::cout << "Using generic max allowed error value of " << epsilon << " (problem specific scaled)" << std::endl;


	std::complex<double> beta(1.0, 0.0);

	for (int use_complex_valued_solver = 1; use_complex_valued_solver < 2; use_complex_valued_solver++)
	{
		SphereDataComplex prog_phi_cplx_setup(sphereDataConfig);
		SphereDataComplex prog_u_cplx_setup(sphereDataConfig);
		SphereDataComplex prog_v_cplx_setup(sphereDataConfig);

		if (simVars.setup.benchmark_scenario_id <= 0)
		{
			std::cout << "SETUP: Computing solution based on time stepping scheme." << std::endl;
			SphereOperators op(sphereDataConfig, 1);

			GenerateConsistentGradDivSphereData g_real(
					simVars,
					sphereDataConfig,
					op,
					M_PI/3, M_PI/3,
					std::string("gen_divgrad_data_re_")+sphereDataConfig->getUniqueIDString()
			);

			GenerateConsistentGradDivSphereData g_imag(
					simVars,
					sphereDataConfig,
					op,
					M_PI/4, M_PI/2,
					std::string("gen_divgrad_data_im_")+sphereDataConfig->getUniqueIDString()
			);

			g_real.generate();
			g_imag.generate();

			prog_phi_cplx_setup = Convert_SphereData_To_SphereDataComplex::physical_convert(g_real.prog_h*simVars.sim.gravitation);
			prog_u_cplx_setup = Convert_SphereData_To_SphereDataComplex::physical_convert(g_real.prog_u);
			prog_v_cplx_setup = Convert_SphereData_To_SphereDataComplex::physical_convert(g_real.prog_v);


			// Combine real and imaginary data
			prog_phi_cplx_setup = prog_phi_cplx_setup + Convert_SphereData_To_SphereDataComplex::physical_convert(g_imag.prog_h*simVars.sim.gravitation)*std::complex<double>(0,1);
			prog_u_cplx_setup = prog_u_cplx_setup + Convert_SphereData_To_SphereDataComplex::physical_convert(g_imag.prog_u)*std::complex<double>(0,1);
			prog_v_cplx_setup = prog_v_cplx_setup + Convert_SphereData_To_SphereDataComplex::physical_convert(g_imag.prog_v)*std::complex<double>(0,1);

			prog_phi_cplx_setup.spectral_truncate();
			prog_u_cplx_setup.spectral_truncate();
			prog_v_cplx_setup.spectral_truncate();
		}

		for (int i = rexi.alpha.size()-1; i >= 0; i--)
		{
			std::complex<double> &alpha = rexi.alpha[i];

			std::cout << std::endl;
			std::cout << "******************************************************" << std::endl;
			std::cout << "* ALPHA: " << alpha << std::endl;
			std::cout << "******************************************************" << std::endl;

			if (!use_complex_valued_solver)
			{
				// the real value of alpha is constant!
				if (i > 0)
					break;

				// eliminate imaginary value to generate real-valued RHS
				alpha.imag(0);
			}

			/*
			 * Manufactured solution
			 */
			SphereDataComplex prog_phi_cplx(sphereDataConfig);
			SphereDataComplex prog_u_cplx(sphereDataConfig);
			SphereDataComplex prog_v_cplx(sphereDataConfig);


			if (simVars.setup.benchmark_scenario_id <= 0)
			{
				prog_phi_cplx = prog_phi_cplx_setup;
				prog_u_cplx = prog_u_cplx_setup;
				prog_v_cplx = prog_v_cplx_setup;
			}
			else if (simVars.setup.benchmark_scenario_id == 1)
			{
				std::cout << "SETUP: Computing steady state solution" << std::endl;

				prog_u_cplx.physical_update_lambda(
						[&](double i_lon, double i_lat, std::complex<double> &io_data)
						{
							io_data = std::cos(i_lat)*std::cos(i_lat);
						}
				);

				prog_v_cplx.spectral_set_zero();

				prog_phi_cplx.physical_update_lambda(
						[&](double i_lon, double i_lat, std::complex<double> &io_data)
						{
							io_data = (simVars.sim.earth_radius*simVars.sim.coriolis_omega)*(std::cos(i_lat)*std::cos(i_lat));///simVars.sim.gravitation;
						}
				);

				prog_phi_cplx *= alpha;
				prog_u_cplx *= alpha;
				prog_v_cplx *= alpha;

				prog_phi_cplx.spectral_truncate();
				prog_u_cplx.spectral_truncate();
				prog_v_cplx.spectral_truncate();
			}
			else
			{
				FatalError("Benchmark scenario not chosen");
			}

			SphereDataComplex prog_phi_cplx_ext = prog_phi_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);
			SphereDataComplex prog_u_cplx_ext = prog_u_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);
			SphereDataComplex prog_v_cplx_ext = prog_v_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);


			/*
			 * Computed right hand side
			 */
			double phi_bar = simVars.sim.gravitation*simVars.sim.h0;

			SphereDataComplex prog_phi0_cplx_ext(sphereDataConfigExt);
			SphereDataComplex prog_u0_cplx_ext(sphereDataConfigExt);
			SphereDataComplex prog_v0_cplx_ext(sphereDataConfigExt);


			/*
			 * Compute RHS for REXI term
			 */
			{
				alpha.imag(0);
				prog_phi0_cplx_ext =
							+ alpha*prog_phi_cplx_ext
							- phi_bar*ir*opComplexExt.robert_div_lon(prog_u_cplx_ext)
							- phi_bar*ir*opComplexExt.robert_div_lat(prog_v_cplx_ext);

				prog_u0_cplx_ext =
							- ir*opComplexExt.robert_grad_lon(prog_phi_cplx_ext)
							+ alpha*prog_u_cplx_ext
							+ two_omega*opComplexExt.mu(prog_v_cplx_ext);

				prog_v0_cplx_ext =
							- ir*opComplexExt.robert_grad_lat(prog_phi_cplx_ext)
							- two_omega*opComplexExt.mu(prog_u_cplx_ext)
							+ alpha*prog_v_cplx_ext;
			}

#if 0
			prog_phi0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfig);
			prog_u0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfig);
			prog_v0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfig);
#endif

			SphereDataComplex prog_phi0_cplx = prog_phi0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereDataComplex prog_u0_cplx = prog_u0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereDataComplex prog_v0_cplx = prog_v0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);

			if (simVars.setup.benchmark_scenario_id == 1)
			{
				SphereDataComplex zero(sphereDataConfig);
				zero.physical_set_zero();

				// Check for geostrophic balance
				ErrorCheck::checkTruncated(prog_phi0_cplx_ext, prog_phi_cplx_ext*alpha, sphereDataConfig, "ERROR Geostrophic balance phi", epsilon, false);
				ErrorCheck::checkTruncated(prog_u0_cplx_ext, prog_u_cplx_ext*alpha, sphereDataConfig, "ERROR Geostrophic balance u", epsilon, false);
				ErrorCheck::checkTruncated(prog_v0_cplx_ext, prog_v_cplx_ext*alpha, sphereDataConfig, "ERROR Geostrophic balance v", epsilon, false);

				ErrorCheck::checkTruncated(prog_v0_cplx_ext, zero, sphereDataConfig, "ERROR Geostrophic balance v", epsilon, false);
			}

			if (simVars.setup.benchmark_scenario_id == 1)
			{
				SphereDataComplex zero(sphereDataConfig);
				zero.physical_set_zero();


				SphereDataComplex f(sphereDataConfig);
				f.physical_update_lambda_gaussian_grid(
						[&](double lon, double mu, std::complex<double> &o_data)
						{
							o_data = mu*two_omega;
						}
					);


				SphereDataComplex data_inv_f(sphereDataConfig);
				data_inv_f.physical_update_lambda_gaussian_grid(
						[&](double lon, double mu, std::complex<double> &o_data)
						{
							o_data = 1.0/(mu*two_omega);
						}
					);
				auto inv_f = [&](const SphereDataComplex &i_data)	-> SphereDataComplex
				{
					return data_inv_f*i_data;
				};


				SphereDataComplex data_inv_cos2phi(sphereDataConfig);
				data_inv_cos2phi.physical_update_lambda(
						[&](double lon, double phi, std::complex<double> &o_data)
						{
							o_data = 1.0/(cos(phi)*cos(phi));
						}
					);
				auto inv_cos2phi = [&](const SphereDataComplex &i_data)	-> SphereDataComplex
				{
					return data_inv_f*i_data;
				};


				{
					/*
					 * Test equations in geostrophic balance test case
					 * Section 5.2, SPREXI ver 15
					 */
					{
						SphereDataComplex lhs = -opComplexExt.robert_div_lon(prog_u_cplx) - opComplexExt.robert_div_lat(prog_v_cplx);
						ErrorCheck::checkTruncated(lhs, zero, sphereDataConfig, "ERROR Geostrophic balance test a", epsilon, false, false);
					}

					{
						// Equation (1)
						SphereDataComplex lhs = -opComplexExt.robert_grad_lon(prog_phi_cplx) + f*prog_v_cplx;
						SphereDataComplex rhs = zero;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 1", epsilon, false, false);
					}

					{
						// Equation (1b)
						SphereDataComplex lhs = prog_v_cplx;
						SphereDataComplex rhs = inv_f(opComplexExt.robert_grad_lon(prog_phi_cplx));
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 1b", epsilon, false, false);
					}

					{
						// Equation (2)
						SphereDataComplex lhs = -opComplexExt.robert_grad_lat(prog_phi_cplx) - f*prog_u_cplx;
						SphereDataComplex rhs = zero;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 2", epsilon, false, false);
					}

					{
						// Equation (2b)
						SphereDataComplex lhs = prog_u_cplx;
						SphereDataComplex rhs = -inv_f(opComplexExt.robert_grad_lat(prog_phi_cplx));
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 2b", epsilon, false, false);
					}

					{
						SphereDataComplex lhs =
								opComplexExt.robert_div_lon(-inv_f(opComplexExt.robert_grad_lat(prog_phi_cplx)))
								+ opComplexExt.robert_div_lat(inv_f(opComplexExt.robert_grad_lon(prog_phi_cplx)));
						SphereDataComplex rhs = zero;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 3", epsilon, false, false);
					}

					{
						SphereDataComplex lhs =
								-opComplexExt.robert_grad_lat(prog_phi_cplx)*inv_cos2phi(opComplexExt.robert_grad_lon(data_inv_f))
								+opComplexExt.robert_grad_lon(prog_phi_cplx)*inv_cos2phi(opComplexExt.robert_grad_lat(data_inv_f));
						SphereDataComplex rhs = zero;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 4", epsilon, false, false);
					}
#if 0
					{
						SphereDataComplex lhs =
								-inv_cos2phi(opComplexExt.robert_grad_lon(prog_phi_cplx))
								+f*prog_v_cplx;
						SphereDataComplex rhs = zero;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 5a", epsilon, false, false);
					}

					{
						SphereDataComplex lhs = inv_cos2phi(opComplexExt.robert_grad_lat(prog_phi_cplx));
						SphereDataComplex rhs = -two_omega*opComplexExt.mu(prog_u_cplx);
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test 5b", epsilon, false, false);
					}
#endif

					{
						SphereDataComplex lhs =
								//inv_cos2phi
								(
										opComplexExt.robert_grad_lat(prog_phi_cplx)
								)
								;
						SphereDataComplex rhs = -f*prog_u_cplx;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test robert", epsilon, false, false);
					}
				}
			}

			{
				/*
				 * Test correct mode extensions for sphere data with complex physical data
				 */
				SphereDataComplex lhs = prog_u0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);
				SphereDataComplex rhs = prog_u0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfigExt).spectral_returnWithDifferentModes(sphereDataConfig);

				SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
				SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

				std::cout << "+ error for mode extensions: " << (lhsr-rhsr).physical_reduce_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_rms() << std::endl;

				if ((lhsr-rhsr).physical_reduce_max_abs() > 10e-10)
					FatalError("Threshold exceeded");
			}

			SphereData prog_phi0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert_real(prog_phi0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigExt);
			SphereData prog_u0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert_real(prog_u0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigExt);
			SphereData prog_v0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert_real(prog_v0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigExt);

			SphereData prog_phi0 = prog_phi0_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereData prog_u0 = prog_u0_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereData prog_v0 = prog_v0_ext.spectral_returnWithDifferentModes(sphereDataConfig);



			SphereDataComplex f(sphereDataConfigExt);
			f.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, std::complex<double> &o_data)
					{
						o_data = mu*two_omega;
					}
				);


			/*
			 * Output of REXI (ext modes)
			 */

			SphereData rexi_prog_phi_ext(sphereDataConfigExt);
			SphereData rexi_prog_u_ext(sphereDataConfigExt);
			SphereData rexi_prog_v_ext(sphereDataConfigExt);


			SphereDataComplex rexi_prog_phi_cplx_ext(sphereDataConfigExt);
			SphereDataComplex rexi_prog_u_cplx_ext(sphereDataConfigExt);
			SphereDataComplex rexi_prog_v_cplx_ext(sphereDataConfigExt);

			if (simVars.misc.sphere_use_robert_functions)
			{
				{
					SWERexiTerm_SPHRobert rexiSPHRobert;

					rexiSPHRobert.setup_vectorinvariant_progphivortdiv(
							sphereDataConfigExt,

							alpha,
							beta,

							simVars.sim.earth_radius,
							simVars.sim.coriolis_omega,
							simVars.sim.f0,
							phi_bar,
							timestep_size,

							simVars.sim.f_sphere,
							false
					);

					if (use_complex_valued_solver)
					{
#if 1
						rexiSPHRobert.solve_vectorinvariant_progphivortdiv(
								prog_phi0_ext,
								prog_u0_ext,
								prog_v0_ext,

								rexi_prog_phi_ext,
								rexi_prog_u_ext,
								rexi_prog_v_ext
							);
#else

						rexiSPHRobert.solve_complex(
								prog_phi0_cplx_ext,
								prog_u0_cplx_ext,
								prog_v0_cplx_ext,

								rexi_prog_phi_cplx_ext,
								rexi_prog_u_cplx_ext,
								rexi_prog_v_cplx_ext
							);
						ErrorCheck::checkTruncated(rexi_prog_phi_cplx_ext, prog_phi_cplx_ext, sphereDataConfig, "prog_phi SSSSSSSSSSS", epsilon*1e+1);
#endif
					}
					else
					{
						FatalError("");

						rexiSPHRobert.solve_vectorinvariant_progphivortdiv(
								prog_phi0_ext,
								prog_u0_ext,
								prog_v0_ext,

								rexi_prog_phi_ext,
								rexi_prog_u_ext,
								rexi_prog_v_ext
							);
					}
				}

			}


			/*
			 * Output of REXI solver (no ext modes)
			 */
			SphereData rexi_prog_phi = rexi_prog_phi_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereData rexi_prog_u = rexi_prog_u_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereData rexi_prog_v = rexi_prog_v_ext.spectral_returnWithDifferentModes(sphereDataConfig);


			/*
			 * Test solution
			 */
			SphereData prog_phi = Convert_SphereDataComplex_To_SphereData::physical_convert_real(prog_phi_cplx);
			SphereData prog_u = Convert_SphereDataComplex_To_SphereData::physical_convert_real(prog_u_cplx);
			SphereData prog_v = Convert_SphereDataComplex_To_SphereData::physical_convert_real(prog_v_cplx);


			/*
			 * REXI results stored in rexi_prog_*
			 * These should match prog_*
			 */
// commented out since this is not how REXI is computed
//			ErrorCheck::checkTruncated(rexi_prog_phi, prog_phi, sphereDataConfig, "prog_phi", epsilon*1e+5, false, simVars.setup.benchmark_scenario_id != 1);
//			ErrorCheck::checkTruncated(rexi_prog_u, prog_u, sphereDataConfig, "prog_u", epsilon*1e+5, false, simVars.setup.benchmark_scenario_id != 1);
//			ErrorCheck::checkTruncated(rexi_prog_v, prog_v, sphereDataConfig, "prog_v", epsilon*1e+5, false, simVars.setup.benchmark_scenario_id != 1);
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

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 0;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{

#if SWEET_PARAREAL
		simVars.parareal.printOptions();
#endif
		return -1;
	}

	if (simVars.disc.res_spectral[0] == 0)
		FatalError("Set number of spectral modes to use SPH!");

	sphereDataConfig->setupAutoPhysicalSpace(
					simVars.disc.res_spectral[0],
					simVars.disc.res_spectral[1],
					&simVars.disc.res_physical[0],
					&simVars.disc.res_physical[1]
			);

	run_tests();
}



#endif /* SRC_TESTSPHSOLVERS_HPP_ */
