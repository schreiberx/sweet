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
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>

#include <sweet/sphere/app_swe/SWERexiTerm_SPH.hpp>
#include <sweet/sphere/app_swe/SWERexiTerm_SPHRobert.hpp>

#include <sweet/sphere/GenerateConsistentGradDivSphereData.hpp>

#include <rexi/REXI.hpp>
#include <rexi/swe_sphere_rexi/SWE_Sphere_REXI.hpp>

#include <sweet/sphere/ErrorCheck.hpp>


SimulationVariables simVars;

SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

SphereDataConfig sphereDataConfigRexiAddedModes;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigRexiAddedModes;


bool param_use_coriolis_formulation = true;




/**
 * Run with
 *
 * 	$ ./build/sh_example T32 P2
 */
void run_tests()
{
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


	if (simVars.rexi.rexi_use_extended_modes == 0)
	{
		sphereDataConfigExt = sphereDataConfig;
	}
	else
	{
		// Add modes only along latitude since these are the "problematic" modes
		sphereDataConfigRexiAddedModes.setupAdditionalModes(
				sphereDataConfig,
				simVars.rexi.rexi_use_extended_modes,	// TODO: Extend SPH wrapper to also support m != n to set this guy to 0
				simVars.rexi.rexi_use_extended_modes
		);
		sphereDataConfigExt = &sphereDataConfigRexiAddedModes;
	}

	SphereOperatorsComplex opComplex(sphereDataConfig);
	SphereOperatorsComplex opComplexExt(sphereDataConfigExt);

	REXI<> rexi(0, simVars.rexi.rexi_h, simVars.rexi.rexi_M);

	if (!simVars.misc.sphere_use_robert_functions)
	{
		std::cerr << "WARNING: TESTS WILL FAIL WITH NON-ROBERT FUNCTIONS DUE TO SINGULARITY AT POLES!" << std::endl;
	}

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
			SphereOperators op(sphereDataConfig);

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

				if (simVars.misc.sphere_use_robert_functions)
				{
					prog_u_cplx.physical_update_lambda(
							[&](double i_lon, double i_lat, std::complex<double> &io_data)
							{
								io_data = std::cos(i_lat)*std::cos(i_lat);
							}
					);
				}
				else
				{
					prog_u_cplx.physical_update_lambda(
							[&](double i_lon, double i_lat, std::complex<double> &io_data)
							{
								io_data = std::cos(i_lat);
							}
					);
				}

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
			if (simVars.misc.sphere_use_robert_functions)
			{
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
			else
			{
				prog_phi0_cplx_ext =
							+ alpha*prog_phi_cplx_ext
							- phi_bar*ir*opComplexExt.div_lon(prog_u_cplx_ext)
							- phi_bar*ir*opComplexExt.div_lat(prog_v_cplx_ext);

				prog_u0_cplx_ext =
							- ir*opComplexExt.grad_lon(prog_phi_cplx_ext)
							+ alpha*prog_u_cplx_ext
							+ two_omega*opComplexExt.mu(prog_v_cplx_ext);

				prog_v0_cplx_ext =
							- ir*opComplexExt.grad_lat(prog_phi_cplx_ext)
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


				if (simVars.misc.sphere_use_robert_functions)
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

			SphereData prog_phi0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_phi0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigExt);
			SphereData prog_u0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_u0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigExt);
			SphereData prog_v0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_v0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigExt);

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


#if 1
			if (!simVars.misc.sphere_use_robert_functions)
			{
				SphereDataComplex phi(sphereDataConfigExt);
				SphereDataComplex u(sphereDataConfigExt);
				SphereDataComplex v(sphereDataConfigExt);


				phi = prog_phi_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);
				u = prog_u_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);
				v = prog_v_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);


//				SphereDataComplex &phi0 = prog_phi0_cplx_ext;
				SphereDataComplex &u0 = prog_u0_cplx_ext;
				SphereDataComplex &v0 = prog_v0_cplx_ext;


				SphereDataComplex div0 = ir*opComplexExt.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplexExt.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplexExt.robert_div(u, v);
				SphereDataComplex eta = ir*opComplexExt.robert_vort(u, v);

				SphereDataComplex fi = ir*opComplexExt.robert_grad_lon(f);
				SphereDataComplex fj = ir*opComplexExt.robert_grad_lat(f);

				SphereDataComplex lhs(sphereDataConfigExt);
				SphereDataComplex rhs(sphereDataConfigExt);

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Check Identity relation
					 *
					 * D.(fA) = f(D.A) + A.(Gf)
					 */

					lhs = ir*opComplexExt.div(f*u, f*v);

					rhs =	  f*ir*opComplexExt.div(u, v)
							+ u*ir*opComplexExt.grad_lon(f)
							+ v*ir*opComplexExt.grad_lat(f);

					lhs.physical_truncate();
					rhs.physical_truncate();

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Identity for non-Robert formulation", epsilon);
				}
			}

			/*
			 * TEST formulations
			 *
			 * assume that 2\Omega = 1.0
			 */
			if (simVars.misc.sphere_use_robert_functions)
			{
				SphereDataComplex &phi0 = prog_phi0_cplx_ext;
				SphereDataComplex &u0 = prog_u0_cplx_ext;
				SphereDataComplex &v0 = prog_v0_cplx_ext;

				SphereDataComplex phi(sphereDataConfigExt);
				SphereDataComplex u(sphereDataConfigExt);
				SphereDataComplex v(sphereDataConfigExt);


				phi = prog_phi_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);
				u = prog_u_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);
				v = prog_v_cplx.spectral_returnWithDifferentModes(sphereDataConfigExt);


				SphereDataComplex div0 = ir*opComplexExt.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplexExt.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplexExt.robert_div(u, v);
				SphereDataComplex eta = ir*opComplexExt.robert_vort(u, v);


				if (simVars.setup.benchmark_scenario_id == 1)
				{
					/*
					 * Test for geostrophic balance
					 */
					{
						SphereDataComplex lhs = prog_phi_cplx*alpha;
						SphereDataComplex rhs = prog_phi0_cplx;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test alpha*phi", epsilon, false);
					}
					{
						SphereDataComplex lhs = prog_u_cplx*alpha;
						SphereDataComplex rhs = prog_u0_cplx;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test alpha*u", epsilon, false);
					}
					{
						SphereDataComplex lhs = prog_v_cplx*alpha;
						SphereDataComplex rhs = prog_v0_cplx;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test alpha*v", epsilon, false);
					}
					{
						SphereDataComplex lhs = eta*alpha;
						SphereDataComplex rhs = eta0;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test alpha*eta", epsilon, false);
					}
					{
						SphereDataComplex lhs = div*alpha;
						SphereDataComplex rhs = div0;
						ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "ERROR Geostrophic balance test alpha*div", epsilon, false);
					}
				}


				/*
				 * M(\mu) term for
				 * 	D.(f A) = f(D.A) + M(\my)*A.grad(f)
				 * fix
				 */
				SphereDataComplex one_over_cos2phi(sphereDataConfigExt);
				one_over_cos2phi.physical_update_lambda_cosphi_grid(
						[&](double lon, double cosphi, std::complex<double> &o_data)
						{
							o_data = 1.0/(cosphi*cosphi);
						}
				);

#if 0
				SphereDataComplex fi(sphereDataConfigExt);
				fi.physical_set_zero();

				SphereDataComplex fj(sphereDataConfigExt);
				fj.physical_update_lambda_cosphi_grid(
						[&](double lon, double cosphi, std::complex<double> &o_data)
						{
							o_data = ir*two_omega;
						}
					);

#else
				double fi = 0;
				double fj = ir*two_omega;
#endif

				SphereDataComplex lhs(sphereDataConfigExt);
				SphereDataComplex rhs(sphereDataConfigExt);

#if 1

#if 0
				{
					lhs = v;
/*
					lhs.physical_update_lambda_gaussian_grid(
							[](double lambda, double mu, std::complex<double> &o_data)
							{
								o_data *= (1.0-mu*mu);
							}
					);
*/

					lhs = opComplexExt.sphSolver_inv_one_minus_mu2.solve(lhs);


					rhs = v;
/*
					rhs.physical_update_lambda_gaussian_grid(
							[](double lambda, double mu, std::complex<double> &o_data)
							{
								o_data *= (1.0-mu*mu);
							}
					);
*/
					rhs.physical_update_lambda_gaussian_grid(
							[](double lambda, double mu, std::complex<double> &o_data)
							{
								o_data /= (1.0-mu*mu);
							}
					);
//					rhs.spectral_truncate();

					ErrorCheck::checkTruncatedSpectral(lhs, rhs, sphereDataConfig, "TEST A", epsilon);
				}
#endif

#if 0
				{
					lhs = u;
/*
					lhs.physical_update_lambda_gaussian_grid(
							[](double lambda, double mu, std::complex<double> &o_data)
							{
								o_data *= (1.0-mu*mu);
							}
					);
*/

					lhs = opComplexExt.sphSolver_inv_one_minus_mu2.solve(lhs);


					rhs = u;
/*
					rhs.physical_update_lambda_gaussian_grid(
							[](double lambda, double mu, std::complex<double> &o_data)
							{
								o_data *= (1.0-mu*mu);
							}
					);
*/
					rhs.physical_update_lambda_gaussian_grid(
							[](double lambda, double mu, std::complex<double> &o_data)
							{
								o_data /= (1.0-mu*mu);
							}
					);

					ErrorCheck::checkTruncatedSpectral(lhs, rhs, sphereDataConfig, "TEST A", epsilon);
				}
				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Check Identity relation
					 *
					 * D.(fA) = f(D.A) + A.(Df)
					 */
					lhs = opComplexExt.robert_div(f*u, f*v);

					rhs =	f*opComplexExt.robert_div(u, v)
							+ opComplexExt.robert_grad_M(f, u, v);

					ErrorCheck::checkTruncatedSpectral(lhs, rhs, sphereDataConfig, "ERROR IDENTITY D.(fA)", epsilon, simVars.rexi.rexi_use_extended_modes <= 2);
				}
#endif

				{
					lhs = opComplexExt.laplace(phi);

					rhs =	  opComplexExt.robert_div_lat(opComplexExt.robert_grad_lat(phi))
							+ opComplexExt.robert_div_lon(opComplexExt.robert_grad_lon(phi))
						;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "Laplace", epsilon, simVars.rexi.rexi_use_extended_modes <= 2);
				}


				if (simVars.rexi.rexi_use_extended_modes == 0)
				{
					lhs = opComplexExt.robert_div_lat(v.spectral_returnWithDifferentModes(sphereDataConfig)).spectral_returnWithDifferentModes(sphereDataConfigExt);

					rhs = (
							opComplexExt.robert_div_lat(prog_v_cplx)
						).spectral_returnWithDifferentModes(sphereDataConfigExt);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "0", epsilon, simVars.rexi.rexi_use_extended_modes <= 2);
				}

				{
					/*
					 * REXI SPH document ver 15, Section 5.3, eq (2)
					 */

					lhs = -ir*opComplexExt.robert_grad_lon(phi) + alpha*u + f*v;
					rhs = u0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "0b", epsilon, false);
				}

				{
					/*
					 * REXI SPH document ver 15, Section 5.3, eq (3)
					 */

					lhs = -ir*opComplexExt.robert_grad_lat(phi) - f*u + alpha*v;
					rhs = v0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "0c", epsilon, false);
				}

				{
					lhs = div;
					rhs = ir*opComplexExt.robert_div_lon(u) + ir*opComplexExt.robert_div_lat(v);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "0d", epsilon, false);
				}

				{
					lhs = eta;
					rhs = ir*opComplexExt.robert_div_lon(v) - ir*opComplexExt.robert_div_lat(u);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "0e", epsilon, false);
				}

#if 0
				{
					/*
					 * REXI SPH document ver 15, Section 5.3.1
					 * Check Identity relation
					 *
					 * D.(fA) = f(D.A) + A.one_over_cos2phi*(Df)
					 */

					lhs = ir*opComplexExt.robert_div(opComplexExt.mu(u)*two_omega, opComplexExt.mu(v)*two_omega);

					rhs =	opComplexExt.mu(ir*opComplexExt.robert_div(u, v))*two_omega
							+ u*fi + v*fj;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "0f", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}
#endif

				{
					/*
					 * REXI SPH document ver 15, Section 5.3.1
					 * div lon (2)
					 */
					lhs = ir*opComplexExt.robert_div_lon(
								-ir*opComplexExt.robert_grad_lon(phi) + alpha*u + two_omega*opComplexExt.mu(v)
							  );

					rhs = ir*opComplexExt.robert_div_lon(u0);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1aa", epsilon, false);
				}

				{
					/*
					 * REXI SPH document ver 15, Section 5.3.1
					 * div lat (3)
					 */
#if 0
					lhs = ir*opComplexExt.robert_div_lat(v0);

					rhs = ir*opComplexExt.robert_div_lat(
							-ir*opComplexExt.robert_grad_lat(phi)
							-f*u
							+ alpha*v
						  );
#else
					// we use this formulation for the error checks since div by cos2phi would amplify errors
					lhs = ir*opComplexExt.robert_cos2phi_div_lat(v0);

					rhs = ir*opComplexExt.robert_cos2phi_div_lat(
							-ir*opComplexExt.robert_grad_lat(phi)
							-two_omega*opComplexExt.mu(u)
							+ alpha*v
						  );
#endif
					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1ab", epsilon, false);
				}

				{
					/*
					 * REXI SPH document ver 15, Section 5.3.1
					 */

					lhs = div0;

					rhs =	ir*opComplexExt.robert_div_lon(
								(-ir*opComplexExt.robert_grad_lon(phi) + f*v) + alpha*u
							  )
							+ ir*opComplexExt.robert_div_lat(
								(-ir*opComplexExt.robert_grad_lat(phi) - f*u) + alpha*v
							  );

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1a", epsilon, false);
				}

#if 0
				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st equation in this section
					 */

					lhs = div0;

					rhs =	  ir*opComplexExt.robert_div_lon(	-ir*opComplexExt.robert_grad_lon(phi)	)
							+ ir*opComplexExt.robert_div_lon(	alpha*u 	)
							+ ir*opComplexExt.robert_div_lon(	f*v		)

							+ ir*opComplexExt.robert_div_lat(	-ir*opComplexExt.robert_grad_lat(phi)	)
							+ ir*opComplexExt.robert_div_lat(	+ alpha*v 	)
							+ ir*opComplexExt.robert_div_lat(	-f*u	)
							;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1ba", epsilon, false, true);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st equation in this section
					 */

					lhs = div0;

					rhs =	 -ir*ir*opComplexExt.laplace(phi)

							+ alpha*ir*opComplexExt.robert_div(u, v)

							+ ir*opComplexExt.robert_div(	f*v	, -f*u	)
							;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1bb", epsilon, simVars.rexi.rexi_use_extended_modes < 2, true);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Test identity of D.(fA) with A=(v,-u) formulation
					 */

					lhs = opComplexExt.robert_div(	opComplexExt.mu(v)	, -opComplexExt.mu(u)	);

					rhs =	  opComplexExt.mu(opComplexExt.robert_div(v, -u))
								+
								opComplexExt.inv_one_minus_mu2
								(
									v*opComplexExt.robert_grad_lon(f)
									- u*opComplexExt.robert_grad_lat(f)
								);

//					ErrorCheck::check(lhs, rhs, "1ca", epsilon, simVars.rexi.rexi_use_extended_modes < 2, true);
					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1ca", epsilon, simVars.rexi.rexi_use_extended_modes < 2, true);
				}


				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Test identity of D.(fA) with A=(v,-u) formulation
					 */

					lhs = ir*opComplexExt.robert_div(	f*v	, -f*u	);

					rhs =	  f*ir*opComplexExt.robert_div(v, -u)
								+ opComplexExt.inv_one_minus_mu2
								(
									v*ir*opComplexExt.robert_grad_lon(f)
									- u*ir*opComplexExt.robert_grad_lat(f)
								);

//					ErrorCheck::check(lhs, rhs, "1ca", epsilon, simVars.rexi.rexi_use_extended_modes < 2, true);
					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1ca", epsilon, simVars.rexi.rexi_use_extended_modes < 2, true);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * eq. (4)
					 */

					lhs =	- ir*ir*opComplexExt.laplace(phi)
							+ alpha*div

							+f*eta
							+ v*fi
							- u*fj
					;

					rhs = div0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "1c", epsilon, simVars.rexi.rexi_use_extended_modes < 2, simVars.setup.benchmark_scenario_id != 1);
				}
#endif

#if 0
				{
					/*
					 * REXI SPH document ver 14, Section 5.2.2
					 * eq. (5)
					 */

					lhs = -alpha*eta + f*div + fj*v + fi*u;
					rhs = -eta0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "2", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.3
					 * eq. (6)
					 */

					lhs = div;
					rhs = 1.0/phi_bar*(phi*alpha-phi0);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "3", epsilon*10e+2, false);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.4
					 * eq. (7)
					 */

					lhs =	- ir*ir*opComplexExt.laplace(phi)
							+ 1.0/phi_bar * phi * alpha*alpha
							+ f * eta
							+ fi*v
							- fj*u;

					rhs =	div0 + 1.0/phi_bar * alpha * phi0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "4a", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					/*
					 * REXI SPH document ver 15, Section 5.3.4
					 * between eq. (7) and (8)
					 */

					lhs = -alpha*eta + 1.0/phi_bar * alpha * f * phi + fj*v + fi*u;
					rhs = -eta0 + 1.0/phi_bar * f * phi0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "4ba", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * REXI SPH document ver 15, Section 5.3.4
					 * eq. (8)
					 */

					lhs = eta;
					rhs =	1.0/alpha * eta0
							- 1.0/(phi_bar*alpha)*f*phi0
							+ 1.0/phi_bar*f*phi
							+ 1.0/alpha * fj*v
							+ 1.0/alpha * fi*u;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "5", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * REXI SPH document ver 15, Section 5.3.5
					 */

					// eq. (8)
					SphereDataComplex ceta =
							1.0/alpha * eta0
							- 1.0/(phi_bar*alpha)*f*phi0
							+ 1.0/phi_bar*f*phi
							+ 1.0/alpha * fj*v
							+ 1.0/alpha * fi*u;

					// eq. (7)
					lhs =	- ir*ir*opComplexExt.laplace(phi)
							+ 1.0/phi_bar * phi * alpha*alpha
							+ f * ceta	/* we use the computed eta here */
							+ fi*v
							- fj*u;

					rhs =	div0 + 1.0/phi_bar * alpha * phi0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "5b", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * REXI SPH document ver 15, Section 5.3.5
					 */

					lhs =
							- ir*ir*opComplexExt.laplace(phi)

							+ 1.0/phi_bar * phi * alpha*alpha

							+ f*(1.0/alpha) * eta0
							- f*(1.0/(phi_bar*alpha))*f*phi0
							+ f*(1.0/phi_bar)*f*phi
							+ f*(1.0/alpha) * fj*v
							+ f*(1.0/alpha) * fi*u

							+ fi*v
							- fj*u;

					rhs =	div0
							+ 1.0/phi_bar * alpha * phi0;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "5c", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * REXI SPH document ver 15, Section 5.3.5 at the end
					 * eq. (9)
					 */
#if 0
					SphereDataComplex F = f*(fj*v + fi*u) + alpha*(fi*v - fj*u);

					lhs = alpha*alpha*phi + f*f*phi - phi_bar*ir*ir*opComplexExt.laplace(phi) + (phi_bar/alpha) * F;
					rhs = phi_bar*(div0 - f*(1.0/alpha)*eta0) + (alpha+f*f*(1.0/alpha))*phi0;

#else

					SphereDataComplex F = two_omega*opComplexExt.mu(fj*v + fi*u) + alpha*(fi*v - fj*u);


					lhs = alpha*alpha*phi + two_omega*two_omega*opComplexExt.mu2(phi) - phi_bar*ir*ir*opComplexExt.laplace(phi) + (phi_bar/alpha) * F;
					rhs = phi_bar*(div0 - two_omega*opComplexExt.mu((1.0/alpha)*eta0)) + (alpha*phi0+(1.0/alpha)*two_omega*two_omega*opComplexExt.mu2(phi0));

#endif

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "5d", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}
#endif

				{
					/*
					 * REXI SPH document ver 15, Section 5.3.6
					 *
					 * eq. V= A^{-1} (V0 + GRAD(PHI)) = A^{-1} T
					 */

					SphereDataComplex Ti = u0 + ir*opComplexExt.robert_grad_lon(phi);
					SphereDataComplex Tj = v0 + ir*opComplexExt.robert_grad_lat(phi);

					lhs =	  alpha*alpha*u
							+ two_omega*two_omega*opComplexExt.mu2(u);

					rhs =	  alpha*Ti
							- two_omega*opComplexExt.mu(Tj);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "6 lemma u", epsilon*1e+2, simVars.rexi.rexi_use_extended_modes < 2);


					lhs =	  alpha*alpha*v
							+ two_omega*two_omega*opComplexExt.mu2(v);

					rhs =	  two_omega*opComplexExt.mu(Ti)
							+ alpha*Tj;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "6 lemma v", epsilon*1e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					/*
					 * REXI SPH document ver 15, Section 5.3.6
					 *
					 * eq. V= A^{-1} (V0 + GRAD(PHI)) = A^{-1} T
					 */
					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
					SphereDataComplex inv_kappa = one/kappa;

					SphereDataComplex Ti = u0 + ir*opComplexExt.robert_grad_lon(phi);
					SphereDataComplex Tj = v0 + ir*opComplexExt.robert_grad_lat(phi);

					lhs = u;

					rhs =	inv_kappa*(
								alpha*Ti
								- two_omega*opComplexExt.mu(Tj)
							);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "6 lemma inv_kappa u", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);


					lhs = v;

					rhs =	inv_kappa*(
							  two_omega*opComplexExt.mu(Ti)
							+ alpha*Tj
							);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "6 lemma inv_kappa v", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * REXI SPH document ver 15, Section 5.3.6
					 *
					 * eq. V= A^{-1} (V0 + GRAD(PHI)) = A^{-1} T
					 */
					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
					SphereDataComplex inv_kappa = one/kappa;

					SphereDataComplex Ti = u0 + ir*opComplexExt.robert_grad_lon(phi);
					SphereDataComplex Tj = v0 + ir*opComplexExt.robert_grad_lat(phi);

					lhs = u;

					rhs =	inv_kappa*(
								alpha*Ti
								- f*Tj
							);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "6 lemma inv_kappa f u", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);


					lhs = v;

					rhs =	inv_kappa*(
							  f*Ti
							+ alpha*Tj
							);

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "6 lemma inv_kappa f v", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					/*
					 * Test F=... formulation
					 *
					 * rexi sph ver 15, sec. 5.3.6
					 * F_lhs: first equation in this section, based on (u,v)
					 * F_rhs: other equation only based on (phi, u0, v0)
					 */
					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
					SphereDataComplex inv_kappa = one/kappa;

					SphereDataComplex F_lhs =
									  f*(fj*v + fi*u)
									+ alpha*(fi*v - fj*u);

					SphereDataComplex Fp_i = fj*(-(alpha*alpha-f*f));
					SphereDataComplex Fp_j = fj*(2.0*alpha*f);

					SphereDataComplex F_rhs = inv_kappa*(
							  Fp_i*(u0 + ir*opComplexExt.robert_grad_lon(phi) )
							+ Fp_j*(v0 + ir*opComplexExt.robert_grad_lat(phi) )
						);

					lhs = F_lhs;
					rhs = F_rhs;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "6 F(u,v) = F(u0,v0,phi)", epsilon*10e+1, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * Test F=... formulation
					 *
					 * rexi sph ver 15, sec. 5.3.7
					 * Solution for geopotential
					 */

					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
					SphereDataComplex inv_kappa = one/kappa;

					SphereDataComplex Fp_i = inv_kappa*fj*(-(alpha*alpha-f*f));
					SphereDataComplex Fp_j = inv_kappa*fj*(2.0*alpha*f);

					SphereDataComplex Fc = Fp_i*u0 + Fp_j*v0;

					SphereDataComplex F_lhs =
							  kappa*phi
							+ phi_bar/alpha*(Fp_i*ir*opComplexExt.robert_grad_lon(phi) + Fp_j*ir*opComplexExt.robert_grad_lat(phi))
							- phi_bar*ir*ir*opComplexExt.laplace(phi);

					SphereDataComplex F_rhs =
							  phi_bar*(div0 - f*(1.0/alpha)*eta0)
							+ (alpha + f*f*(1.0/alpha))*phi0
							- phi_bar/alpha*Fc;

					lhs = F_lhs;
					rhs = F_rhs;

					if (ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "7a", epsilon, true, simVars.setup.benchmark_scenario_id != 1))
						std::cout << "ERRORS IGNORED WITH 1/kappa formulation!" << std::endl;
				}

#if 0
				{
					/*
					 * Test F=... formulation
					 *
					 * rexi sph ver 15, sec. 5.3.7
					 * Solution for geopotential
					 */

					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
					SphereDataComplex inv_kappa = one/kappa;

					SphereDataComplex Fp_i = inv_kappa*fj*(-(alpha*alpha-f*f));
					SphereDataComplex Fp_j = inv_kappa*fj*(2.0*alpha*f);

					SphereDataComplex Fc = Fp_i*u0 + Fp_j*v0;

					SphereDataComplex F_lhs =
							  kappa*phi
							+ phi_bar/alpha*(Fp_i*ir*opComplexExt.robert_grad_lon(phi) + Fp_j*ir*opComplexExt.robert_grad_lat(phi))
							- phi_bar*ir*ir*opComplexExt.laplace(phi);

					SphereDataComplex F_rhs =
							  phi_bar*(div0 - f*(1.0/alpha)*eta0)
							+ (alpha + f*f*(1.0/alpha))*phi0
							- phi_bar/alpha*Fc;

					// Multiply with mu2 to generate modes which are not representable
					lhs = opComplexExt.mu2(F_lhs);
					rhs = opComplexExt.mu2(F_rhs);

					std::cout << "This should trigger an error with ext_modes <= 3!" << std::endl;
					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "7ab", epsilon, true);
				}
#endif

#if 0
				{
					/*
					 * Test F=... formulation, |*kappa scaled lhs and rhs
					 *
					 * rexi sph ver 15, sec. 5.3.7
					 * Solution for geopotential
					 */

					SphereDataComplex one(sphereDataConfigExt);
					one.physical_update_lambda(
							[&](double i_lon, double i_lat, std::complex<double> &o_data)
							{
								o_data = 1.0;
							}
					);

					SphereDataComplex kappa = alpha*alpha+f*f;

					SphereDataComplex Fkp_i = fj*(-(alpha*alpha-f*f));
					SphereDataComplex Fkp_j = fj*(2.0*alpha*f);

					SphereDataComplex Fkc = Fkp_i*u0 + Fkp_j*v0;

					SphereDataComplex lhs =
							  kappa*kappa*phi
							+ phi_bar/alpha*(Fkp_i*ir*opComplexExt.robert_grad_lon(phi) + Fkp_j*ir*opComplexExt.robert_grad_lat(phi))
							- kappa*phi_bar*ir*ir*opComplexExt.laplace(phi);

					SphereDataComplex rhs =
							kappa*phi_bar*(div0 - f*(1.0/alpha)*eta0)
							+ kappa*(alpha + f*f*(1.0/alpha))*phi0
							- phi_bar/alpha*Fkc;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "7b", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 4);
				}
#endif

				{
					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
					SphereDataComplex inv_kappa = one/kappa;

					/*
					 * REXI SPH document ver 15 Section 5.3.8 at the end
					 * Computing velocities
					 */
					lhs = kappa*u;

					SphereDataComplex a = u0 + ir*opComplexExt.robert_grad_lon(phi);
					SphereDataComplex b = v0 + ir*opComplexExt.robert_grad_lat(phi);

					rhs = alpha*a -f*b;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "8a", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
					SphereDataComplex inv_kappa = one/kappa;

					/*
					 * REXI SPH document ver 15 Section 5.3.8 at the end
					 * Computing velocities
					 */
					lhs = kappa*v;

					SphereDataComplex a = u0 + ir*opComplexExt.robert_grad_lon(phi);
					SphereDataComplex b = v0 + ir*opComplexExt.robert_grad_lat(phi);

					rhs = f*a + alpha*b;

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "8b", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}
#endif

				{
					/*
					 * Test F=... formulation
					 *
					 * rexi sph ver 15, sec. 5.3.7
					 * Solution for geopotential
					 */

					SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverPhi;
					SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverVel;
					double avg_geopotential = phi_bar;
					double inv_r = 1.0/ir;

					SphereDataComplex phi0 = prog_phi0_cplx_ext;
					SphereDataComplex u0 = prog_u0_cplx_ext;
					SphereDataComplex v0 = prog_v0_cplx_ext;

					SphereDataComplex phi = prog_phi_cplx_ext;
					SphereDataComplex u = prog_u_cplx_ext;
					SphereDataComplex v = prog_v_cplx_ext;

					{
						SphereDataComplex mu(sphereDataConfigExt);
						mu.physical_update_lambda_gaussian_grid(
								[&](double lon, double mu, std::complex<double> &o_data)
								{
									o_data = mu;
								}
							);


						/*
						 * lhs_direct contains the directly (forward) computed values
						 */
						SphereDataComplex lhs_direct(sphereDataConfigExt);
						lhs_direct.physical_set_zero();

						// Using sphereDataConfig results in accurate results!!!!!!!
						sphSolverPhi.setup(sphereDataConfig, 4);
						sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
						lhs_direct += (alpha*alpha)*(alpha*alpha)*phi;

						//if (use_formulation_with_coriolis_effect)
						{
							sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega*alpha*alpha, r);
							lhs_direct += (2.0*two_omega*two_omega*alpha*alpha)*opComplexExt.mu2(phi);

							sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);
							lhs_direct += (two_omega*two_omega)*(two_omega*two_omega)*opComplexExt.mu2(opComplexExt.mu2(phi));

							sphSolverPhi.solver_component_rexi_z4robert(	-avg_geopotential*alpha*two_omega, r);
							lhs_direct += (-avg_geopotential*alpha*two_omega)*(1.0/(r*r))*/* 1/cos^2phi opComplexExt.robert_grad_lat(mu)*/ opComplexExt.robert_grad_lon(phi);

							sphSolverPhi.solver_component_rexi_z5robert(	avg_geopotential/alpha*two_omega*two_omega*two_omega, r);
							lhs_direct += (avg_geopotential/alpha*two_omega*two_omega*two_omega)*(1.0/(r*r))*opComplexExt.mu2(opComplexExt.robert_grad_lon(phi));

							sphSolverPhi.solver_component_rexi_z6robert(	avg_geopotential*2.0*two_omega*two_omega, r);
							lhs_direct += (avg_geopotential*2.0*two_omega*two_omega)*(1.0/(r*r))*opComplexExt.mu(opComplexExt.robert_grad_lat(phi));
						}

						sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential*alpha*alpha, r);
						lhs_direct += (-avg_geopotential*alpha*alpha)*(1.0/(r*r))*opComplexExt.laplace(phi);

						//if (use_formulation_with_coriolis_effect)
						{
							sphSolverPhi.solver_component_rexi_z8(	-avg_geopotential*two_omega*two_omega, r);
							lhs_direct += (-avg_geopotential*two_omega*two_omega)*(1.0/(r*r))*opComplexExt.mu2(opComplexExt.laplace(phi));
						}


						SphereDataComplex rhs_direct(sphereDataConfigExt);
						{

							auto kappa = [&](
									const SphereDataComplex &i_data
							) -> SphereDataComplex
							{
								return (alpha*alpha)*i_data + two_omega*two_omega*opComplexExt.mu2(i_data);
							};

							SphereDataComplex div0 = inv_r*opComplexExt.robert_div(u0, v0);
							SphereDataComplex eta0 = inv_r*opComplexExt.robert_vort(u0, v0);

							double fj = inv_r*two_omega;
							double phi_bar = avg_geopotential;

							SphereDataComplex f(sphereDataConfigExt);
							f.physical_update_lambda_gaussian_grid(
									[&](double lon, double mu, std::complex<double> &o_data)
									{
										o_data = mu*two_omega;
									}
								);

							SphereDataComplex Fp_i = fj*(-(alpha*alpha-f*f));
							SphereDataComplex Fp_j = fj*(2.0*alpha*f);

							SphereDataComplex Fck = Fp_i*u0 + Fp_j*v0;

							rhs_direct =
									kappa(
											phi_bar*(div0 - f*(1.0/alpha)*eta0)
											+ (alpha + f*f*(1.0/alpha))*phi0
									)
									- phi_bar/alpha*Fck;
						}
//						rhs_direct = rhs_direct.spectral_returnWithTruncatedModes(sphereDataConfig);

//						sphSolverVel.setup(sphereDataConfig, 2);
//						sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
						//if (use_formulation_with_coriolis_effect)
						{
//							sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);
						}

						SphereDataComplex computed_phi_lhs_direct = sphSolverPhi.solve(lhs_direct.spectral_returnWithDifferentModes(sphereDataConfig));
						SphereDataComplex computed_phi_rhs_direct = sphSolverPhi.solve(rhs_direct.spectral_returnWithDifferentModes(sphereDataConfig));

#if 0
						SphereDataComplex lhs_direct =
								  kappa*kappa*phi
								+ phi_bar/alpha*(Fp_i*ir*opComplexExt.robert_grad_lon(phi) + Fp_j*ir*opComplexExt.robert_grad_lat(phi))
								- kappa*phi_bar*ir*ir*opComplexExt.laplace(phi);
#endif

						SphereDataComplex phi_reduced = phi.spectral_returnWithDifferentModes(sphereDataConfig);

						ErrorCheck::checkTruncated(computed_phi_lhs_direct, phi_reduced, sphereDataConfig, "test Z solvers LHS phi", epsilon, true);
						ErrorCheck::checkTruncated(computed_phi_rhs_direct, phi_reduced, sphereDataConfig, "test Z solvers RHS phi", epsilon, true);

						std::cout << "OK" << std::endl;
					}

#if 0
					auto kappa = [&](
							const SphereDataComplex &i_data
					)	-> SphereDataComplex
					{
						return (alpha*alpha)*i_data + two_omega*two_omega*opComplexExt.mu2(i_data);
					};

					{
						double fj = inv_r*two_omega;
						double phi_bar = avg_geopotential;

						SphereDataComplex f(sphereDataConfig);
						f.physical_update_lambda_gaussian_grid(
								[&](double lon, double mu, std::complex<double> &o_data)
								{
									o_data = mu*two_omega;
								}
							);

						SphereDataComplex Fp_i = fj*(-(alpha*alpha-f*f));
						SphereDataComplex Fp_j = fj*(2.0*alpha*f);

						SphereDataComplex Fck = Fp_i*u0 + Fp_j*v0;

						SphereDataComplex rhs =
								kappa(
										phi_bar*(div0 - f*(1.0/alpha)*eta0)
										+ (alpha + f*f*(1.0/alpha))*phi0
								)
								- phi_bar/alpha*Fck;

						SphereDataComplex phi = sphSolverPhi.solve(rhs);

						ErrorCheck::checkTruncated(phi, prog_u_cplx, sphereDataConfig, "VEL (phi) u_cplx", epsilon, simVars.rexi.rexi_use_extended_modes < 2);

						SphereDataComplex a = u0 + inv_r*opComplexExt.robert_grad_lon(phi);
						SphereDataComplex b = v0 + inv_r*opComplexExt.robert_grad_lat(phi);

						SphereDataComplex rhsa = alpha*a - two_omega*opComplexExt.mu(b);
						SphereDataComplex rhsb = two_omega*opComplexExt.mu(a) + alpha*b;

						u = sphSolverVel.solve(rhsa);
						v = sphSolverVel.solve(rhsb);
					}
					{
						SphereDataComplex one(sphereDataConfigExt);
						one.physical_set_zero();
						one = one+1.0;

						SphereDataComplex kappa = alpha*alpha+f*f;

						SphereDataComplex Fp_i = fj*(-(alpha*alpha-f*f));
						SphereDataComplex Fp_j = fj*(2.0*alpha*f);

						SphereDataComplex Fc = Fp_i*u0 + Fp_j*v0;

						SphereDataComplex lhs =
								  kappa*kappa*phi
								+ phi_bar/alpha*(Fp_i*ir*opComplexExt.robert_grad_lon(phi) + Fp_j*ir*opComplexExt.robert_grad_lat(phi))
								- kappa*phi_bar*ir*ir*opComplexExt.laplace(phi);

						SphereDataComplex rhs =
								kappa*phi_bar*(div0 - f*(1.0/alpha)*eta0)
								+ kappa*(alpha + f*f*(1.0/alpha))*phi0
								- phi_bar/alpha*Fc;
					}

					ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "combined a", epsilon, true);
#endif
				}
			}
#endif

#if 0
			if (false)
			{

				SphereDataComplex f(sphereDataConfig);
				f.physical_update_lambda_gaussian_grid(
						[&](double lon, double mu, std::complex<double> &o_data)
						{
							o_data = mu*two_omega;
						}
					);

				SphereDataComplex &phi0 = prog_phi0_cplx;
				SphereDataComplex &u0 = prog_u0_cplx;
				SphereDataComplex &v0 = prog_v0_cplx;

				SphereDataComplex phi = prog_phi_cplx;
				SphereDataComplex u = prog_u_cplx;
				SphereDataComplex v = prog_v_cplx;

				SphereDataComplex div0 = ir*opComplexExt.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplexExt.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplexExt.robert_div(u, v);
				SphereDataComplex eta = ir*opComplexExt.robert_vort(u, v);


				/*
				 * REXI SPH document ver 15, Section 5.3.6
				 *
				 * eq. V= A^{-1} (V0 + GRAD(PHI)) = A^{-1} T
				 */
				SphereDataComplex one(sphereDataConfig);
				one.physical_set_zero();
				one = one+1.0;

				SphereDataComplex kappa = alpha*alpha+f*f;
				SphereDataComplex inv_kappa = one/kappa;

#if 1
				SphereDataComplex Ti = u0 + ir*opComplexExt.robert_grad_lon(prog_phi_cplx);
				SphereDataComplex Tj = v0 + ir*opComplexExt.robert_grad_lat(prog_phi_cplx);
#else
				SphereDataComplex Ti = u0 + ir*opComplexExt.robert_grad_lon(Convert_SphereData_To_SphereDataComplex::physical_convert(rexi_prog_phi));
				SphereDataComplex Tj = v0 + ir*opComplexExt.robert_grad_lat(Convert_SphereData_To_SphereDataComplex::physical_convert(rexi_prog_phi));
#endif

				SphereDataComplex lhs = u;
				SphereDataComplex rhs =	inv_kappa*(
							alpha*Ti
							- two_omega*opComplexExt.mu(Tj)
						);

				ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "VEL (phi) REXI u", epsilon, simVars.rexi.rexi_use_extended_modes < 2);

				lhs = v;
				rhs =	inv_kappa*(
						  two_omega*opComplexExt.mu(Ti)
						+ alpha*Tj
						);

				ErrorCheck::checkTruncated(lhs, rhs, sphereDataConfig, "VEL (phi) REXI v", epsilon, simVars.rexi.rexi_use_extended_modes < 2);

				{
					SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverVel;

					sphSolverVel.setup(sphereDataConfig, 2);
					sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
					sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);

					double inv_r = ir;
					SphereDataComplex a = u0 + inv_r*opComplexExt.robert_grad_lon(prog_phi_cplx);
					SphereDataComplex b = v0 + inv_r*opComplexExt.robert_grad_lat(prog_phi_cplx);

					SphereDataComplex rhsa = alpha*a - two_omega*opComplexExt.mu(b);
					SphereDataComplex rhsb = two_omega*opComplexExt.mu(a) + alpha*b;

					SphereDataComplex u_cplx = sphSolverVel.solve(rhsa);
					SphereDataComplex v_cplx = sphSolverVel.solve(rhsb);


					ErrorCheck::checkTruncated(u_cplx, prog_u_cplx, sphereDataConfig, "VEL (phi) u_cplx", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
					ErrorCheck::checkTruncated(v_cplx, prog_v_cplx, sphereDataConfig, "VEL (phi) v", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}
			}
#endif

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

					rexiSPHRobert.setup(
							sphereDataConfigExt,
							sphereDataConfig,
							alpha,
							beta,
							simVars.sim.earth_radius,
							simVars.sim.coriolis_omega,
							phi_bar,
							timestep_size,
							param_use_coriolis_formulation
					);

					if (use_complex_valued_solver)
					{
#if 1
						rexiSPHRobert.solve_complexRHS(
								prog_phi0_cplx_ext,
								prog_u0_cplx_ext,
								prog_v0_cplx_ext,

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

						rexiSPHRobert.solve(
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
			else
			{

				{
					SWERexiTerm_SPH rexiSPH;

					rexiSPH.setup(
							sphereDataConfigExt,
							alpha,
							beta,
							simVars.sim.earth_radius,
							simVars.sim.coriolis_omega,
							phi_bar,
							timestep_size,
							param_use_coriolis_formulation
					);

					if (use_complex_valued_solver)
					{
						rexiSPH.solve_complexRHS(
								prog_phi0_cplx_ext,
								prog_u0_cplx_ext,
								prog_v0_cplx_ext,

								rexi_prog_phi_ext,
								rexi_prog_u_ext,
								rexi_prog_v_ext
							);
					}
					else
					{
						rexiSPH.solve(
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
			SphereData prog_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_phi_cplx);
			SphereData prog_u = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_u_cplx);
			SphereData prog_v = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_v_cplx);


			/*
			 * REXI results stored in rexi_prog_*
			 * These should match prog_*
			 */
			ErrorCheck::checkTruncated(rexi_prog_phi, prog_phi, sphereDataConfig, "prog_phi", epsilon*1e+5, false, simVars.setup.benchmark_scenario_id != 1);
			ErrorCheck::checkTruncated(rexi_prog_u, prog_u, sphereDataConfig, "prog_u", epsilon*1e+5, false, simVars.setup.benchmark_scenario_id != 1);
			ErrorCheck::checkTruncated(rexi_prog_v, prog_v, sphereDataConfig, "prog_v", epsilon*1e+5, false, simVars.setup.benchmark_scenario_id != 1);
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
			"rexi-use-coriolis-formulation",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 1;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << "	--rexi-use-coriolis-formulation [0/1]	Use REXI formulation with coriolis effect" << std::endl;

#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}

	param_use_coriolis_formulation = simVars.bogus.var[0];
	assert (param_use_coriolis_formulation == 0 || param_use_coriolis_formulation == 1);

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
