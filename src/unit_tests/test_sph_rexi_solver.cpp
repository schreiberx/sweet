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
#include <rexi/SWERexi_SPHRobert.hpp>
#include <rexi/SWERexi_SPH.hpp>
#include <rexi/REXI.hpp>
#include <../programs/swe_sphere_rexi/SWE_Sphere_REXI.hpp>

#include <sweet/sphere/GenerateConsistentGradDivSphereData.hpp>


SimulationVariables simVars;
SphereOperatorsComplex opComplex;

SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

SphereDataConfig sphereDataConfigRexiAddedModes;
SphereDataConfig *sphereDataConfigRexi = &sphereDataConfigRexiAddedModes;


bool param_rexi_use_coriolis_formulation = true;


/**
 * Run with
 *
 * 	$ ./build/sh_example T32 P2
 */
void run_tests()
{

	std::cout << "Using time step size dt = " << simVars.timecontrol.current_timestep_size << std::endl;
	std::cout << "Running simulation until t_end = " << simVars.timecontrol.max_simulation_time << std::endl;
	std::cout << "Parameters:" << std::endl;
	std::cout << " + Gravity: " << simVars.sim.gravitation << std::endl;
	std::cout << " + Earth_radius: " << simVars.sim.earth_radius << std::endl;
	std::cout << " + Average height: " << simVars.sim.h0 << std::endl;
	std::cout << " + Coriolis_omega: " << simVars.sim.coriolis_omega << std::endl;
	std::cout << " + Viscosity D: " << simVars.sim.viscosity << std::endl;
	std::cout << " + use_nonlinear: " << simVars.misc.use_nonlinear_equations << std::endl;
	std::cout << " + Use REXI Coriolis formulation: " << (param_rexi_use_coriolis_formulation ? "true" : "false") << std::endl;
	std::cout << std::endl;
	std::cout << " + Benchmark scenario id: " << simVars.setup.benchmark_scenario_id << std::endl;
	std::cout << " + Use robert functions: " << simVars.misc.sphere_use_robert_functions << std::endl;
	std::cout << " + Use REXI: " << simVars.rexi.use_rexi << std::endl;
	std::cout << " + REXI h: " << simVars.rexi.rexi_h << std::endl;
	std::cout << " + REXI M: " << simVars.rexi.rexi_M << std::endl;
	std::cout << " + REXI use half poles: " << simVars.rexi.rexi_use_half_poles << std::endl;
	std::cout << " + REXI additional modes: " << simVars.rexi.rexi_use_extended_modes << std::endl;
	std::cout << std::endl;
	std::cout << " + RK order: " << simVars.disc.timestepping_runge_kutta_order << std::endl;
	std::cout << " + timestep size: " << simVars.timecontrol.current_timestep_size << std::endl;
	std::cout << " + output timestep size: " << simVars.misc.output_each_sim_seconds << std::endl;


	double epsilon = 1e-10;
	epsilon *= (sphereDataConfig->spectral_modes_n_max);

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
		sphereDataConfigRexi = sphereDataConfig;
	}
	else
	{
		// Add modes only along latitude since these are the "problematic" modes
		sphereDataConfigRexiAddedModes.setupAdditionalModes(
				sphereDataConfig,
				simVars.rexi.rexi_use_extended_modes,	// TODO: Extend SPH wrapper to also support m != n to set this guy to 0
				simVars.rexi.rexi_use_extended_modes
		);
		sphereDataConfigRexi = &sphereDataConfigRexiAddedModes;
	}

	REXI rexi;
	rexi.setup(simVars.rexi.rexi_h, simVars.rexi.rexi_M);

	if (!simVars.misc.sphere_use_robert_functions)
	{
		std::cerr << "WARNING: TESTS WILL FAIL WITH NON-ROBERT FUNCTIONS DUE TO SINGULARITY AT POLES!" << std::endl;
	}

	double timestep_size = 1.0;

	std::cout << "Using generic max allowed error value of " << epsilon << " (problem specific scaled)" << std::endl;

	std::complex<double> beta(1.0, 0.0);

	for (int use_complex_valued_solver = 0; use_complex_valued_solver < 2; use_complex_valued_solver++)
	{
		for (std::size_t i = 0; i < rexi.alpha.size(); i++)
		{
			std::complex<double> &alpha = rexi.alpha[i];

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


			{
				SphereOperators op;

				GenerateConsistentGradDivSphereData g_real(
						simVars,
						sphereDataConfig,
						op,
						M_PI/3, M_PI/3,
						std::string("gen_divgrad_data_re_")+sphereDataConfig->getUniqueIDString()
				);

				g_real.generate();

				prog_phi_cplx = Convert_SphereData_To_SphereDataComplex::physical_convert(g_real.prog_h*simVars.sim.gravitation);
				prog_u_cplx = Convert_SphereData_To_SphereDataComplex::physical_convert(g_real.prog_u);
				prog_v_cplx = Convert_SphereData_To_SphereDataComplex::physical_convert(g_real.prog_v);


				GenerateConsistentGradDivSphereData g_imag(
						simVars,
						sphereDataConfig,
						op,
						M_PI/4, M_PI/2,
						std::string("gen_divgrad_data_im_")+sphereDataConfig->getUniqueIDString()
				);

				g_imag.generate();

				prog_phi_cplx = prog_phi_cplx + Convert_SphereData_To_SphereDataComplex::physical_convert(g_imag.prog_h*simVars.sim.gravitation)*std::complex<double>(0,1);
				prog_u_cplx = prog_u_cplx + Convert_SphereData_To_SphereDataComplex::physical_convert(g_imag.prog_u)*std::complex<double>(0,1);
				prog_v_cplx = prog_v_cplx + Convert_SphereData_To_SphereDataComplex::physical_convert(g_imag.prog_v)*std::complex<double>(0,1);


				prog_phi_cplx.spectral_truncate();
				prog_u_cplx.spectral_truncate();
				prog_v_cplx.spectral_truncate();
			}

			SphereDataComplex prog_phi_cplx_ext = prog_phi_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
			SphereDataComplex prog_u_cplx_ext = prog_u_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
			SphereDataComplex prog_v_cplx_ext = prog_v_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);


			/*
			 * Computed right hand side
			 */
			double phi_bar = simVars.sim.gravitation*simVars.sim.h0;

			SphereDataComplex prog_phi0_cplx_ext(sphereDataConfigRexi);
			SphereDataComplex prog_u0_cplx_ext(sphereDataConfigRexi);
			SphereDataComplex prog_v0_cplx_ext(sphereDataConfigRexi);


			/*
			 * Compute RHS for REXI term
			 */
			if (simVars.misc.sphere_use_robert_functions)
			{
				prog_phi0_cplx_ext =
							+ alpha*prog_phi_cplx_ext
							- phi_bar*ir*opComplex.robert_div_lon(prog_u_cplx_ext)
							- phi_bar*ir*opComplex.robert_div_lat(prog_v_cplx_ext);

				prog_u0_cplx_ext =
							- ir*opComplex.robert_grad_lon(prog_phi_cplx_ext)
							+ alpha*prog_u_cplx_ext
							+ two_omega*opComplex.mu(prog_v_cplx_ext);

				prog_v0_cplx_ext =
							- ir*opComplex.robert_grad_lat(prog_phi_cplx_ext)
							- two_omega*opComplex.mu(prog_u_cplx_ext)
							+ alpha*prog_v_cplx_ext;

			}
			else
			{
				prog_phi0_cplx_ext =
							+ alpha*prog_phi_cplx_ext
							- phi_bar*ir*opComplex.div_lon(prog_u_cplx_ext)
							- phi_bar*ir*opComplex.div_lat(prog_v_cplx_ext);

				prog_u0_cplx_ext =
							- ir*opComplex.grad_lon(prog_phi_cplx_ext)
							+ alpha*prog_u_cplx_ext
							+ two_omega*opComplex.mu(prog_v_cplx_ext);

				prog_v0_cplx_ext =
							- ir*opComplex.grad_lat(prog_phi_cplx_ext)
							- two_omega*opComplex.mu(prog_u_cplx_ext)
							+ alpha*prog_v_cplx_ext;
			}



			SphereData prog_phi0_ext(sphereDataConfigRexi);
			SphereData prog_u0_ext(sphereDataConfigRexi);
			SphereData prog_v0_ext(sphereDataConfigRexi);

			{
				/*
				 * Test correct mode extensions for sphere data with complex physical data
				 */
				SphereDataComplex lhs = prog_u0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);
				SphereDataComplex rhs = prog_u0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfigRexi).spectral_returnWithDifferentModes(sphereDataConfig);

				SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
				SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

				std::cout << "+ error for mode extensions: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() << std::endl;

				if ((lhsr-rhsr).physical_reduce_error_max_abs() > 10e-10)
					FatalError("Threshold exceeded");
			}

			if (!use_complex_valued_solver)
			{
				prog_phi0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_phi0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigRexi);
				prog_u0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_u0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigRexi);
				prog_v0_ext = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_v0_cplx_ext).spectral_returnWithDifferentModes(sphereDataConfigRexi);
			}


#if 1
			if (!simVars.misc.sphere_use_robert_functions)
			{
				SphereDataComplex phi(sphereDataConfigRexi);
				SphereDataComplex u(sphereDataConfigRexi);
				SphereDataComplex v(sphereDataConfigRexi);


				phi = prog_phi_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
				u = prog_u_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
				v = prog_v_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);



//				SphereDataComplex &phi0 = prog_phi0_cplx_ext;
				SphereDataComplex &u0 = prog_u0_cplx_ext;
				SphereDataComplex &v0 = prog_v0_cplx_ext;


				SphereDataComplex f(sphereDataConfigRexi);
				f.physical_update_lambda_gaussian_grid(
						[&](double lon, double mu, std::complex<double> &o_data)
						{
							o_data = mu*two_omega;
						}
					);

				SphereDataComplex div0 = ir*opComplex.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplex.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplex.robert_div(u, v);
				SphereDataComplex eta = ir*opComplex.robert_vort(u, v);

				SphereDataComplex fi = ir*opComplex.robert_grad_lon(f);
				SphereDataComplex fj = ir*opComplex.robert_grad_lat(f);

				SphereDataComplex lhs(sphereDataConfigRexi);
				SphereDataComplex rhs(sphereDataConfigRexi);

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Check Identity relation
					 *
					 * D.(fA) = f(D.A) + A.(Gf)
					 */

					lhs = ir*opComplex.div(f*u, f*v);

					rhs =	  f*ir*opComplex.div(u, v)
							+ u*ir*opComplex.grad_lon(f)
							+ v*ir*opComplex.grad_lat(f);

					lhs.physical_truncate();
					rhs.physical_truncate();

					Convert_SphereDataComplex_To_SphereData::physical_convert(lhs).physical_file_write("test_lhs.csv");
					Convert_SphereDataComplex_To_SphereData::physical_convert(rhs).physical_file_write("test_rhs.csv");



					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR Identity for non-Robert formulation: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon)
						FatalError("Error too large");
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

				SphereDataComplex phi(sphereDataConfigRexi);
				SphereDataComplex u(sphereDataConfigRexi);
				SphereDataComplex v(sphereDataConfigRexi);


				phi = prog_phi_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
				u = prog_u_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
				v = prog_v_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);

				//
				SphereDataComplex f(sphereDataConfigRexi);
				f.physical_update_lambda_gaussian_grid(
						[&](double lon, double mu, std::complex<double> &o_data)
						{
							o_data = mu*two_omega;
						}
					);

				SphereDataComplex div0 = ir*opComplex.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplex.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplex.robert_div(u, v);
				SphereDataComplex eta = ir*opComplex.robert_vort(u, v);


				/*
				 * M(\mu) term for
				 * 	D.(f A) = f(D.A) + M(\my)*A.grad(f)
				 * fix
				 */
				SphereDataComplex one_over_cos2phi(sphereDataConfigRexi);
				one_over_cos2phi.physical_update_lambda_cosphi_grid(
						[&](double lon, double cosphi, std::complex<double> &o_data)
						{
							o_data = 1.0/(cosphi*cosphi);
						}
				);
#if 0
				SphereDataComplex fi(sphereDataConfigRexi);
				fi.physical_set_zero();

				SphereDataComplex fj(sphereDataConfigRexi);
	#if 0

				fj.physical_update_lambda_cosphi_grid(
						[&](double lon, double cosphi, std::complex<double> &o_data)
						{
							o_data = ir*two_omega*cosphi*cosphi;
						}
					);

				fj = fj*one_over_cos2phi;
	#else
				fj.physical_update_lambda_cosphi_grid(
						[&](double lon, double cosphi, std::complex<double> &o_data)
						{
							o_data = ir*two_omega;
						}
					);
	#endif

#else
				double fi = 0;
				double fj = ir*two_omega;
#endif

				SphereDataComplex lhs(sphereDataConfigRexi);
				SphereDataComplex rhs(sphereDataConfigRexi);


				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Check Identity relation
					 *
					 * D.(fA) = f(D.A) + A.(Df)
					 */

					lhs = opComplex.robert_div(f*u, f*v);

					rhs =	  f*opComplex.robert_div(u, v)
								+ one_over_cos2phi*
								(
									u*opComplex.robert_grad_lon(f)
									+ v*opComplex.robert_grad_lat(f)
								);

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR IDENTITY D.(fA): " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						if (simVars.rexi.rexi_use_extended_modes >= 2)
						{
							Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
							Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
							Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
							FatalError("Error too large");
						}
						else
						{
							std::cerr << "Error ignored since REXI extended modes < 2" << std::endl;
						}
					}
				}


				{
					lhs = opComplex.laplace(phi);

					rhs =	  opComplex.robert_div_lat(opComplex.robert_grad_lat(phi))
							+ opComplex.robert_div_lon(opComplex.robert_grad_lon(phi))
						;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR laplace: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}


				if (simVars.rexi.rexi_use_extended_modes == 0)
				{
					lhs = opComplex.robert_div_lat(v);

					rhs = (
							opComplex.robert_div_lat(prog_v_cplx)
						).spectral_returnWithDifferentModes(sphereDataConfigRexi);

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 0:\t" << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2, eq (2)
					 */

					lhs = -ir*opComplex.robert_grad_lon(phi) + alpha*u + f*v;
					rhs = u0;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 0b: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}
				{
					/*
					 * REXI SPH document ver 14, Section 5.2, eq (3)
					 */

					lhs = -ir*opComplex.robert_grad_lat(phi) - f*u + alpha*v;
					rhs = v0;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 0c: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}
				{
					lhs = div;
					rhs = ir*opComplex.robert_div_lon(u) + ir*opComplex.robert_div_lat(v);

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 0d: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}
				{
					lhs = eta;
					rhs = ir*opComplex.robert_div_lon(v) - ir*opComplex.robert_div_lat(u);

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 0e: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Check Identity relation
					 *
					 * D.(fA) = f(D.A) + A.one_over_cos2phi*(Df)
					 */

					lhs = ir*opComplex.robert_div(f*u, f*v);

					rhs =	f*ir*opComplex.robert_div(u, v)
							+ u*fi + v*fj;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 0f: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						if (simVars.rexi.rexi_use_extended_modes >= 2)
						{
							Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
							Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
							Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
							FatalError("Error too large");
						}
						else
						{
							std::cerr << "Error ignored since REXI extended modes < 2" << std::endl;
						}
					}
				}



				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st/1st equation in this section
					 */
					lhs = ir*opComplex.robert_div_lon(
								-ir*opComplex.robert_grad_lon(phi) + alpha*u + two_omega*opComplex.mu(v)
							  );

					rhs = ir*opComplex.robert_div_lon(u0);


					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 1aa: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+3)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st/2nd equation in this section
					 */

					lhs = ir*opComplex.robert_div_lat(v0);

					rhs = ir*opComplex.robert_div_lat(
							-ir*opComplex.robert_grad_lat(phi) - f*u + alpha*v
						  );

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 1ab: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st equation in this section
					 */

					lhs = div0;

					rhs =	ir*opComplex.robert_div_lon(
								-ir*opComplex.robert_grad_lon(phi) + alpha*u + f*v
							  )
							+ ir*opComplex.robert_div_lat(
								-ir*opComplex.robert_grad_lat(phi) + alpha*v - f*u
							  );

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 1a: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st equation in this section
					 */

					lhs = div0;

					rhs =	  ir*opComplex.robert_div_lon(	-ir*opComplex.robert_grad_lon(phi)	)
							+ ir*opComplex.robert_div_lon(	alpha*u 	)
							+ ir*opComplex.robert_div_lon(	f*v		)

							+ ir*opComplex.robert_div_lat(	-ir*opComplex.robert_grad_lat(phi)	)
							+ ir*opComplex.robert_div_lat(	+ alpha*v 	)
							+ ir*opComplex.robert_div_lat(	-f*u	)
							;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 1ba: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st equation in this section
					 */

					lhs = div0;

					rhs =	 -ir*ir*opComplex.laplace(phi)

							+ alpha*ir*opComplex.robert_div(u, v)

							+ ir*opComplex.robert_div(	f*v	, -f*u	)
							;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 1bb: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Test identity of D.(fA) with A=(v,-u) formulation
					 */

					lhs = ir*opComplex.robert_div(	f*v	, -f*u	);

					rhs =	  f*ir*opComplex.robert_div(v, -u)
								+ one_over_cos2phi*
								(
									v*ir*opComplex.robert_grad_lon(f)
									- u*ir*opComplex.robert_grad_lat(f)
								);


					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 1ca: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");

						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * eq. (4)
					 */

					lhs =	- ir*ir*opComplex.laplace(phi)
							+ alpha*div

							+f*eta
							+ v*fi
							- u*fj
					;

					rhs = div0;


					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 1c: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.2
					 * eq. (5)
					 */

					lhs = -alpha*eta + f*div + fj*v + fi*u;
					rhs = -eta0;


					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 2: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+3)
					{
						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}



				{
					/*
					 * REXI SPH document ver 14, Section 5.2.3
					 * eq. (6)
					 */

					lhs = div;
					rhs = 1.0/phi_bar*(phi*alpha-phi0);


					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 3: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+3)
						FatalError("Error too large");
				}


				{
					/*
					 * REXI SPH document ver 14, Section 5.2.4
					 * eq. (7)
					 */

					lhs =	- ir*ir*opComplex.laplace(phi)
							+ 1.0/phi_bar * phi * alpha*alpha
							+ f * eta
							+ fi*v
							- fj*u;

					rhs =	div0 + 1.0/phi_bar * alpha * phi0;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 4a: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+3)
					{

						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}


				{
					/*
					 * REXI SPH document ver 14, Section 5.2.4
					 * between eq. (7) and (8)
					 */

					lhs = -alpha*eta + 1.0/phi_bar * f * phi * alpha + fj*v + fi*u;
					rhs = -eta0 + 1.0/phi_bar * f * phi0;


					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);


					std::cout << "ERROR 4ba: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+1)
					{

						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}


				{
					/*
					 * REXI SPH document ver 14, Section 5.2.4
					 * eq. (8)
					 */

					lhs = eta;
					rhs =	1.0/alpha * eta0
							- 1.0/(phi_bar*alpha)*f*phi0
							+ 1.0/phi_bar*f*phi
							+ 1.0/alpha * fj*v
							+ 1.0/alpha * fi*u;

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 5: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}



				{
					/*
					 * REXI SPH document ver 14, Section 5.2.5 at the end
					 * eq. (9)
					 */
#if 0
					SphereDataComplex F = f*(fj*v + fi*u) + alpha*(fi*v - fj*u);

					lhs = alpha*alpha*phi + f*f*phi - phi_bar*ir*ir*opComplex.laplace(phi) + (phi_bar/alpha) * F;
					rhs = phi_bar*(div0 - f*(1.0/alpha)*eta0) + (alpha+f*f*(1.0/alpha))*phi0;

#else

					SphereDataComplex F = two_omega*opComplex.mu(fj*v + fi*u) + alpha*(fi*v - fj*u);


					lhs = alpha*alpha*phi + two_omega*two_omega*opComplex.mu2(phi) - phi_bar*ir*ir*opComplex.laplace(phi) + (phi_bar/alpha) * F;
					rhs = phi_bar*(div0 - two_omega*opComplex.mu((1.0/alpha)*eta0)) + (alpha*phi0+(1.0/alpha)*two_omega*two_omega*opComplex.mu2(phi0));

#endif

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 5_blarg: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}
#if 0
				{
					/*
					 * REXI SPH document ver 15, Section 8 at the end
					 */
					lhs = (
							(alpha*alpha+f*f)*phi
							+phi_var/alpha*Fp*opComplex.grad/*nonsense*/(phi)
							+phi_bar*opComplex.laplace(phi)
						);

					rhs = phi_bar*(div0 - two_omega*opComplex.mu((1.0/alpha)*eta0)) + (alpha*phi0+(1.0/alpha)*two_omega*two_omega*opComplex.mu2(phi0));

					SphereDataComplex lhsr = lhs.spectral_returnWithDifferentModes(sphereDataConfig);
					SphereDataComplex rhsr = rhs.spectral_returnWithDifferentModes(sphereDataConfig);

					std::cout << "ERROR 8: " << (lhsr-rhsr).physical_reduce_error_max_abs() << "\t" << (lhsr-rhsr).physical_reduce_error_rms() <<std::endl;

					if ((lhsr-rhsr).physical_reduce_error_max_abs() > epsilon*10e+4)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");
						if (simVars.rexi.rexi_use_extended_modes < 2)
							std::cerr << "Ignoring error since extended modes < 2" << std::endl;
						else
							FatalError("Error too large");
					}
				}
#endif

			}
#endif
			////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////
			continue;
			////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////

			/*
			 * Output of REXI (ext modes)
			 */
			SphereData rexi_prog_phi_ext(sphereDataConfigRexi);
			SphereData rexi_prog_u_ext(sphereDataConfigRexi);
			SphereData rexi_prog_v_ext(sphereDataConfigRexi);

			if (simVars.misc.sphere_use_robert_functions)
			{
				{
					SWERexi_SPHRobert rexiSPHRobert;

					rexiSPHRobert.setup(
							sphereDataConfigRexi,
							alpha,
							beta,
							simVars.sim.earth_radius,
							simVars.sim.coriolis_omega,
							phi_bar,
							timestep_size,
							param_rexi_use_coriolis_formulation
					);

					if (use_complex_valued_solver)
					{
						rexiSPHRobert.solve_complexRHS(
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
					SWERexi_SPH rexiSPH;

					rexiSPH.setup(
							sphereDataConfigRexi,
							alpha,
							beta,
							simVars.sim.earth_radius,
							simVars.sim.coriolis_omega,
							phi_bar,
							timestep_size,
							param_rexi_use_coriolis_formulation
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
			SphereData rexi_prog_phi(sphereDataConfig);
			SphereData rexi_prog_u(sphereDataConfig);
			SphereData rexi_prog_v(sphereDataConfig);


			rexi_prog_phi = rexi_prog_phi_ext.spectral_returnWithDifferentModes(rexi_prog_phi.sphereDataConfig);
			rexi_prog_u = rexi_prog_u_ext.spectral_returnWithDifferentModes(rexi_prog_phi.sphereDataConfig);
			rexi_prog_v = rexi_prog_v_ext.spectral_returnWithDifferentModes(rexi_prog_phi.sphereDataConfig);


			SphereData prog_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_phi_cplx);
			SphereData prog_u = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_u_cplx);
			SphereData prog_v = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_v_cplx);


#if 1
			{
				std::ostringstream obuf;
				obuf << "prog_alphaid" << i << "_phi.csv";
				prog_phi.physical_file_write_lon_pi_shifted(obuf.str().c_str());
			}
			{
				std::ostringstream obuf;
				obuf << "prog_alphaid" << i << "_phi_rexi.csv";
				rexi_prog_phi.physical_file_write_lon_pi_shifted(obuf.str().c_str());
			}
			{
				std::ostringstream obuf;
				obuf << "prog_alphaid" << i << "_diff.csv";
				(rexi_prog_phi-prog_phi).physical_file_write_lon_pi_shifted(obuf.str().c_str());
			}
#endif
			/*
			 * REXI results stored in rexi_prog_*
			 * These should match prog_*
			 */
			double err_phi = (rexi_prog_phi - prog_phi).physical_reduce_max_abs();
			double err_u = (rexi_prog_u - prog_u).physical_reduce_max_abs();
			double err_v = (rexi_prog_v - prog_v).physical_reduce_max_abs();

			std::cout << i << ":\t";
			std::cout << alpha << "\t";
			std::cout << err_phi << "\t";
			std::cout << err_u << "\t";
			std::cout << err_v << std::endl;

			if (err_phi > epsilon || err_u > epsilon || err_v > epsilon)
				FatalError("Error too large");


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
	simVars.bogus.var[0] = 0;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}

	param_rexi_use_coriolis_formulation = simVars.bogus.var[0];
	assert (param_rexi_use_coriolis_formulation == 0 || param_rexi_use_coriolis_formulation == 1);

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
