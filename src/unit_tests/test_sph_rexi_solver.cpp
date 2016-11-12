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
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigRexiAddedModes;


bool param_rexi_use_coriolis_formulation = true;


void errorCheck(
		const SphereDataComplex &i_lhs,
		const SphereDataComplex &i_rhs,
		const std::string &i_id,
		double i_error_threshold = 1.0,
		double i_ignore_error = false
)
{
	SphereDataComplex lhsr = i_lhs.spectral_returnWithDifferentModes(sphereDataConfig);
	SphereDataComplex rhsr = i_rhs.spectral_returnWithDifferentModes(sphereDataConfig);

	double normalize_fac = std::min(lhsr.physical_reduce_max_abs(), rhsr.physical_reduce_max_abs());

	SphereDataComplex diff = lhsr-rhsr;
	diff.physical_reduce_max_abs();

	if (normalize_fac == 0)
	{
		std::cout << "Error computation for '" << i_id << "' ignored since both fields are Zero" << std::endl;
		return;
	}
	double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
	double rel_rms = diff.physical_reduce_rms() / normalize_fac;

	std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tThreshold: " << i_error_threshold << std::endl;

	if (rel_max_abs > i_error_threshold)
	{
		Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr).physical_file_write("o_error_lhs.csv");
		Convert_SphereDataComplex_To_SphereData::physical_convert(rhsr).physical_file_write("o_error_rhs.csv");
		Convert_SphereDataComplex_To_SphereData::physical_convert(lhsr-rhsr).physical_file_write("o_error_diff.csv");

		if (i_ignore_error)
			std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
		else
			FatalError("Error too large");
	}
}


void errorCheck(
		const SphereData &i_lhs,
		const SphereData &i_rhs,
		const std::string &i_id,
		double i_error_threshold = 1.0,
		double i_ignore_error = false
)
{
	SphereData lhsr = i_lhs.spectral_returnWithDifferentModes(sphereDataConfig);
	SphereData rhsr = i_rhs.spectral_returnWithDifferentModes(sphereDataConfig);

	double normalize_fac = std::min(lhsr.physical_reduce_max_abs(), rhsr.physical_reduce_max_abs());

	SphereData diff = lhsr-rhsr;
	diff.physical_reduce_max_abs();

	if (normalize_fac == 0)
	{
		std::cout << "Error computation for '" << i_id << "' ignored since both fields are Zero" << std::endl;
		return;
	}
	double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
	double rel_rms = diff.physical_reduce_rms() / normalize_fac;

	std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tThreshold: " << i_error_threshold << std::endl;

	if (rel_max_abs > i_error_threshold)
	{
		lhsr.physical_file_write("o_error_lhs.csv");
		rhsr.physical_file_write("o_error_rhs.csv");
		(lhsr-rhsr).physical_file_write("o_error_diff.csv");

		if (i_ignore_error)
			std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
		else
			FatalError("Error too large");
	}
}




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
	std::cout << std::endl;
	std::cout << " + Benchmark scenario id: " << simVars.setup.benchmark_scenario_id << std::endl;
	std::cout << " + Use robert functions: " << simVars.misc.sphere_use_robert_functions << std::endl;
	std::cout << " + Use REXI: " << simVars.rexi.use_rexi << std::endl;
	std::cout << " + REXI h: " << simVars.rexi.rexi_h << std::endl;
	std::cout << " + REXI M: " << simVars.rexi.rexi_M << std::endl;
	std::cout << " + REXI use half poles: " << simVars.rexi.rexi_use_half_poles << std::endl;
	std::cout << " + REXI additional modes: " << simVars.rexi.rexi_use_extended_modes << std::endl;
	std::cout << " + Use REXI Coriolis formulation: " << (param_rexi_use_coriolis_formulation ? "true" : "false") << std::endl;
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

	REXI rexi(0, simVars.rexi.rexi_h, simVars.rexi.rexi_M);

	if (!simVars.misc.sphere_use_robert_functions)
	{
		std::cerr << "WARNING: TESTS WILL FAIL WITH NON-ROBERT FUNCTIONS DUE TO SINGULARITY AT POLES!" << std::endl;
	}

	double timestep_size = 1.0;

	std::cout << "Using generic max allowed error value of " << epsilon << " (problem specific scaled)" << std::endl;

	std::complex<double> beta(1.0, 0.0);

	for (int use_complex_valued_solver = 1; use_complex_valued_solver < 2; use_complex_valued_solver++)
	{
		for (int i = rexi.alpha.size()-1; i >= 0; i--)
//		for (int i = 0; i < rexi.alpha.size(); i++)
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

//			prog_phi0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfig);
//			prog_u0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfig);
//			prog_v0_cplx_ext.spectral_returnWithTruncatedModes(sphereDataConfig);

			SphereDataComplex prog_phi0_cplx = prog_phi0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereDataComplex prog_u0_cplx = prog_u0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);
			SphereDataComplex prog_v0_cplx = prog_v0_cplx_ext.spectral_returnWithDifferentModes(sphereDataConfig);


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


				SphereDataComplex div0 = ir*opComplex.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplex.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplex.robert_div(u, v);
				SphereDataComplex eta = ir*opComplex.robert_vort(u, v);

				SphereDataComplex fi = ir*opComplex.robert_grad_lon(f);
				SphereDataComplex fj = ir*opComplex.robert_grad_lat(f);

				SphereDataComplex lhs(sphereDataConfigExt);
				SphereDataComplex rhs(sphereDataConfigExt);

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

					errorCheck(lhs, rhs, "ERROR Identity for non-Robert formulation", epsilon);
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


				SphereDataComplex div0 = ir*opComplex.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplex.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplex.robert_div(u, v);
				SphereDataComplex eta = ir*opComplex.robert_vort(u, v);


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

					errorCheck(lhs, rhs, "ERROR IDENTITY D.(fA)", epsilon, simVars.rexi.rexi_use_extended_modes <= 2);
				}

				{
					lhs = opComplex.laplace(phi);

					rhs =	  opComplex.robert_div_lat(opComplex.robert_grad_lat(phi))
							+ opComplex.robert_div_lon(opComplex.robert_grad_lon(phi))
						;

					errorCheck(lhs, rhs, "laplace", epsilon, simVars.rexi.rexi_use_extended_modes <= 2);
				}


				if (simVars.rexi.rexi_use_extended_modes == 0)
				{
					lhs = opComplex.robert_div_lat(v.spectral_returnWithDifferentModes(sphereDataConfig)).spectral_returnWithDifferentModes(sphereDataConfigExt);

					rhs = (
							opComplex.robert_div_lat(prog_v_cplx)
						).spectral_returnWithDifferentModes(sphereDataConfigExt);

					errorCheck(lhs, rhs, "0", epsilon, simVars.rexi.rexi_use_extended_modes <= 2);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2, eq (2)
					 */

					lhs = -ir*opComplex.robert_grad_lon(phi) + alpha*u + f*v;
					rhs = u0;

					errorCheck(lhs, rhs, "0b", epsilon);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2, eq (3)
					 */

					lhs = -ir*opComplex.robert_grad_lat(phi) - f*u + alpha*v;
					rhs = v0;

					errorCheck(lhs, rhs, "0c", epsilon);
				}

				{
					lhs = div;
					rhs = ir*opComplex.robert_div_lon(u) + ir*opComplex.robert_div_lat(v);

					errorCheck(lhs, rhs, "0d", epsilon);
				}

				{
					lhs = eta;
					rhs = ir*opComplex.robert_div_lon(v) - ir*opComplex.robert_div_lat(u);

					errorCheck(lhs, rhs, "0e", epsilon);
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

					errorCheck(lhs, rhs, "0f", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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

					errorCheck(lhs, rhs, "1aa", epsilon);
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

					errorCheck(lhs, rhs, "1ab", epsilon);
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

					errorCheck(lhs, rhs, "1a", epsilon);
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

					errorCheck(lhs, rhs, "1ba", epsilon);
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

					errorCheck(lhs, rhs, "1bb", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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

					errorCheck(lhs, rhs, "1ca", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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

					errorCheck(lhs, rhs, "1c", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.2
					 * eq. (5)
					 */

					lhs = -alpha*eta + f*div + fj*v + fi*u;
					rhs = -eta0;

					errorCheck(lhs, rhs, "2", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
				}



				{
					/*
					 * REXI SPH document ver 14, Section 5.2.3
					 * eq. (6)
					 */

					lhs = div;
					rhs = 1.0/phi_bar*(phi*alpha-phi0);

					errorCheck(lhs, rhs, "3", epsilon*10e+2);
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

					errorCheck(lhs, rhs, "4a", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * REXI SPH document ver 14, Section 5.2.4
					 * between eq. (7) and (8)
					 */

					lhs = -alpha*eta + 1.0/phi_bar * f * phi * alpha + fj*v + fi*u;
					rhs = -eta0 + 1.0/phi_bar * f * phi0;

					errorCheck(lhs, rhs, "4ba", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
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

					errorCheck(lhs, rhs, "5", epsilon*10e+2, simVars.rexi.rexi_use_extended_modes < 2);
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
					lhs =	- ir*ir*opComplex.laplace(phi)
							+ 1.0/phi_bar * phi * alpha*alpha
							+ f * ceta	/* we use the computed eta here */
							+ fi*v
							- fj*u;

					rhs =	div0 + 1.0/phi_bar * alpha * phi0;

					errorCheck(lhs, rhs, "5b", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}



				{
					/*
					 * REXI SPH document ver 15, Section 5.3.5
					 */

					lhs =
							- ir*ir*opComplex.laplace(phi)

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

					errorCheck(lhs, rhs, "5c", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}


				{
					/*
					 * REXI SPH document ver 15, Section 5.3.5 at the end
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

					errorCheck(lhs, rhs, "5d", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					/*
					 * REXI SPH document ver 15, Section 5.3.6
					 *
					 * eq. V= A^{-1} (V0 + GRAD(PHI)) = A^{-1} T
					 */

					SphereDataComplex Ti = u0 + ir*opComplex.robert_grad_lon(phi);
					SphereDataComplex Tj = v0 + ir*opComplex.robert_grad_lat(phi);

					lhs =	  alpha*alpha*u
							+ two_omega*two_omega*opComplex.mu2(u);

					rhs =	  alpha*Ti
							- two_omega*opComplex.mu(Tj);

					errorCheck(lhs, rhs, "6 lemma u", epsilon, simVars.rexi.rexi_use_extended_modes < 2);


					lhs =	  alpha*alpha*v
							+ two_omega*two_omega*opComplex.mu2(v);

					rhs =	  two_omega*opComplex.mu(Ti)
							+ alpha*Tj;

					errorCheck(lhs, rhs, "6 lemma v", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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

					SphereDataComplex Ti = u0 + ir*opComplex.robert_grad_lon(phi);
					SphereDataComplex Tj = v0 + ir*opComplex.robert_grad_lat(phi);

					lhs = u;

					rhs =	inv_kappa*(
								alpha*Ti
								- two_omega*opComplex.mu(Tj)
							);

					errorCheck(lhs, rhs, "6 lemma inv_kappa u", epsilon, simVars.rexi.rexi_use_extended_modes < 2);


					lhs = v;

					rhs =	inv_kappa*(
							  two_omega*opComplex.mu(Ti)
							+ alpha*Tj
							);

					errorCheck(lhs, rhs, "6 lemma inv_kappa v", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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

					SphereDataComplex Ti = u0 + ir*opComplex.robert_grad_lon(phi);
					SphereDataComplex Tj = v0 + ir*opComplex.robert_grad_lat(phi);

					lhs = u;

					rhs =	inv_kappa*(
								alpha*Ti
								- f*Tj
							);

					errorCheck(lhs, rhs, "6 lemma inv_kappa f u", epsilon, simVars.rexi.rexi_use_extended_modes < 2);


					lhs = v;

					rhs =	inv_kappa*(
							  f*Ti
							+ alpha*Tj
							);

					errorCheck(lhs, rhs, "6 lemma inv_kappa f v", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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
							  Fp_i*(u0 + ir*opComplex.robert_grad_lon(phi) )
							+ Fp_j*(v0 + ir*opComplex.robert_grad_lat(phi) )
						);

					lhs = F_lhs;
					rhs = F_rhs;

					errorCheck(lhs, rhs, "6 F(u,v) = F(u0,v0,phi)", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
				}

				{
					/*
					 * Test F=... formulation, |*kappa scaled lhs and rhs
					 *
					 * rexi sph ver 15, sec. 5.3.7
					 * Solution for geopotential
					 */

					SphereDataComplex one(sphereDataConfigExt);
					one.physical_set_zero();
					one = one+1.0;

					SphereDataComplex kappa = alpha*alpha+f*f;
//					SphereDataComplex inv_kappa = one/kappa;

					SphereDataComplex Fp_i = fj*(-(alpha*alpha-f*f));
					SphereDataComplex Fp_j = fj*(2.0*alpha*f);

					SphereDataComplex Fc = Fp_i*u0 + Fp_j*v0;

					SphereDataComplex lhs =
							  kappa*kappa*phi
							+ phi_bar/alpha*(Fp_i*ir*opComplex.robert_grad_lon(phi) + Fp_j*ir*opComplex.robert_grad_lat(phi))
							- kappa*phi_bar*ir*ir*opComplex.laplace(phi);

					SphereDataComplex rhs =
							kappa*phi_bar*(div0 - f*(1.0/alpha)*eta0)
							+ kappa*(alpha + f*f*(1.0/alpha))*phi0
							- phi_bar/alpha*Fc;

					errorCheck(lhs, rhs, "7a", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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
							+ phi_bar/alpha*(Fp_i*ir*opComplex.robert_grad_lon(phi) + Fp_j*ir*opComplex.robert_grad_lat(phi))
							- phi_bar*ir*ir*opComplex.laplace(phi);

					SphereDataComplex F_rhs =
							  phi_bar*(div0 - f*(1.0/alpha)*eta0)
							+ (alpha + f*f*(1.0/alpha))*phi0
							- phi_bar/alpha*Fc;

					lhs = F_lhs;
					rhs = F_rhs;

					std::cout << "ERRORS IGNORED WITH 1/kappa formulation!" << std::endl;
					errorCheck(lhs, rhs, "7b", epsilon, true);
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
					lhs = kappa*u;

					SphereDataComplex a = u0 + ir*opComplex.robert_grad_lon(phi);
					SphereDataComplex b = v0 + ir*opComplex.robert_grad_lat(phi);

					rhs = alpha*a -f*b;

					errorCheck(lhs, rhs, "8a", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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

					SphereDataComplex a = u0 + ir*opComplex.robert_grad_lon(phi);
					SphereDataComplex b = v0 + ir*opComplex.robert_grad_lat(phi);

					rhs = f*a + alpha*b;

					errorCheck(lhs, rhs, "8b", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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
							lhs_direct += (2.0*two_omega*two_omega*alpha*alpha)*opComplex.mu2(phi);

							sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);
							lhs_direct += (two_omega*two_omega)*(two_omega*two_omega)*opComplex.mu2(opComplex.mu2(phi));

							sphSolverPhi.solver_component_rexi_z4robert(	-avg_geopotential*alpha*two_omega, r);
							lhs_direct += (-avg_geopotential*alpha*two_omega)*(1.0/(r*r))*/* 1/cos^2phi opComplex.robert_grad_lat(mu)*/ opComplex.robert_grad_lon(phi);

							sphSolverPhi.solver_component_rexi_z5robert(	avg_geopotential/alpha*two_omega*two_omega*two_omega, r);
							lhs_direct += (avg_geopotential/alpha*two_omega*two_omega*two_omega)*(1.0/(r*r))*opComplex.mu2(opComplex.robert_grad_lon(phi));

							sphSolverPhi.solver_component_rexi_z6robert(	avg_geopotential*2.0*two_omega*two_omega, r);
							lhs_direct += (avg_geopotential*2.0*two_omega*two_omega)*(1.0/(r*r))*opComplex.mu(opComplex.robert_grad_lat(phi));
						}

						sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential*alpha*alpha, r);
						lhs_direct += (-avg_geopotential*alpha*alpha)*(1.0/(r*r))*opComplex.laplace(phi);

						//if (use_formulation_with_coriolis_effect)
						{
							sphSolverPhi.solver_component_rexi_z8(	-avg_geopotential*two_omega*two_omega, r);
							lhs_direct += (-avg_geopotential*two_omega*two_omega)*(1.0/(r*r))*opComplex.mu2(opComplex.laplace(phi));
						}


						SphereDataComplex rhs_direct(sphereDataConfigExt);
						{

							auto kappa = [&](
									const SphereDataComplex &i_data
							) -> SphereDataComplex
							{
								return (alpha*alpha)*i_data + two_omega*two_omega*SphereOperatorsComplex::mu2(i_data);
							};

							SphereDataComplex div0 = inv_r*SphereOperatorsComplex::robert_div(u0, v0);
							SphereDataComplex eta0 = inv_r*SphereOperatorsComplex::robert_vort(u0, v0);

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
								+ phi_bar/alpha*(Fp_i*ir*opComplex.robert_grad_lon(phi) + Fp_j*ir*opComplex.robert_grad_lat(phi))
								- kappa*phi_bar*ir*ir*opComplex.laplace(phi);
#endif

						SphereDataComplex phi_reduced = phi.spectral_returnWithDifferentModes(sphereDataConfig);

						errorCheck(computed_phi_lhs_direct, phi_reduced, "test Z solvers LHS phi", epsilon, true);
						errorCheck(computed_phi_rhs_direct, phi_reduced, "test Z solvers RHS phi", epsilon, true);

						std::cout << "OK" << std::endl;
					}

#if 0
					auto kappa = [&](
							const SphereDataComplex &i_data
					)	-> SphereDataComplex
					{
						return (alpha*alpha)*i_data + two_omega*two_omega*SphereOperatorsComplex::mu2(i_data);
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

						errorCheck(phi, prog_u_cplx, "VEL (phi) u_cplx", epsilon, simVars.rexi.rexi_use_extended_modes < 2);

						SphereDataComplex a = u0 + inv_r*SphereOperatorsComplex::robert_grad_lon(phi);
						SphereDataComplex b = v0 + inv_r*SphereOperatorsComplex::robert_grad_lat(phi);

						SphereDataComplex rhsa = alpha*a - two_omega*SphereOperatorsComplex::mu(b);
						SphereDataComplex rhsb = two_omega*SphereOperatorsComplex::mu(a) + alpha*b;

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
								+ phi_bar/alpha*(Fp_i*ir*opComplex.robert_grad_lon(phi) + Fp_j*ir*opComplex.robert_grad_lat(phi))
								- kappa*phi_bar*ir*ir*opComplex.laplace(phi);

						SphereDataComplex rhs =
								kappa*phi_bar*(div0 - f*(1.0/alpha)*eta0)
								+ kappa*(alpha + f*f*(1.0/alpha))*phi0
								- phi_bar/alpha*Fc;
					}

					errorCheck(lhs, rhs, "combined a", epsilon, true);
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

				SphereDataComplex div0 = ir*opComplex.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplex.robert_vort(u0, v0);

				SphereDataComplex div = ir*opComplex.robert_div(u, v);
				SphereDataComplex eta = ir*opComplex.robert_vort(u, v);


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
				SphereDataComplex Ti = u0 + ir*opComplex.robert_grad_lon(prog_phi_cplx);
				SphereDataComplex Tj = v0 + ir*opComplex.robert_grad_lat(prog_phi_cplx);
#else
				SphereDataComplex Ti = u0 + ir*opComplex.robert_grad_lon(Convert_SphereData_To_SphereDataComplex::physical_convert(rexi_prog_phi));
				SphereDataComplex Tj = v0 + ir*opComplex.robert_grad_lat(Convert_SphereData_To_SphereDataComplex::physical_convert(rexi_prog_phi));
#endif

				SphereDataComplex lhs = u;
				SphereDataComplex rhs =	inv_kappa*(
							alpha*Ti
							- two_omega*opComplex.mu(Tj)
						);

				errorCheck(lhs, rhs, "VEL (phi) REXI u", epsilon, simVars.rexi.rexi_use_extended_modes < 2);

				lhs = v;
				rhs =	inv_kappa*(
						  two_omega*opComplex.mu(Ti)
						+ alpha*Tj
						);

				errorCheck(lhs, rhs, "VEL (phi) REXI v", epsilon, simVars.rexi.rexi_use_extended_modes < 2);

				{
					SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverVel;

					sphSolverVel.setup(sphereDataConfig, 2);
					sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
					sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);

					double inv_r = ir;
					SphereDataComplex a = u0 + inv_r*SphereOperatorsComplex::robert_grad_lon(prog_phi_cplx);
					SphereDataComplex b = v0 + inv_r*SphereOperatorsComplex::robert_grad_lat(prog_phi_cplx);

					SphereDataComplex rhsa = alpha*a - two_omega*SphereOperatorsComplex::mu(b);
					SphereDataComplex rhsb = two_omega*SphereOperatorsComplex::mu(a) + alpha*b;

					SphereDataComplex u_cplx = sphSolverVel.solve(rhsa);
					SphereDataComplex v_cplx = sphSolverVel.solve(rhsb);


					errorCheck(u_cplx, prog_u_cplx, "VEL (phi) u_cplx", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
					errorCheck(v_cplx, prog_v_cplx, "VEL (phi) v", epsilon, simVars.rexi.rexi_use_extended_modes < 2);
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
					SWERexi_SPHRobert rexiSPHRobert;

					rexiSPHRobert.setup(
							sphereDataConfigExt,
							sphereDataConfig,
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
						errorCheck(rexi_prog_phi_cplx_ext, prog_phi_cplx_ext, "prog_phi SSSSSSSSSSS", epsilon*1e+1);
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
					SWERexi_SPH rexiSPH;

					rexiSPH.setup(
							sphereDataConfigExt,
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
			errorCheck(rexi_prog_phi, prog_phi, "prog_phi", epsilon*1e+1);
			errorCheck(rexi_prog_u, prog_u, "prog_u", epsilon*1e+1);
			errorCheck(rexi_prog_v, prog_v, "prog_v", epsilon*1e+1);
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
