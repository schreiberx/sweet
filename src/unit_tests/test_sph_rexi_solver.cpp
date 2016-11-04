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
	std::cout << "Using max allowed error of " << epsilon << std::endl;

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

	SphereTestSolutions_Gaussian testSolutions_phi_re(M_PI/1.0, M_PI/2.0), testSolutions_phi_im(M_PI/3.0, M_PI/4.0);
	SphereTestSolutions_Gaussian testSolutions_u_re(M_PI/5.0, M_PI/7.0), testSolutions_u_im(M_PI/9.0, -M_PI/12.0);
	SphereTestSolutions_Gaussian testSolutions_v_re(M_PI/6.0, -M_PI/8.0), testSolutions_v_im(M_PI/10.0, M_PI/11.0);


	REXI rexi;
	rexi.setup(simVars.rexi.rexi_h, simVars.rexi.rexi_M);

	if (!simVars.misc.sphere_use_robert_functions)
	{
		std::cerr << "WARNING: TESTS WILL FAIL WITH NON-ROBERT FUNCTIONS DUE TO SINGULARITY AT POLES!" << std::endl;
	}

	double timestep_size = 1.0;

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

			prog_phi_cplx.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
						double tmp;
						testSolutions_phi_re.test_function__grid_gaussian(lat,mu,tmp);
						io_data.real(tmp);

						if (!use_complex_valued_solver)
						{
							io_data.imag(0);
						}
						else
						{
							double tmp;
							testSolutions_phi_im.test_function__grid_gaussian(lat,mu,tmp);
							io_data.imag(tmp);
						}

//						if (simVars.misc.sphere_use_robert_functions)
//							io_data *= std::sqrt(1.0-mu*mu);
					}
			);

			prog_u_cplx.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
#if 0
						io_data = {0,0};
#else
						double tmp;
						testSolutions_u_re.test_function__grid_gaussian(lat,mu,tmp);
						io_data.real(tmp);

						if (!use_complex_valued_solver)
						{
							io_data.imag(0);
						}
						else
						{
							double tmp;
							testSolutions_u_im.test_function__grid_gaussian(lat,mu,tmp);
							io_data.imag(tmp);
						}
#endif

						if (simVars.misc.sphere_use_robert_functions)
							io_data *= std::sqrt(1.0-mu*mu);
					}
			);

			prog_v_cplx.physical_update_lambda_gaussian_grid(
					[&](double lat, double mu, std::complex<double> &io_data)
					{
#if 0
						io_data = {0,0};
#else
						double tmp;
						testSolutions_v_re.test_function__grid_gaussian(lat,mu,tmp);
						io_data.real(tmp);

						if (!use_complex_valued_solver)
						{
							io_data.imag(0);
						}
						else
						{
							double tmp;
							testSolutions_v_im.test_function__grid_gaussian(lat,mu,tmp);
							io_data.imag(tmp);
						}
#endif

						if (simVars.misc.sphere_use_robert_functions)
							io_data *= std::sqrt(1.0-mu*mu);
					}
			);

			prog_phi_cplx.spectral_truncate();
			prog_u_cplx.spectral_truncate();
			prog_v_cplx.spectral_truncate();

			/*
			 * Computed right hand side
			 */
			SphereDataComplex prog_phi0_cplx(sphereDataConfig);
			SphereDataComplex prog_u0_cplx(sphereDataConfig);
			SphereDataComplex prog_v0_cplx(sphereDataConfig);

			double phi_bar = simVars.sim.gravitation*simVars.sim.h0;


			/*
			 * Compute RHS for REXI term
			 */
			if (simVars.misc.sphere_use_robert_functions)
			{
				prog_phi0_cplx =
							+ alpha*prog_phi_cplx
							- phi_bar*ir*opComplex.robert_div_lon(prog_u_cplx)
							- phi_bar*ir*opComplex.robert_div_lat(prog_v_cplx);

				prog_u0_cplx =
							- ir*opComplex.robert_grad_lon(prog_phi_cplx)
							+ alpha*prog_u_cplx
							+ two_omega*opComplex.mu(prog_v_cplx);

				prog_v0_cplx =
							- ir*opComplex.robert_grad_lat(prog_phi_cplx)
							- two_omega*opComplex.mu(prog_u_cplx)
							+ alpha*prog_v_cplx;
			}
			else
			{
				prog_phi0_cplx =
							+ alpha*prog_phi_cplx
							- phi_bar*ir*opComplex.div_lon(prog_u_cplx)
							- phi_bar*ir*opComplex.div_lat(prog_v_cplx);

				prog_u0_cplx =
							- ir*opComplex.grad_lon(prog_phi_cplx)
							+ alpha*prog_u_cplx
							+ two_omega*opComplex.mu(prog_v_cplx);

				prog_v0_cplx =
							- ir*opComplex.grad_lat(prog_phi_cplx)
							- two_omega*opComplex.mu(prog_u_cplx)
							+ alpha*prog_v_cplx;
			}

			/*
			 * Forward calculated output (reference solution)
			 */
			SphereDataComplex prog_phi0_cplx_ext(sphereDataConfigRexi);
			SphereDataComplex prog_u0_cplx_ext(sphereDataConfigRexi);
			SphereDataComplex prog_v0_cplx_ext(sphereDataConfigRexi);

			SphereData prog_phi0_ext(sphereDataConfigRexi);
			SphereData prog_u0_ext(sphereDataConfigRexi);
			SphereData prog_v0_ext(sphereDataConfigRexi);


			prog_phi0_cplx_ext = prog_phi0_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
			prog_u0_cplx_ext = prog_u0_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);
			prog_v0_cplx_ext = prog_v0_cplx.spectral_returnWithDifferentModes(sphereDataConfigRexi);

			if (!use_complex_valued_solver)
			{
				SphereData prog_phi0 = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_phi0_cplx);
				SphereData prog_u0 = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_u0_cplx);
				SphereData prog_v0 = Convert_SphereDataComplex_To_SphereData::physical_convert(prog_v0_cplx);

				prog_phi0_ext = prog_phi0.spectral_returnWithDifferentModes(sphereDataConfigRexi);
				prog_u0_ext = prog_u0.spectral_returnWithDifferentModes(sphereDataConfigRexi);
				prog_v0_ext = prog_v0.spectral_returnWithDifferentModes(sphereDataConfigRexi);
			}



#if 1
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


				//
				SphereDataComplex grad_lat_f(sphereDataConfigRexi);
				grad_lat_f.physical_update_lambda_gaussian_grid(
						[&](double lon, double mu, std::complex<double> &o_data)
						{
							o_data = std::sqrt(1-mu*mu)*two_omega;
						}
					);


				SphereDataComplex div0 = ir*opComplex.robert_div(u0, v0);
				SphereDataComplex eta0 = ir*opComplex.robert_vort(u0, v0);
				SphereDataComplex etad0 = ir*opComplex.robert_div_lon(v0) - ir*opComplex.robert_div_lat(u0);

				SphereDataComplex div = ir*opComplex.robert_div(u, v);
				SphereDataComplex eta = ir*opComplex.robert_vort(u, v);
				SphereDataComplex etad = ir*opComplex.robert_div_lon(v) - ir*opComplex.robert_div_lat(u);

#if 0
				SphereDataComplex fi = ir*opComplex.robert_div_lon(f);
				SphereDataComplex fj = ir*opComplex.robert_div_lat(f);
#else
				SphereDataComplex fi = ir*opComplex.robert_grad_lon(f);
				SphereDataComplex fj = ir*opComplex.robert_grad_lat(f);
#endif
				SphereDataComplex lhs(sphereDataConfigRexi);
				SphereDataComplex rhs(sphereDataConfigRexi);


				{
					/*
					 * REXI SPH document ver 14, Section 5.2, eq (1)
					 */

					lhs = alpha*phi - phi_bar*ir*opComplex.robert_div_lon(u) - phi_bar*1.0/r*opComplex.robert_div_lat(v);
					rhs = phi0;

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 0a: " << error << std::endl;

					if (error > epsilon*10e+4)
						FatalError("Error too large");
				}


				{
					/*
					 * REXI SPH document ver 14, Section 5.2, eq (2)
					 */

					lhs = -ir*opComplex.robert_grad_lon(phi) + alpha*u + f*v;
					rhs = u0;

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 0b: " << error << std::endl;

					if (error > epsilon)
						FatalError("Error too large");
				}
				{
					/*
					 * REXI SPH document ver 14, Section 5.2, eq (3)
					 */

					lhs = -ir*opComplex.robert_grad_lat(phi) - f*u + alpha*v;
					rhs = v0;

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 0c: " << error << std::endl;

					if (error > epsilon)
						FatalError("Error too large");
				}
				{
					lhs = div;
					rhs = ir*opComplex.robert_div_lon(u) + ir*opComplex.robert_div_lat(v);

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 0d: " << error << std::endl;

					if (error > epsilon)
						FatalError("Error too large");
				}
				{
					lhs = eta;
					rhs = ir*opComplex.robert_div_lon(v) - ir*opComplex.robert_div_lat(u);

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 0e: " << error << std::endl;

					if (error > epsilon)
						FatalError("Error too large");
				}
#if 1
				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Identity relation
					 */
					lhs = ir*opComplex.robert_div(phi*u, phi*v);
					rhs = phi*ir*opComplex.robert_div(u, v) + u*ir*opComplex.robert_grad_lon(phi) + v*ir*opComplex.robert_grad_lat(phi);


					double error_max = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					double error_rms = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_rms();

					std::cout << "ERROR 0f robert GRAD f: " << error_max << ", " << error_rms << std::endl;

					if (error_max > epsilon*10e+5)
						std::cerr << "Error too large - ignoring" << std::endl;;
				}
#endif
#if 0
				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * Identity relation
					 */
					lhs = ir*opComplex.robert_div(f*u, f*v);
					rhs = f*ir*opComplex.robert_div(u, v)
							+ u*fi	/// TODO: UNCLEAR WHY WE HAVE TO USE DIV HERE!!!!!!!!!!
							+ v*fj;	/// TODO: UNCLEAR WHY WE HAVE TO USE DIV HERE!!!!!!!!!!


					double error_max = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					double error_rms = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_rms();

					std::cout << "ERROR 0f robert DIV f: " << error_max << ", " << error_rms << std::endl;

					if (error_max > epsilon*10e+5)
						FatalError("Error too large");
				}
#endif
				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st equation in this section
					 */

					lhs = ir*opComplex.robert_div_lon(
								-ir*opComplex.robert_grad_lon(phi) + alpha*u + f*v
							  );

					rhs = ir*opComplex.robert_div_lon(u0);



					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 1aa: " << error << std::endl;

					if (error > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhs).physical_write_file("o_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhs).physical_write_file("o_rhs.csv");
						FatalError("Error too large");
					}
				}

				{
					/*
					 * REXI SPH document ver 14, Section 5.2.1
					 * 1st equation in this section
					 */

					lhs = ir*opComplex.robert_div_lon(v0);

					rhs = ir*opComplex.robert_div_lon(
							-ir*opComplex.robert_grad_lon(phi) - f*u + alpha*v
						  );

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 1ab: " << error << std::endl;

					if (error > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhs).physical_write_file("o_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhs).physical_write_file("o_rhs.csv");
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

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 1a: " << error << std::endl;

					if (error > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhs).physical_write_file("o_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhs).physical_write_file("o_rhs.csv");
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

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 1ba: " << error << std::endl;

					if (error > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhs).physical_write_file("o_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhs).physical_write_file("o_rhs.csv");
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

							+ ir*opComplex.robert_div_lon(	f*v		)
							+ ir*opComplex.robert_div_lat(	-f*u	)
							;

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 1ba: " << error << std::endl;

					if (error > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhs).physical_write_file("o_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhs).physical_write_file("o_rhs.csv");
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
							+ ir*f*opComplex.robert_vort(u, v)
							- u*fj + v*fi;

					rhs = div0;

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 1c: " << error << std::endl;

					if (error > epsilon*10e+2)
					{
						Convert_SphereDataComplex_To_SphereData::physical_convert(lhs).physical_write_file("o_lhs.csv");
						Convert_SphereDataComplex_To_SphereData::physical_convert(rhs).physical_write_file("o_rhs.csv");
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

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 2: " << error << std::endl;

					if (error > epsilon*10e+3)
						FatalError("Error too large");
				}



				{
					/*
					 * REXI SPH document ver 14, Section 5.2.3
					 * eq. (6)
					 */

					lhs = div;
					rhs = 1.0/phi_bar*(phi*alpha-phi0);

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 3: " << error << std::endl;

					if (error > epsilon*10e+3)
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

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 4a: " << error << std::endl;

					if (error > epsilon*10e+3)
						FatalError("Error too large");
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

					double error = (
										lhs.spectral_returnWithDifferentModes(sphereDataConfig)
										- rhs.spectral_returnWithDifferentModes(sphereDataConfig)
									).physical_reduce_error_max();

					std::cout << "ERROR 4b: " << error << std::endl;

					if (error > epsilon*10e+3)
						FatalError("Error too large");
				}



				{
					/*
					 * REXI SPH document ver 14, Section 5.2.5 at the end
					 * eq. (9)
					 */
					SphereDataComplex F =
							f*(opComplex.robert_grad_lon(f)*v0)
							- alpha*(opComplex.robert_grad_lon(f)*u0);


					lhs = (alpha*alpha*phi + f*f*phi) - phi_bar*opComplex.laplace(phi) + (phi_bar/alpha) * F;
					rhs = phi_bar*(div0 - f*(1.0/alpha)*eta0) + (alpha+f*f*(1.0/alpha))*phi0;

					lhs = alpha*alpha*phi - phi_bar*opComplex.laplace(phi);
					rhs = phi_bar*div0 + alpha*phi0;

					double error = (lhs.spectral_returnWithDifferentModes(sphereDataConfig) - rhs.spectral_returnWithDifferentModes(sphereDataConfig)).physical_reduce_error_max();

					std::cout << "ERROR 5: " << error << std::endl;


					if (error > epsilon*10e+3)
						FatalError("Error too large");
				}

			}
#endif


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
				prog_phi.file_physical_writeFile_lon_pi_shifted(obuf.str().c_str());
			}
			{
				std::ostringstream obuf;
				obuf << "prog_alphaid" << i << "_phi_rexi.csv";
				rexi_prog_phi.file_physical_writeFile_lon_pi_shifted(obuf.str().c_str());
			}
			{
				std::ostringstream obuf;
				obuf << "prog_alphaid" << i << "_diff.csv";
				(rexi_prog_phi-prog_phi).file_physical_writeFile_lon_pi_shifted(obuf.str().c_str());
			}
#endif
			/*
			 * REXI results stored in rexi_prog_*
			 * These should match prog_*
			 */
			double err_phi = (rexi_prog_phi - prog_phi).reduce_abs_max();
			double err_u = (rexi_prog_u - prog_u).reduce_abs_max();
			double err_v = (rexi_prog_v - prog_v).reduce_abs_max();

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
