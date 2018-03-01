/*
 * BenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>
//#include <benchmarks_sphere/BenchmarkGalewsky.hpp>
#include <benchmarks_sphere/BenchmarkGaussianDam.hpp>
#include <benchmarks_sphere/BenchmarkFlowOverMountain.hpp>

class SphereBenchmarksCombined
{

public:
	static void printAvailableBenchmarks()
	{
		std::cout << std::endl;
		std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
		std::cout << "Available benchmark scenarios:" << std::endl;
		std::cout << "	0: Dummy" << std::endl;
		std::cout << "	2: Use Gaussian bump initial conditions (0, pi/3, exp=10)" << std::endl;
		std::cout << "	3: Use Gaussian bump initial conditions (pi/3, pi/3, exp=10)" << std::endl;
		std::cout << "	4: Use Gaussian bump initial conditions (-pi/2, pi/4, exp=10)" << std::endl;
		std::cout << "	4: Use Gaussian bump initial conditions (-pi/2, pi/4, exp=10)" << std::endl;
		std::cout << "	5: Use Gaussian bump initial conditions (-pi/2, pi/4, exp=10)" << std::endl;
		std::cout << "	6: Use Gaussian bump initial conditions (pi/3, pi/3, exp=20)" << std::endl;
		std::cout << "	7: Use Gaussian bump initial conditions (-pi, pi/4, exp=100)" << std::endl;
		std::cout << "	9: Combination of Gaussian bumps" << std::endl;

		std::cout << "	10: Williamson Geostrophic balance test case, unit sphere" << std::endl;
		std::cout << "	11: Williamson Geostrophic balance test case, earth scale" << std::endl;

		std::cout << "	21: Williamson Advection test benchmark 1 (DIV formulation)" << std::endl;
		std::cout << "	22: Williamson Advection test benchmark 1 (GRAD formulation)" << std::endl;

		std::cout << "	50: Swartztrauber 2004 Gaussian breaking dam" << std::endl;
		std::cout << "	100: Galweski" << std::endl;
		std::cout << "	101: Galweski - geostrophic case including non-linear parts" << std::endl;
		std::cout << "	200: h=h0, u=0, v=0" << std::endl;
		std::cout << std::endl;
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	static void computeGeostrophicBalance(
			SphereData &o_h,
			SphereDataPhysical &i_u,
			SphereDataPhysical &i_v,

			SimulationVariables &i_simVars,
			SphereOperators &i_op
	)
	{
		/*
		 * Setup Coriolis effect
		 */
		SphereDataPhysical f(o_h.sphereDataConfig);
		f.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = 2.0*i_simVars.sim.coriolis_omega*mu;
			}
		);

		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereData vrtspec(o_h.sphereDataConfig);
		SphereData divspec(o_h.sphereDataConfig);

		if (i_simVars.misc.sphere_use_robert_functions)
		  i_op.robert_uv_to_vortdiv(i_u, i_v, vrtspec, divspec);
		else
		  i_op.uv_to_vortdiv(i_u, i_v, vrtspec, divspec);

		SphereDataPhysical vrtg = vrtspec.getSphereDataPhysical();

		SphereDataPhysical tmpg1 = i_u*(vrtg+f);
		SphereDataPhysical tmpg2 = i_v*(vrtg+f);

		SphereData tmpspec1(o_h.sphereDataConfig);
		SphereData tmpspec2(o_h.sphereDataConfig);

		if (i_simVars.misc.sphere_use_robert_functions)
		  i_op.robert_uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);
		else
		  i_op.uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);
		
		tmpspec2 = 0.5*(i_u*i_u+i_v*i_v);

		SphereData phispec = i_op.inv_laplace(tmpspec1) - tmpspec2 - tmpspec1;

		phispec = phispec + i_simVars.sim.gravitation*i_simVars.sim.h0;
		phispec.spectral_truncate();

		o_h = phispec/i_simVars.sim.gravitation;
	}


public:

        static 
	void setupTopography(
			     SimulationVariables &io_simVars,
			     SphereOperators &i_op
			     )
         {
	   if (io_simVars.setup.benchmark_scenario_name == "flow_over_mountain") 
	     {
	       // set the topography flag to true
	       io_simVars.sim.use_topography = true;

	       // setup the parameters for the flow-over-mountain test case
	       const double R            = M_PI/9.;
	       const double h_topo_0     = 2000.;
	       const double i_center_lon = 3.*M_PI/2.;
	       const double i_center_lat = M_PI/6.;

	       // initialize the topography
	       io_simVars.sim.h_topo.physical_set_zero();

	       // setup the topography vector
	       BenchmarkFlowOverMountain::setup_topography(io_simVars.sim.h_topo,
							   io_simVars,
							   R,
							   h_topo_0,
							   i_center_lon,
							   i_center_lat
							   );

	       io_simVars.sim.h_topo.spectral_truncate();
	     }
	   else 
	     {
	       // set the topography flag to false
	       io_simVars.sim.use_topography = false;
	     }
         }



	static
	void setupInitialConditions(
			SphereData &o_h,
			SphereData &o_u,
			SphereData &o_v,

			SimulationVariables &io_simVars,
			SphereOperators &i_op
	)
	{
		if (io_simVars.setup.benchmark_scenario_name == "")
		{
			if (io_simVars.setup.benchmark_scenario_id < 0)
			{
				printAvailableBenchmarks();
				FatalError("Benchmark scenario not selected");
			}

			if (io_simVars.setup.benchmark_scenario_id == 0)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 2)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, 0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 3)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 4)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/4.0, -M_PI);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 5)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/4.0, -M_PI, 100.0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 6)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0, 20.0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 7)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/4.0, -M_PI, 100.0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 9)
			{
				SphereData tmp(o_h.sphereDataConfig);

				o_u.physical_set_all_value(0);
				o_v.physical_set_all_value(0);

				o_h.physical_set_all_value(io_simVars.sim.h0);

				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, io_simVars, 2.0*M_PI*0.1, M_PI/3, 20.0);
				o_h += (tmp-io_simVars.sim.h0);
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, io_simVars, 2.0*M_PI*0.6, M_PI/5.0, 50.0);
				o_h += (tmp-io_simVars.sim.h0);
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, io_simVars, 2.0*M_PI*0.8, -M_PI/4, 100.0);
				o_h += (tmp-io_simVars.sim.h0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 10 || io_simVars.setup.benchmark_scenario_id == 11)
			{
				/*
				 * Williamson test case 2 for geostrophic balance.
				 *
				 * See Williamson paper for accurate setup
				 */

				double u0 = 1.0;
				if (io_simVars.setup.benchmark_scenario_id == 11)
				{
					if (io_simVars.timecontrol.current_simulation_time == 0)
					{
						std::cout << "!!! WARNING !!!" << std::endl;
						std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
						std::cout << "!!! WARNING !!!" << std::endl;
					}

					io_simVars.sim.coriolis_omega = 7.292e-5;
					io_simVars.sim.gravitation = 9.80616;
					io_simVars.sim.earth_radius = 6.37122e6;
					io_simVars.sim.h0 = 29400.0/io_simVars.sim.gravitation;

					u0 = (2.0*M_PI*io_simVars.sim.earth_radius)/(12.0*24.0*60.0*60.0);
				}
				double a = io_simVars.sim.earth_radius;

				double phi0 = io_simVars.sim.h0*io_simVars.sim.gravitation;

				o_h.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						io_data = (phi0 + (a*io_simVars.sim.coriolis_omega*u0*std::cos(i_lat)*std::cos(i_lat)))/io_simVars.sim.gravitation;
					}
				);

				if (io_simVars.misc.sphere_use_robert_functions)
				{
					o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							io_data = u0*std::cos(i_lat)*std::cos(i_lat);
						}
					);
				}
				else
				{
					o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							io_data = std::cos(i_lat);
						}
					);
				}

				o_v.spectral_set_zero();

#if 0
				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					double inv_r = 1.0/io_simVars.sim.earth_radius;
					double h_max_error, u_max_error, v_max_error;

					if (io_simVars.misc.sphere_use_robert_functions)
					{
						h_max_error = (-inv_r*i_op.robert_div_lon(o_u) -inv_r*i_op.robert_div_lat(o_v)).physical_reduce_max_abs();
						u_max_error = (-inv_r*i_op.robert_grad_lon(o_h*io_simVars.sim.gravitation) + 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_v)).physical_reduce_max_abs();
						v_max_error = (-inv_r*i_op.robert_grad_lat(o_h*io_simVars.sim.gravitation) - 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_u)).physical_reduce_max_abs();
					}
					else
					{
						h_max_error = (-inv_r*i_op.div_lon(o_u) -inv_r*i_op.div_lat(o_v)).physical_reduce_max_abs();
						u_max_error = (-inv_r*i_op.grad_lon(o_h*io_simVars.sim.gravitation) + 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_v)).physical_reduce_max_abs();
						v_max_error = (-inv_r*i_op.grad_lat(o_h*io_simVars.sim.gravitation) - 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_u)).physical_reduce_max_abs();
					}

					std::cout << "h_max_error for geostrophic balance case: " << h_max_error << std::endl;
					std::cout << "u_max_error for geostrophic balance case: " << u_max_error << std::endl;
					std::cout << "v_max_error for geostrophic balance case: " << v_max_error << std::endl;
				}
#endif
			}
			else if (io_simVars.setup.benchmark_scenario_id == 21)
			{
				/*
				 * Advection test case
				 * See Williamson test case, eq. (75), (76)
				 *
				 * phi_t = DIV (u phi)
				 */

				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				io_simVars.sim.coriolis_omega = 7.292e-5;
				io_simVars.sim.gravitation = 9.80616;
				io_simVars.sim.earth_radius = 6.37122e6;
				io_simVars.sim.h0 = 1000.0;

				double lambda_c = 3.0*M_PI/2.0;
				double theta_c = M_PI/4.0;
				double a = io_simVars.sim.earth_radius;

				double R = a/3.0;
				double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

				o_h.physical_update_lambda(
					[&](double i_lambda, double i_theta, double &io_data)
				{
						double r = a * std::acos(
								std::sin(theta_c)*std::sin(i_theta) +
								std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
						);

						if (r < R)
							io_data = io_simVars.sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
						else
							io_data = 0;
					}
				);


				if (io_simVars.misc.sphere_use_robert_functions)
				{
					o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							double i_theta = i_lat;
							double i_lambda = i_lon;
							io_data =
								u0*(
										std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
										std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
									);

							io_data *= std::cos(i_lat);
						}
					);

					o_v.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
	//						double i_phi = i_lat;
							double i_lambda = i_lon;
							io_data =
								-u0*(
										std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
									);

							io_data *= std::cos(i_lat);
						}
					);
				}
				else
				{
					o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							double i_theta = i_lat;
							double i_lambda = i_lon;
							io_data =
								u0*(
									std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
									);
						}
					);

					o_v.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
	//						double i_phi = i_lat;
							double i_lambda = i_lon;
							io_data =
								-u0*(
									std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
								);
						}
					);
				}

				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Storing advection in output velocities !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);

				std::cout << "advection_rotation_angle: " << io_simVars.setup.advection_rotation_angle << std::endl;
			}
			else if (io_simVars.setup.benchmark_scenario_id == 22)
			{
				/**
				 * Advection test case
				 * See Williamson test case, eq. (75), (76)
				 *
				 * phi_t = phi GRAD u
				 */

				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				io_simVars.sim.coriolis_omega = 7.292e-5;
				io_simVars.sim.gravitation = 9.80616;
				io_simVars.sim.earth_radius = 6.37122e6;
				io_simVars.sim.h0 = 1000.0;

				double lambda_c = 3.0*M_PI/2.0;
				double theta_c = 0;
				double a = io_simVars.sim.earth_radius;

				double R = a/3.0;
				double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

				o_h.physical_update_lambda(
					[&](double i_lambda, double i_theta, double &io_data)
					{
						double r = a * std::acos(
								std::sin(theta_c)*std::sin(i_theta) +
								std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
						);

						if (r < R)
							io_data = io_simVars.sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
						else
							io_data = 0;
					}
				);


				if (io_simVars.misc.sphere_use_robert_functions)
				{
					o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							double i_theta = i_lat;
							double i_lambda = i_lon;
							io_data =
									u0*(
										std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
										std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
								);

							io_data *= std::cos(i_lat);
						}
					);

					o_v.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
	//						double i_phi = i_lat;
							double i_lambda = i_lon;
							io_data =
								-u0*(
										std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
								);

							io_data *= std::cos(i_lat);
						}
					);
				}
				else
				{
					o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							double i_theta = i_lat;
							double i_lambda = i_lon;
							io_data =
									u0*(
										std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
										std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
								);
						}
					);

					o_v.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
	//						double i_phi = i_lat;
							double i_lambda = i_lon;
							io_data =
								-u0*(
										std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
								);
						}
					);
				}

				o_h.physical_truncate();
				o_u.physical_truncate();
				o_v.physical_truncate();

				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Storing advection in output velocities !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}


				io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);

				std::cout << "advection_rotation_angle: " << io_simVars.setup.advection_rotation_angle << std::endl;
			}
			else if (io_simVars.setup.benchmark_scenario_id == 40)
			{
				o_h.physical_set_all_value(io_simVars.sim.h0);
				o_u.physical_set_all_value(0);
				o_v.physical_set_all_value(0);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 50)
			{
				/*
				 * PAUL N. SWARZTRAUBER
				 * Shallow Water Flow on the Sphere
				 * 2004
				 */
				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				if (io_simVars.sim.coriolis_omega != 0)
					io_simVars.sim.coriolis_omega = 7.292e-5;

				io_simVars.sim.gravitation = 9.80616;
				io_simVars.sim.earth_radius = 6.37122e6;
				io_simVars.sim.h0 = 29400.0;

				o_u.spectral_set_zero();
				o_v.spectral_set_zero();

				const double a = io_simVars.sim.earth_radius;
				const double A = 6000.0;
				const double alpha = 10;

				const double center_lat = M_PI/4;
				const double center_lon = M_PI;

				o_h.physical_update_lambda(
					[&](double i_lambda, double i_phi, double &io_data)
					{
					  double x = a* ( std::cos(i_phi)*std::cos(i_lambda) - std::cos(center_lat)*std::cos(center_lon) );
					  double y = a* ( std::cos(i_phi)*std::sin(i_lambda) - std::cos(center_lat)*std::sin(center_lon) );
					  double z = a* ( std::sin(i_phi)                    - std::sin(center_lat) );

						double d = std::sqrt(x*x+y*y+z*z);

						io_data = io_simVars.sim.h0 + A*std::exp(-alpha*(d/a)*(d/a));
					}
				);
			}
			else if (io_simVars.setup.benchmark_scenario_id == 1 || io_simVars.setup.benchmark_scenario_id == 100 || io_simVars.setup.benchmark_scenario_id == 101)
			{
				//BenchmarkGalewsky benchmarkGalewsky(io_simVars);

				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				/// Setup Galewski parameters
				io_simVars.sim.coriolis_omega = 7.292e-5;
				io_simVars.sim.gravitation = 9.80616;
				io_simVars.sim.earth_radius = 6.37122e6;
				io_simVars.sim.h0 = 10000;

				io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);

#if 0

				benchmarkGalewsky.setup_initial_h(o_h);
	#if 0
				std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
				std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
				std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
				std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
				std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
				std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
				/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
				/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
				/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
				/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
				/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
				/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
				double foo[100] = {
						9055.493745200144986,  9055.465517385822750,  9055.482442174412427,  9055.485353917110842,  9055.470071095474850,  9055.482225177744112,  9055.483738383378295,  9055.470977648797998,  9055.482690202114100,  9055.483017123375248,  9055.470682326214956,  9055.483931642287644,  9055.482260546816178,  9055.469577987099910,  9055.486730172287935,  9055.479593312960787,  9055.470174185988071,  9055.506848781622466,  9056.095270992454971,  9061.214200969297963,  9081.885216784596196,  9134.617053737521928,  9233.437475196804371,  9380.388945759646958,  9560.630706589838155,  9746.223428225062889,  9907.553693800491601,  10025.289612307593416,  10095.839141933991414,  10128.926698702083740,  10140.037739517365480,  10142.276285688618373,  10142.446739097140380,  10142.457887177219163,  10142.460637987251175,  10142.454087103898928,  10142.457439067522500,  10142.458913848522570,  10142.455881685920758,  10142.456781481529106,  10142.458514041243689,  10142.456551124038015,  10142.456497322984433,  10142.458327286207350,  10142.456871132475499,  10142.456378846281950,  10142.458196108580523,  10142.457058206358852,  10142.456330244574929,  10142.458089350835508,  10142.457184941487867,  10142.456312120679286,  10142.457997722302025,  10142.457280603905019,  10142.456307374064636,  10142.457916977255081,  10142.457358891591866,  10142.456308088640071,  10142.457844461576315,  10142.457427072755308,  10142.456310368295817,  10142.457778217069063,  10142.457489475389593,  10142.456312155902197,  10142.457716672483002,  10142.457549000710060,  10142.456312239946783,  10142.457658493196504,  10142.457607867550905,  10142.456309759470969,  10142.457602479804336,  10142.457668042088699,  10142.456303918996127,  10142.457547479152709,  10142.457731553769918,  10142.456293774963342,  10142.457492286770503,  10142.457800817175666,  10142.456278003344778,  10142.457435515911129,  10142.457879081055580,  10142.456254545904812,  10142.457375387479260,  10142.457971214398640,  10142.456219940253504,  10142.457309336608887,  10142.458085327591107,  10142.456167834954613,  10142.457233155662834,  10142.458236684131407,  10142.456085108484331,  10142.457138761870738,  10142.458459248211511,  10142.455939082734403,  10142.457006667080350,  10142.458852741790906,  10142.455614339081876,  10142.456765204698968,  10142.459972075210317,  10142.453909706555351
				};

				o_h.physical_update_lambda_array(
						[&](int i_lon, int i_lat, double &io_data)
						{
							io_data = foo[i_lat];
						}
				);
	#endif

				if (io_simVars.setup.benchmark_scenario_id != 101)
					benchmarkGalewsky.setup_initial_h_add_bump(o_h);

				benchmarkGalewsky.setup_initial_u(o_u);
				benchmarkGalewsky.setup_initial_v(o_v);
	#if 0
				std::cout << o_u.sphereDataConfig->shtns->nlm << " " << o_u.sphereDataConfig->spectral_array_data_number_of_elements << std::endl;
				for (int i = 0; i < o_u.sphereDataConfig->physical_num_lat; i++)
					std::cout << i << ": " << o_u.physical_get(0, i) << std::endl;
				o_u.physical_file_write("o_u_0_physical.csv");
				o_u.request_data_spectral();
				o_u.physical_file_write("o_u_1_spectral_physical.csv");

				std::cout << "----" << std::endl;
				for (int i = 0; i < o_u.sphereDataConfig->physical_num_lat; i++)
					std::cout << i << ": " << o_u.physical_get(0, i) << std::endl;
	#endif

				if (io_simVars.misc.sphere_use_robert_functions)
				{
					o_u.physical_update_lambda_cosphi_grid(
						[&](double i_lon, double i_cosphi, double &io_data)
						{
							io_data *= i_cosphi;
						}
					);

					o_v.physical_update_lambda_cosphi_grid(
						[&](double i_lon, double i_cosphi, double &io_data)
						{
							io_data *= i_cosphi;
						}
					);
				}
	/*
				o_h.spectral_truncate();
				o_u.spectral_truncate();
				o_v.spectral_truncate();
	*/

#else

				/*
				 * Initialization code from Jeffrey Whitaker
				 */

				double umax = 80.;
				double phi0 = M_PI/7.;
				double phi1 = 0.5*M_PI - phi0;
				double phi2 = 0.25*M_PI;
				double en = std::exp(-4.0/std::pow((phi1-phi0), 2.0));
				double alpha = 1./3.;
				double beta = 1./15.;
				double hamp = 120.;

				SphereDataPhysical f(o_h.sphereDataConfig);
				f.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, double &o_data)
					{
						o_data = 2.0*io_simVars.sim.coriolis_omega*mu;
					}
				);

				SphereDataPhysical vg(o_h.sphereDataConfig);
				vg.physical_set_zero();

				//u1 = (umax/en)*np.exp(1./((x.lats-phi0)*(x.lats-phi1)))
				SphereDataPhysical u1(o_h.sphereDataConfig);
				u1.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
					}
				);

				SphereDataPhysical ug = u1;
				ug.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						if (phi >= phi1 || phi <= phi0)
							o_data = 0;
					}
				);

				SphereDataPhysical hbump(o_h.sphereDataConfig);

				if (io_simVars.setup.benchmark_scenario_id != 101)
				{
					hbump.physical_update_lambda(
						[&](double lon, double phi, double &o_data)
						{
							o_data = hamp*std::cos(phi)*std::exp(-(lon/alpha)*(lon/alpha))*std::exp(-(phi2-phi)*(phi2-phi)/beta);
						}
					);
				}
				else
				{
					hbump.physical_set_zero();
				}



				SphereData vrtspec(o_h.sphereDataConfig);
				SphereData divspec(o_h.sphereDataConfig);

				i_op.uv_to_vortdiv(ug, vg, vrtspec, divspec);


				SphereDataPhysical vrtg = vrtspec.getSphereDataPhysical();

				SphereDataPhysical tmpg1 = ug*(vrtg+f);
				SphereDataPhysical tmpg2 = vg*(vrtg+f);

				SphereData tmpspec1(o_h.sphereDataConfig);
				SphereData tmpspec2(o_h.sphereDataConfig);
				i_op.uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

				tmpspec2 = (0.5*(ug*ug+vg*vg));

				SphereData phispec = i_op.inv_laplace(tmpspec1) - tmpspec2;

				io_simVars.sim.h0 = 10000.0;
				SphereDataPhysical phig = io_simVars.sim.gravitation*(hbump+io_simVars.sim.h0) + phispec.getSphereDataPhysical();

				phispec = phig;
				phispec.spectral_truncate();

				////////////////////////////////

				i_op.vortdiv_to_uv(vrtspec, divspec, ug, vg);

				o_h = phispec/io_simVars.sim.gravitation;
				o_u = ug;
				o_v = vg;

				if (io_simVars.misc.sphere_use_robert_functions)
				{
					o_u.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					o_v.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);
				}
#endif

			}
			else if (io_simVars.setup.benchmark_scenario_id == 200)
			{
				o_h.physical_set_all_value(io_simVars.sim.h0);
				o_u.spectral_set_zero();
				o_v.spectral_set_zero();
			}
			else
			{
				FatalError("Unknown scenario id");
			}
		}
		else
		{
			if (io_simVars.setup.benchmark_scenario_name == "gaussian_bumps2")
			{
				SphereData tmp(o_h.sphereDataConfig);

				o_h.physical_set_all_value(io_simVars.sim.h0);

				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, io_simVars, 2.0*M_PI*0.1, M_PI/3, 20.0);
				o_h += (tmp-io_simVars.sim.h0);
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, io_simVars, 2.0*M_PI*0.6, M_PI/5.0, 80.0);
				o_h += (tmp-io_simVars.sim.h0);
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, io_simVars, 2.0*M_PI*0.8, -M_PI/4, 360.0);
				o_h += (tmp-io_simVars.sim.h0);
			}
			else if (	io_simVars.setup.benchmark_scenario_name == "geostrophic_balance"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_1"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_2"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_4"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_8"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_16"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_32"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_64"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_128"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_256"	||
					io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_512"
			)
			{
				/*
				 * geostrophic_balance / geostrophic_balance_1:
				 * Williamson test case 2 for geostrophic balance.
				 *
				 * See Williamson paper for accurate setup
				 *
				 * "geostrophic_balance_N" means that N is the multiplier for the frequency
				 * in the direction of the Latitude
				 */
				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				io_simVars.sim.coriolis_omega = 7.292e-5;
				io_simVars.sim.gravitation = 9.80616;
				io_simVars.sim.earth_radius = 6.37122e6;
				io_simVars.sim.h0 = 29400.0/io_simVars.sim.gravitation;

				//double u0 = (2.0*M_PI*io_simVars.sim.earth_radius)/(12.0*24.0*60.0*60.0);
				//double a = io_simVars.sim.earth_radius;

				//double phi0 = io_simVars.sim.h0*io_simVars.sim.gravitation;

				double freq_multiplier = 1.0;

				if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_2")
					freq_multiplier = 2.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_4")
					freq_multiplier = 4.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_8")
					freq_multiplier = 8.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_16")
					freq_multiplier = 16.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_32")
					freq_multiplier = 32.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_64")
					freq_multiplier = 64.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_128")
					freq_multiplier = 128.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_256")
					freq_multiplier = 256.0;
				else if (io_simVars.setup.benchmark_scenario_name == "geostrophic_balance_512")
					freq_multiplier = 512.0;


				/*
				 * Setup V=0
				 */
				SphereDataPhysical vg(o_h.sphereDataConfig);
				vg.physical_set_zero();

				/*
				 * Setup U=...
				 * initial velocity along longitude
				 */
				SphereDataPhysical ug(o_h.sphereDataConfig);
				ug.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = std::cos(phi*freq_multiplier);
					}
				);

				computeGeostrophicBalance(
						o_h,
						ug,
						vg,
						io_simVars,
						i_op
				);

				o_u = ug;
				o_v = vg;

				if (io_simVars.misc.sphere_use_robert_functions)
				{
					o_u.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					o_v.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);
				}
			}
			else if (io_simVars.setup.benchmark_scenario_name == "flow_over_mountain") 
			{
			        if (io_simVars.timecontrol.current_simulation_time == 0)
				  {
				    std::cout << "!!! WARNING !!!" << std::endl;
				    std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				    std::cout << "!!! WARNING !!!" << std::endl;
				  }
				
				/// Setup Williamson's parameters
				io_simVars.sim.coriolis_omega = 7.292e-5;
				io_simVars.sim.gravitation    = 9.80616;
				io_simVars.sim.earth_radius   = 6.37122e6;
				io_simVars.sim.h0             = 5600;

				const double u0 = 20.0;
					
				/*
				 * Setup V=0
				 */
				SphereDataPhysical vg(o_h.sphereDataConfig);
				vg.physical_set_zero();

				/*
				 * Setup U=...
				 * initial velocity along longitude
				 */
				SphereDataPhysical ug(o_h.sphereDataConfig);
				ug.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = u0 * std::cos(phi);
					}
				);		

				if (io_simVars.misc.sphere_use_robert_functions)
				{
					ug.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					vg.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);
				}

				
				computeGeostrophicBalance(
							  o_h,
							  ug,
							  vg,
							  io_simVars,
							  i_op
							  );
				

				o_u = ug;
				o_v = vg;

			}
			else if (io_simVars.setup.benchmark_scenario_name == "rossby_haurwitz_wave") 
			  {
			    if (io_simVars.timecontrol.current_simulation_time == 0)
			      {
				    std::cout << "!!! WARNING !!!" << std::endl;
				    std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				    std::cout << "!!! WARNING !!!" << std::endl;
			      }
				
			    /// Setup Williamson's parameters
			    io_simVars.sim.coriolis_omega = 7.292e-5;
			    io_simVars.sim.gravitation    = 9.80616;
			    io_simVars.sim.earth_radius   = 6.37122e6;
			    io_simVars.sim.h0             = 8000;

			    const double omega = 7.484e-6;
			    const double K     = omega;
			    const int    R     = 4; // wave number 4
			    
			    
			    // attempt to initialize vort and div directly
			    /*
			    SphereData vort(o_h.sphereDataConfig);
			    vort.physical_update_lambda(
			    			      [&](double lon, double phi, double &o_data)
			    			      {
			    				o_data = 2 * omega * sin(phi) 
							  - K * sin(phi) * pow(cos(phi), R) * (R*R + 3*R + 2) * cos(R*lon) ;
			    			      }
			    			      );

			    SphereData div(o_h.sphereDataConfig);
			    div.physical_set_zero();
			    */ 

			    /*
			     * Setup U=...
			     */
			    SphereDataPhysical ug(o_h.sphereDataConfig);
			    ug.physical_update_lambda(
			    			      [&](double lon, double phi, double &o_data)
			    			      {
			    				o_data = io_simVars.sim.earth_radius * omega * cos(phi) 
			    				       + io_simVars.sim.earth_radius * K * pow(cos(phi), R-1) * (R * sin(phi)*sin(phi) - cos(phi)*cos(phi)) * cos(R*lon) ;
			    			      }
			    			      );
			    

			    /*
			     * Setup V=...
			     */
			    SphereDataPhysical vg(o_h.sphereDataConfig);
			    vg.physical_update_lambda(
			    			      [&](double lon, double phi, double &o_data)
			    			      {
			    				o_data = - io_simVars.sim.earth_radius * K * R * pow(cos(phi), R-1) * sin(phi) * sin(R*lon);
			    			      }
			    			      );
			    

			    if (io_simVars.misc.sphere_use_robert_functions)
			    	{
			    		ug.physical_update_lambda_cosphi_grid(
			    						       [&](double lon, double phi, double &o_data)
			    						       {
			    							 o_data *= phi;
			    						       }
			    						       );
					
			    		vg.physical_update_lambda_cosphi_grid(
			    						       [&](double lon, double phi, double &o_data)
			    						       {
			    							 o_data *= phi;
			    						       }
			    						       );
			    	}
			    
			    
			    computeGeostrophicBalance(
						      o_h,
						      ug,
						      vg,
						      io_simVars,
						      i_op
						      );
			    
			    /*
			    o_h.physical_update_lambda(
			    			       [&](double lon, double phi, double &o_data)
			    			       {
			    				 const double A = 0.5 * omega * (2 * io_simVars.sim.coriolis_omega + omega) * cos(phi)*cos(phi)
  			    				                + 0.25 * K * K * pow(cos(phi), 2*R) * ( (R+1)*cos(phi)*cos(phi) + (2*R*R - R - 2) - 2*R*R/pow(cos(phi),2) );
			    				 const double B = 2 * (io_simVars.sim.coriolis_omega + omega) * K / ((R+1)*(R+2)) * pow(cos(phi),R) 
			    				                * ( (R*R + 2*R + 2) - (R+1)*(R+1)*cos(phi)*cos(phi) );
			    				 const double C = 0.25 * K*K * pow(cos(phi),2*R) 
  			    				                * ( (R+1)*cos(phi)*cos(phi) - (R+2) );

			    				 o_data = io_simVars.sim.h0 
			    				   + io_simVars.sim.earth_radius*io_simVars.sim.earth_radius * ( 
			    												 A
			    											       + B * cos(R*lon)  
			    											       + C * cos(2*R*lon) 
			    											       ) /io_simVars.sim.gravitation;
			    			       }
			    			       );			    
			    
			    */
     		        }
			else if (io_simVars.setup.benchmark_scenario_name == "galewsky" || io_simVars.setup.benchmark_scenario_name == "galewsky_nobump")
			{
				if (io_simVars.timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				/// Setup Galewski parameters
				io_simVars.sim.coriolis_omega = 7.292e-5;
				io_simVars.sim.gravitation = 9.80616;
				io_simVars.sim.earth_radius = 6.37122e6;
				io_simVars.sim.h0 = 10000;

				io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);

				/*
				 * Parameters from Galewsky paper setup
				 */
				double umax = 80.;
				double phi0 = M_PI/7.;
				double phi1 = 0.5*M_PI - phi0;
				double phi2 = 0.25*M_PI;
				double en = std::exp(-4.0/std::pow((phi1-phi0), 2.0));
				double alpha = 1./3.;
				double beta = 1./15.;
				double hamp = 120.;

#if 0
				/*
				 * Setup SphereDataConfig without dealiasing
				 */
				
				io_simVars.disc.res_physical[0] = 784; //2 * io_simVars.disc.res_spectral[0] + 2; //784;
				io_simVars.disc.res_physical[1] = 388; //    io_simVars.disc.res_spectral[1] + 2; //388;
				SphereDataConfig sphereConfigNoDealiasing;
				sphereConfigNoDealiasing.setupAuto(
								   io_simVars.disc.res_physical,
								   io_simVars.disc.res_spectral 
								   );
				

				/*
				 * Setup V=0
				 */
				SphereDataPhysical vg_nodealiasing(&sphereConfigNoDealiasing);
				vg_nodealiasing.physical_set_zero();
							      
				/*
				 * Setup U=...
				 * initial velocity along longitude
				 */
				SphereDataPhysical ug_nodealiasing(&sphereConfigNoDealiasing);
				ug_nodealiasing.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						if (phi >= phi1 || phi <= phi0)
							o_data = 0;
						else
							o_data = umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
					}
				);

				if (io_simVars.misc.sphere_use_robert_functions)
				{
					ug_nodealiasing.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					vg_nodealiasing.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);
				}

				SphereData vort(&sphereConfigNoDealiasing);
				SphereData div(&sphereConfigNoDealiasing);
				SphereDataPhysical ug_nodealiasing_new(&sphereConfigNoDealiasing);

				SphereOperators op_nodealiasing(&sphereConfigNoDealiasing, 
								io_simVars.sim.earth_radius);

				op_nodealiasing.robert_uv_to_vortdiv(ug_nodealiasing,
								     vg_nodealiasing,
								     vort, 
								     div);
				op_nodealiasing.robert_vortdiv_to_uv(vort,
								     div,
								     ug_nodealiasing_new,
								     vg_nodealiasing);

				char buffer[1024];
				const char* filename_template = io_simVars.misc.output_file_name_prefix.c_str();

				sprintf(buffer, 
					filename_template, 
					"u_diff_nodealiasing", 
					io_simVars.timecontrol.current_timestep_size);
				(ug_nodealiasing - ug_nodealiasing_new).physical_file_write(buffer);
								
				sprintf(buffer, 
					filename_template, 
					"vort_nodealiasing", 
					io_simVars.timecontrol.current_timestep_size);
				vort.spectrum_file_write(buffer);
				
				std::cout << (ug_nodealiasing - ug_nodealiasing_new).physical_reduce_max_abs() << std::endl;
				std::cout << (ug_nodealiasing - ug_nodealiasing_new).physical_reduce_rms() << std::endl;

#endif
				/*
				 * Setup V=0
				 */
				SphereDataPhysical vg(o_h.sphereDataConfig);
				vg.physical_set_zero();
							      
				/*
				 * Setup U=...
				 * initial velocity along longitude
				 */
				SphereDataPhysical ug(o_h.sphereDataConfig);
				ug.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						if (phi >= phi1 || phi <= phi0)
							o_data = 0;
						else
							o_data = umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
					}
				);


				if (io_simVars.misc.sphere_use_robert_functions)
				{
					ug.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					vg.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);
				}

				computeGeostrophicBalance(
						o_h,
						ug,
						vg,
						io_simVars,
						i_op
				);

				if (io_simVars.setup.benchmark_scenario_name == "galewsky")
				{
					SphereData hbump(o_h.sphereDataConfig);
					hbump.physical_update_lambda(
						[&](double lon, double phi, double &o_data)
						{
							o_data = hamp*std::cos(phi)*std::exp(-(lon/alpha)*(lon/alpha))*std::exp(-(phi2-phi)*(phi2-phi)/beta);
						}
					);

					o_h += hbump;
				}

				o_u = ug;
				o_v = vg;

			}
			else if (io_simVars.setup.benchmark_scenario_name == "flat")
			{
				o_h.physical_set_all_value(io_simVars.sim.h0);
				o_u.physical_set_all_value(0);
				o_v.physical_set_all_value(0);
			}
			else
			{
				FatalError("Benchmark not implemented");
			}
		}
	}
};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
