/*
 * BenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/BenchmarkGalewsky.hpp>
#include <benchmarks_sphere/BenchmarkGaussianDam.hpp>

class SphereBenchmarksCombined
{
public:
	static
	void setupInitialConditions(
			SphereData &o_h,
			SphereData &o_u,
			SphereData &o_v,
			SimulationVariables &io_simVars,
			SphereOperators &i_op
	)
	{
		BenchmarkGalewsky benchmarkGalewsky(io_simVars);

		if (io_simVars.setup.benchmark_scenario_id <= 0)
		{
			std::cout << std::endl;
			std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
			std::cout << "Available benchmark scenarios:" << std::endl;
			std::cout << "	1: Galweski" << std::endl;
			std::cout << "	2: Use Gaussian bump initial conditions (0, pi/3)" << std::endl;
			std::cout << "	3: Use Gaussian bump initial conditions (pi/3, pi/3)" << std::endl;
			std::cout << "	4: Use geostrophic balance test case" << std::endl;
			std::cout << "	5: Williamson test benchmark 1" << std::endl;
			std::cout << std::endl;
			FatalError("Benchmark scenario not selected");
		}


		if (io_simVars.setup.benchmark_scenario_id == 0)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 1)
		{
			benchmarkGalewsky.setup_initial_h(o_h);
//			o_h.spat_set_zero();
			benchmarkGalewsky.setup_initial_h_add_bump(o_h);

			benchmarkGalewsky.setup_initial_u(o_u);
			benchmarkGalewsky.setup_initial_v(o_v);

			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;


			/// Setup Galewski parameters
			io_simVars.sim.coriolis_omega = 7.292e-5;
			io_simVars.sim.gravitation = 9.80616;
			io_simVars.sim.earth_radius = 6.37122e6;
			io_simVars.sim.h0 = 10000.0;

			io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 2)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, 0, M_PI/3.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 3)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 4)
		{
			double inv_r = 1.0/io_simVars.sim.earth_radius;

			o_v.spectral_set_zero();

			o_h.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						io_data = io_simVars.sim.earth_radius*io_simVars.sim.coriolis_omega*std::cos(i_lat)*std::cos(i_lat)/io_simVars.sim.gravitation;
					}
			);

			if (io_simVars.misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							io_data = std::cos(i_lat)*std::cos(i_lat);
						}
				);
			}
			else
			{
				o_u.physical_update_lambda(
						[&](double i_lon, double i_lat, double &io_data)
						{
							//io_data = simVars.sim.earth_radius*2.0*simVars.sim.coriolis_omega*std::cos(i_lat)/simVars.sim.gravitation;
							io_data = std::cos(i_lat);
						}
				);
			}

			if (io_simVars.timecontrol.current_simulation_time == 0)
			{
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
		}
		else if (io_simVars.setup.benchmark_scenario_id == 5)
		{
			/**
			 * See Williamson test case, eq. (75), (76)
			 */

			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;

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
						double i_phi = i_lat;
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
						double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
									std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
							);
					}
				);
			}

			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Storing advection in output velocities !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;


			io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);

			std::cout << "advection_rotation_angle: " << io_simVars.setup.advection_rotation_angle << std::endl;		}
		else
		{
			FatalError("Unknown scenario id");
		}
	}

};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
