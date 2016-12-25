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
			SimulationVariables &i_simVars,
			SphereOperators &i_op
	)
	{
		BenchmarkGalewsky benchmarkGalewsky(i_simVars);

		if (i_simVars.setup.benchmark_scenario_id <= 0)
		{
			std::cout << std::endl;
			std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
			std::cout << "Available benchmark scenarios:" << std::endl;
			std::cout << "	1: Galweski" << std::endl;
			std::cout << "	2: Use Gaussian bump initial conditions (0, pi/3)" << std::endl;
			std::cout << "	3: Use Gaussian bump initial conditions (pi/3, pi/3)" << std::endl;
			std::cout << "	4: Use geostrophic balance test case" << std::endl;
			std::cout << std::endl;
			FatalError("Benchmark scenario not selected");
		}


		if (i_simVars.setup.benchmark_scenario_id == 0)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, i_simVars, M_PI/3.0, M_PI/3.0);
		}
		else if (i_simVars.setup.benchmark_scenario_id == 1)
		{
			benchmarkGalewsky.setup_initial_h(o_h);
//			o_h.spat_set_zero();
			benchmarkGalewsky.setup_initial_h_add_bump(o_h);

			benchmarkGalewsky.setup_initial_u(o_u);
			benchmarkGalewsky.setup_initial_v(o_v);
		}
		else if (i_simVars.setup.benchmark_scenario_id == 2 || i_simVars.setup.benchmark_scenario_id == 3)
		{
#if 0
			if (i_simVars.timecontrol.current_timestep_size <= 0)
			{
				std::cout << "Timestep size not positive" << std::endl;
				assert(false);
				exit(1);
			}
#endif

			if (i_simVars.setup.benchmark_scenario_id == 2)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, i_simVars, 0, M_PI/3.0);
			}
			else if (i_simVars.setup.benchmark_scenario_id == 3)
			{
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, i_simVars, M_PI/3.0, M_PI/3.0);
				//setup_initial_conditions_gaussian(M_PI, M_PI);
//				setup_initial_conditions_gaussian(M_PI*0.5, M_PI*0.5);
//				setup_initial_conditions_gaussian(-M_PI/3.0);
			}
		}
		else if (i_simVars.setup.benchmark_scenario_id == 4)
		{
			double inv_r = 1.0/i_simVars.sim.earth_radius;

			o_v.spectral_set_zero();

			o_h.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						io_data = i_simVars.sim.earth_radius*i_simVars.sim.coriolis_omega*std::cos(i_lat)*std::cos(i_lat)/i_simVars.sim.gravitation;
					}
			);

			if (i_simVars.misc.sphere_use_robert_functions)
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

			if (i_simVars.timecontrol.current_simulation_time == 0)
			{
				double h_max_error, u_max_error, v_max_error;

				if (i_simVars.misc.sphere_use_robert_functions)
				{
					h_max_error = (-inv_r*i_op.robert_div_lon(o_u) -inv_r*i_op.robert_div_lat(o_v)).physical_reduce_max_abs();
					u_max_error = (-inv_r*i_op.robert_grad_lon(o_h*i_simVars.sim.gravitation) + 2.0*i_simVars.sim.coriolis_omega*i_op.mu(o_v)).physical_reduce_max_abs();
					v_max_error = (-inv_r*i_op.robert_grad_lat(o_h*i_simVars.sim.gravitation) - 2.0*i_simVars.sim.coriolis_omega*i_op.mu(o_u)).physical_reduce_max_abs();
				}
				else
				{
					h_max_error = (-inv_r*i_op.div_lon(o_u) -inv_r*i_op.div_lat(o_v)).physical_reduce_max_abs();
					u_max_error = (-inv_r*i_op.grad_lon(o_h*i_simVars.sim.gravitation) + 2.0*i_simVars.sim.coriolis_omega*i_op.mu(o_v)).physical_reduce_max_abs();
					v_max_error = (-inv_r*i_op.grad_lat(o_h*i_simVars.sim.gravitation) - 2.0*i_simVars.sim.coriolis_omega*i_op.mu(o_u)).physical_reduce_max_abs();
				}
				std::cout << "h_max_error for geostrophic balance case: " << h_max_error << std::endl;
				std::cout << "u_max_error for geostrophic balance case: " << u_max_error << std::endl;
				std::cout << "v_max_error for geostrophic balance case: " << v_max_error << std::endl;
			}
		}
		else
		{
			FatalError("Unknown scenario id");
		}
	}

};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
