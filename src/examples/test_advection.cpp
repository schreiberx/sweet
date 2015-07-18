
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
#endif
#include <sweet/SimulationParameters.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Stopwatch.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>



SimulationParameters parameters;

//
// 0: staggered, up/downwinding
// 1: staggered, finite-difference
// 2: non-staggered, finite-difference
//
#define GRID_LAYOUT_AND_ADVECTION	2


class SimulationAdvection
{
public:
	DataArray<2> prog_h;
	DataArray<2> prog_u;
	DataArray<2> prog_v;

	DataArray<2> hu;
	DataArray<2> hv;

	Operators2D op;

	TimesteppingRK timestepping;


public:
	SimulationAdvection()	:
		prog_h(parameters.res),
		prog_u(parameters.res),
		prog_v(parameters.res),

		hu(parameters.res),
		hv(parameters.res),

		op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs)
	{
		reset();
	}


	void reset()
	{
		parameters.status_timestep_nr = 0;

		prog_u.setAll(parameters.bogus_var0);
		prog_v.setAll(parameters.bogus_var1);

		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
				double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

				prog_h.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
			}
		}
#if 0
		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				if (parameters.bogus_var2 == 0 || parameters.bogus_var2 == 1)
				{
					double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];
				}
				else
				{
					double x = (((double)i)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j)/(double)parameters.res[1])*parameters.sim_domain_size[1];
				}

				prog_u.set(j, i, parameters.bogus_var0+sin(2.0*M_PIl*x/(double)parameters.sim_domain_size[0])*0.01);
			}
		}
#endif
	}



	void p_run_euler_timestep_update(
			const DataArray<2> &i_h,	///< prognostic variables
			const DataArray<2> &i_u,	///< prognostic variables
			const DataArray<2> &i_v,	///< prognostic variables

			DataArray<2> &o_h_t,	///< time updates
			DataArray<2> &o_u_t,	///< time updates
			DataArray<2> &o_v_t,	///< time updates

			double &o_dt,			///< time step restriction
			double i_fixed_dt = 0	///< if this value is not equal to 0, use this time step size instead of computing one
	)
	{
		if (parameters.bogus_var2 == 0)
		{
#if SWEET_USE_SPECTRAL_SPACE
			static bool output_given = false;
			if (!output_given)
			{
				std::cout << "WARNING: upwinding in spectral space not working" << std::endl;
				output_given = true;
			}
#endif
			o_h_t =
				(
					(
						// u is positive
						op.shift_right(i_h)*i_u.return_value_if_positive()	// inflow
						-i_h*op.shift_left(i_u.return_value_if_positive())	// outflow

						// u is negative
						+(i_h*i_u.return_value_if_negative())					// outflow
						-op.shift_left(i_h*i_u.return_value_if_negative())	// inflow
					)*(1.0/parameters.sim_cell_size[0])				// here we see a finite-difference-like formulation
					+
					(
						// v is positive
						op.shift_up(i_h)*i_v.return_value_if_positive()		// inflow
						-i_h*op.shift_down(i_v.return_value_if_positive())	// outflow

						// v is negative
						+(i_h*i_v.return_value_if_negative())					// outflow
						-op.shift_down(i_h*i_v.return_value_if_negative())	// inflow
					)*(1.0/parameters.sim_cell_size[1])
				);
		}
		else if (parameters.bogus_var2 == 1)
		{
			//             |                       |                       |
			// --v---------|-----------v-----------|-----------v-----------|
			//   h-1       u0          h0          u1          h1          u2

			// staggered
			o_h_t = -(
					op.diff_f_x(op.avg_b_x(i_h)*i_u) +
					op.diff_f_y(op.avg_b_y(i_h)*i_v)
				);

		}
		else  if (parameters.bogus_var2 == 2)
		{
			// non-staggered
			o_h_t = -(
					op.diff_c_x(i_h*i_u) +
					op.diff_c_y(i_h*i_v)
				);
		}
		else  if (parameters.bogus_var2 == 3)
		{
			o_h_t.setAll(0);
		}
		else
		{
			std::cerr << "Advection type not specified, use -c option [0: up/downwinding, 1: staggered, 2: non-staggered]" << std::endl;
		}


		if (i_fixed_dt > 0)
		{
			o_dt = i_fixed_dt;
		}
		else
		{
			if (i_fixed_dt < 0)
				o_dt = -i_fixed_dt;
			else
				o_dt = parameters.sim_CFL*std::min(parameters.sim_cell_size[0]/i_u.reduce_maxAbs(), parameters.sim_cell_size[1]/i_v.reduce_maxAbs());
		}

		o_u_t.setAll(0);
		o_v_t.setAll(0);

		parameters.status_timestep_nr++;
	}


	void run()
	{
	}


	void run_timestep()
	{
		double dt;

		timestepping.run_rk_timestep(
				this,
				&SimulationAdvection::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_h, prog_u, prog_v,
				dt,
				parameters.timestepping_timestep_size,
				parameters.timestepping_runge_kutta_order
			);


		// provide information to parameters
		parameters.status_simulation_timestep_size = dt;
		parameters.status_simulation_time += dt;
		parameters.status_timestep_nr++;
	}



	bool should_quit()
	{
		return false;
	}



	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (parameters.run_simulation)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}



	void vis_get_vis_data_array(
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		*o_dataArray = &prog_h;
		*o_aspect_ratio = parameters.sim_domain_size[1] / parameters.sim_domain_size[0];
	}

	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				parameters.status_simulation_time,
				parameters.status_simulation_time/(60.0*60.0*24.0),
				parameters.status_timestep_nr,
				parameters.status_simulation_timestep_size,
				parameters.diagnostics_mass, parameters.diagnostics_energy, parameters.diagnostics_potential_entrophy);
		return title_string;
	}

	void vis_pause()
	{
		parameters.run_simulation = !parameters.run_simulation;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			parameters.vis_id++;
			break;

		case 'V':
			parameters.vis_id--;
			break;
		}
	}

	bool instability_detected()
	{
		return !(	prog_h.reduce_all_finite() &&
					prog_u.reduce_all_finite() &&
					prog_v.reduce_all_finite()
				);
	}
};



int main(
		int i_argc,
		char *i_argv[]
)
{
	parameters.setup(i_argc, i_argv);

	double u = parameters.bogus_var0;
	double v = parameters.bogus_var1;

	double total_speed;
	double turnaround_time;
	if (u == 0 && v == 0)
	{
		std::cerr << "Both velocity components are zero, EXIT" << std::endl;
		exit(1);
	}

	if (u != 0 && v == 0)
	{
		total_speed = u;
		turnaround_time = parameters.sim_domain_size[0]/u;
	}
	else if (u == 0 && v != 0)
	{
		total_speed = v;
		turnaround_time = parameters.sim_domain_size[1]/v;
	}
	else
	{
		total_speed = v;
		if (std::abs(parameters.sim_domain_size[1]/parameters.sim_domain_size[0]-v/u) > 0.000000001)
		{
			std::cerr << "ratio of domain sizes and speed have to be similar" << std::endl;
			exit(1);
		}

		total_speed = std::sqrt(u*u+v*v);
		double diagonal = std::sqrt(parameters.sim_domain_size[0]*parameters.sim_domain_size[0] + parameters.sim_domain_size[1]*parameters.sim_domain_size[1]);
		turnaround_time = diagonal/total_speed;
	}

	if (parameters.verbosity > 1)
	{
		std::cout << "Turnaround time: " << turnaround_time << std::endl;
		std::cout << "Total speed: " << total_speed << std::endl;
	}



#if SWEET_GUI
	if (parameters.gui_enabled)
	{
		SimulationAdvection *simulationAdvection = new SimulationAdvection;
		VisSweet<SimulationAdvection> visSweet(simulationAdvection);
		delete simulationAdvection;
	}
	else
#endif

	{
		/*
		 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
		 */
		double prev_error = 0;

		std::size_t res_x = parameters.res[0];
		std::size_t res_y = parameters.res[1];

		std::size_t max_res = 256;

		if (res_x > max_res || res_y > max_res)
			max_res = std::max(res_x, res_y);

		for (; res_x <= max_res && res_y <= max_res; res_x *= 2, res_y *= 2)
		{
			std::cout << "*************************************************************" << std::endl;
			std::cout << "Testing advection with resolution " << res_x << " x " << res_y << std::endl;
			std::cout << "*************************************************************" << std::endl;
			std::size_t res[2] = {res_x, res_y};

			parameters.res[0] = res[0];
			parameters.res[1] = res[1];
			parameters.reset();

			SimulationAdvection *simulationAdvection = new SimulationAdvection;

			Stopwatch time;
			time.reset();

			while(true)
			{
				if (parameters.verbosity > 0)
					std::cout << "time: " << parameters.status_simulation_time << std::endl;

				simulationAdvection->run_timestep();

				if (simulationAdvection->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}

				if (turnaround_time < parameters.status_simulation_time)
				{
					DataArray<2> benchmark_h(parameters.res);

					for (std::size_t j = 0; j < parameters.res[1]; j++)
					{
						for (std::size_t i = 0; i < parameters.res[0]; i++)
						{
							double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
							double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

							benchmark_h.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
						}
					}

					double error = (simulationAdvection->prog_h-benchmark_h).reduce_rms_quad();
					std::cout << "RMS error in height: " << error << std::endl;

					if (prev_error != 0)
					{
						double conv_rate = prev_error / error;
						std::cout << "          Convergence rate: " << conv_rate << std::endl;

						if (std::abs(conv_rate-4.0) > 0.5)
						{
							std::cerr << "Convergence rate threshold (0.5) exceeded" << std::endl;
							exit(1);
						}
					}

					prev_error = error;
					break;
				}
			}

			time.stop();

			double seconds = time();

			std::cout << "Simulation time: " << seconds << " seconds" << std::endl;
			std::cout << "Time per time step: " << seconds/(double)parameters.status_timestep_nr << " sec/ts" << std::endl;

			delete simulationAdvection;
		}
	}


	return 1;
}
