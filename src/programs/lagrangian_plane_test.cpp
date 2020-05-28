/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SWEET_GUI
#define SWEET_GUI 1
#endif

#include "../include/sweet/plane/PlaneData.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
//#include <benchmarks_plane/SWE_bench_PlaneBenchmarks_DEPRECATED.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarksCombined.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>


// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;


SimulationVariables simVars;

double param_velocity_u;
double param_velocity_v;


class SimulationSWE
{
public:
	PlaneData prog_h_pert;

	PlaneData prog_h0_pert;	// at t0

	PlaneData prog_u, prog_v;
	PlaneData prog_u_prev, prog_v_prev;

	ScalarDataArray posx_a, posy_a;
	ScalarDataArray *input_pos_arrival[2];

	PlaneData h_t;

	PlaneOperators op;

	PlaneDataSampler sampler2D;
	PlaneDataSemiLagrangian semiLagrangian;



public:
	SimulationSWE()	:
		prog_h_pert(planeDataConfig),
		prog_h0_pert(planeDataConfig),

		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		prog_u_prev(planeDataConfig),
		prog_v_prev(planeDataConfig),

		posx_a(planeDataConfig->physical_array_data_number_of_elements),
		posy_a(planeDataConfig->physical_array_data_number_of_elements),

		h_t(planeDataConfig),

		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)
	{
		reset();
	}



	void reset()
	{
		simVars.timecontrol.current_timestep_nr = 0;

		SWEPlaneBenchmarksCombined swePlaneBenchmarks;
		swePlaneBenchmarks.setupInitialConditions(prog_h_pert, prog_u, prog_v, simVars, op);

#if 0
		prog_h_pert.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double x = (((double)i)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
				double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

				io_data = SWEPlaneBenchmarks_DEPRECATED::return_h(simVars, x, y);
			}
		);
#endif

		prog_h0_pert = prog_h_pert;

		prog_u = param_velocity_u;
		prog_v = param_velocity_v;

		prog_u_prev = prog_u;
		prog_v_prev = prog_v;

		input_pos_arrival[0] = &posx_a;
		input_pos_arrival[1] = &posy_a;

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position

		posx_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
			{
				int i = idx % planeDataConfig->physical_data_size[0];
				//int j = idx / planeDataConfig->physical_data_size[0];

				io_data = ((double)i)*((double)simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0]);
			}
		);
		posy_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				//int i = idx % planeDataConfig->physical_data_size[0];
				int j = idx / planeDataConfig->physical_data_size[0];

				io_data = ((double)j)*((double)simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1]);
			}
		);

		sampler2D.setup(simVars.sim.plane_domain_size, planeDataConfig);

		//PXT- This just calls sampler2D.setup, so any reason for having it?
		semiLagrangian.setup(simVars.sim.plane_domain_size, planeDataConfig);
	}



	void run_timestep()
	{
		double dt = simVars.timecontrol.current_timestep_size;

		// update velocities
		prog_u = param_velocity_u;
		prog_v = param_velocity_v;

		// position of departure points at t
		ScalarDataArray posx_d(planeDataConfig->physical_array_data_number_of_elements);
		ScalarDataArray posy_d(planeDataConfig->physical_array_data_number_of_elements);

		semiLagrangian.semi_lag_departure_points_settls(
				prog_u_prev, prog_v_prev,
				prog_u, prog_v,
				posx_a, posy_a,
				dt,
				posx_d, posy_d,
				simVars.sim.plane_domain_size,
				nullptr,
				simVars.disc.timestepping_order,

				simVars.disc.semi_lagrangian_max_iterations,
				simVars.disc.semi_lagrangian_convergence_threshold
		);

		prog_u_prev = prog_u;
		prog_v_prev = prog_v;

		PlaneData new_prog_h(planeDataConfig);
		sampler2D.bicubic_scalar(
				prog_h_pert,
				posx_d,
				posy_d,
				new_prog_h
		);

		prog_h_pert = new_prog_h;


		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;
	}


	void compute_error()
	{

#if 0
		double t = simVars.timecontrol.current_simulation_time;

		PlaneData prog_testh(planeDataConfig);
		prog_testh.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double x = (((double)i)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
				double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

				x -= param_velocity_u*t;
				y -= param_velocity_v*t;

				while (x < 0)
					x += simVars.sim.plane_domain_size[0];

				while (y < 0)
					y += simVars.sim.plane_domain_size[1];

				x = std::fmod(x, simVars.sim.plane_domain_size[0]);
				y = std::fmod(y, simVars.sim.plane_domain_size[1]);

				io_data = SWEPlaneBenchmarks_DEPRECATED::return_h(simVars, x, y);
			}
		);
#endif
		std::cerr << "TODO: This was removed due to removal of the benchmark_id parameter. A new benchmark must be generated for this" << std::endl;

		std::cout << "Lmax Error: " << (prog_h_pert-prog_h0_pert).reduce_maxAbs() << std::endl;
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
		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();

		compute_error();
	}


	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		switch (simVars.misc.vis_id)
		{
		case 0:
			*o_dataArray = &prog_h_pert;
			break;

		case 1:
			*o_dataArray = &prog_u;
			break;

		case 2:
			*o_dataArray = &prog_v;
			break;
		}
		*o_aspect_ratio = simVars.sim.plane_domain_size[1] / simVars.sim.plane_domain_size[0];
	}



	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Timestep: %i, timestep size: %e", simVars.timecontrol.current_timestep_nr, simVars.timecontrol.current_timestep_size);
		return title_string;
	}

	void vis_pause()
	{
		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			simVars.misc.vis_id++;
			break;

		case 'V':
			simVars.misc.vis_id--;
			break;
		}
	}
};



int main(int i_argc, char *i_argv[])
{
	const char *bogus_var_names[] = {
			"velocity-u",
			"velocity-v",
			nullptr
	};

	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Program-specific options:" << std::endl;
		std::cout << "	--velocity-u [advection velocity u]" << std::endl;
		std::cout << "	--velocity-v [advection velocity v]" << std::endl;
		return -1;
	}

	if (simVars.bogus.var[0] != "")
		param_velocity_u = atof(simVars.bogus.var[0].c_str());

	if (simVars.bogus.var[1] != "")
		param_velocity_v = atof(simVars.bogus.var[1].c_str());

	if (param_velocity_u == 0 && param_velocity_v == 0)
	{
		std::cout << "At least one velocity has to be set, see parameters --velocity-u, --velocity-v" << std::endl;
		return -1;
	}

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

	SimulationSWE *simulationSWE = new SimulationSWE;

#if SWEET_GUI
	if (simVars.misc.gui_enabled)
	{
		VisSweet<SimulationSWE> visSweet(simulationSWE);
	}
	else
#endif
	{
		simulationSWE->reset();
		while (!simulationSWE->should_quit())
		{
			simulationSWE->run_timestep();

			if (simVars.misc.verbosity > 2)
				std::cout << simVars.timecontrol.current_simulation_time << std::endl;

			if (simVars.timecontrol.max_simulation_time != -1)
				if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
					break;

			if (simVars.timecontrol.max_timesteps_nr != -1)
				if (simVars.timecontrol.current_timestep_nr > simVars.timecontrol.max_timesteps_nr)
					break;
		}
	}

	delete simulationSWE;

	return 0;
}
