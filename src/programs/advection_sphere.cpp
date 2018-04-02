/*
 * lagrangian_sphere_test.cpp
 *
 *  Created on: 3 Dec 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SWEET_GUI
#define SWEET_GUI 1
#endif

#include "../include/sweet/sphere/SphereData.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <benchmarks_sphere/SphereBenchmarksCombined.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/Convert_SphereData_To_PlaneData.hpp>
#include <sweet/Convert_SphereDataPhysical_To_PlaneData.hpp>

#include "advection_sphere/Adv_Sphere_TimeSteppers.hpp"



// Sphere data config
SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif

SimulationVariables simVars;


class SimulationInstance
{
public:
	SphereData prog_h;
	SphereData prog_h0;	// at t0
	SphereData prog_vort, prog_div;

	Adv_Sphere_TimeSteppers timeSteppers;


	SphereOperators op;

#if SWEET_GUI
	PlaneData viz_plane_data;

	int render_primitive_id = 1;
#endif


public:
	SimulationInstance()	:
		prog_h(sphereDataConfig),
		prog_h0(sphereDataConfig),

		prog_vort(sphereDataConfig),
		prog_div(sphereDataConfig),

		op(sphereDataConfig, simVars.sim.earth_radius)

#if SWEET_GUI
		,
		viz_plane_data(planeDataConfig)
#endif
	{
		reset();
	}



	void reset()
	{
		simVars.timecontrol.current_timestep_nr = 0;

		SphereData tmp_vort(sphereDataConfig);
		SphereData tmp_div(sphereDataConfig);

		SphereBenchmarksCombined::setupInitialConditions(prog_h, prog_vort, prog_div, simVars, op);

		timeSteppers.setup(simVars.disc.timestepping_method, op, simVars);

		simVars.outputConfig();
	}



	void run_timestep()
	{

		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_h, prog_vort, prog_div,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);

#if 0
		double dt = simVars.timecontrol.current_timestep_size;

		op.robert_vortdiv_to_uv(prog_vort, prog_div, diag_u, diag_v);


#if 0
		// velocities at t
		SphereDataPhysical* vel[2] = {&diag_u, &diag_v};

		// velocities at t-1
		SphereDataPhysical* vel_prev[2] = {&diag_u_prev, &diag_v_prev};

		// OUTPUT: position of departure points at t
		ScalarDataArray posx_d(sphereDataConfig->physical_array_data_number_of_elements);
		ScalarDataArray posy_d(sphereDataConfig->physical_array_data_number_of_elements);
		ScalarDataArray* output_pos_departure[2] = {&posx_d, &posy_d};

		//return;
		semiLagrangian.compute_departure_points_settls(
				vel_prev,
				vel,
				input_pos_arrival,
				dt,
				output_pos_departure
		);

		diag_u_prev = diag_u;
		diag_v_prev = diag_v;
#endif

		return;
		SphereData new_prog_h(sphereDataConfig);
		sampler2D.bicubic_scalar(
				prog_h,
#if 0
				posx_d,
				posy_d,
#else
				posx_a,
				posy_a,
#endif
				new_prog_h
		);

		prog_h = new_prog_h;


		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

#endif
	}


	void compute_error()
	{
#if 0
		double t = simVars.timecontrol.current_simulation_time;

		SphereData prog_testh(sphereDataConfig);
		prog_testh.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double x = (((double)i)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
				double y = (((double)j)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

				x -= param_velocity_u*t;
				y -= param_velocity_v*t;

				while (x < 0)
					x += simVars.sim.domain_size[0];

				while (y < 0)
					y += simVars.sim.domain_size[1];

				x = std::fmod(x, simVars.sim.domain_size[0]);
				y = std::fmod(y, simVars.sim.domain_size[1]);

				io_data = SWESphereBenchmarks::return_h(simVars, x, y);
			}
		);

		std::cout << "Lmax Error: " << (prog_h-prog_testh).reduce_maxAbs() << std::endl;
#endif
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
			int *o_render_primitive_id,
			void **o_bogus_data
	)
	{
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		int id = simVars.misc.vis_id % 5;
		switch (id)
		{
		case 0:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_h, planeDataConfig);
			break;

		case 1:
			{
				SphereDataPhysical u(sphereDataConfig);
				SphereDataPhysical v(sphereDataConfig);

//				op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);
				op.vortdiv_to_uv(prog_vort, prog_div, u, v);
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(u, planeDataConfig);
			}
			break;

		case 2:
			{
				SphereDataPhysical u(sphereDataConfig);
				SphereDataPhysical v(sphereDataConfig);

//				op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);
				op.vortdiv_to_uv(prog_vort, prog_div, u, v);
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(v, planeDataConfig);
			}
			break;

		case 3:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_vort, planeDataConfig);
			break;

		case 4:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_div, planeDataConfig);
			break;
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
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

		case 'b':
			render_primitive_id = (render_primitive_id + 1) % 2;
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
		return -1;
	}

	if (simVars.timecontrol.current_timestep_size < 0)
		FatalError("Timestep size not set");

	sphereDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	SimulationInstance *simulationSWE = new SimulationInstance;

#if SWEET_GUI

	if (simVars.misc.gui_enabled)
	{
		planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical);

		VisSweet<SimulationInstance> visSweet(simulationSWE);
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
