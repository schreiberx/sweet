/*
 * lagrangian_test.cpp
 *
 *  Created on: 3 Dec 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <schreiberx@gmail.com>
 */


#include "../include/sweet/plane/PlaneData.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <array>


// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;


SimulationVariables simVars;

double param_velocity_u;
double param_velocity_v;


class SimulationSWE
{
public:
	PlaneData prog_h;

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
		prog_h(planeDataConfig),

		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		prog_u_prev(planeDataConfig),
		prog_v_prev(planeDataConfig),

		posx_a(planeDataConfig->physical_array_data_number_of_elements),
		posy_a(planeDataConfig->physical_array_data_number_of_elements),

		h_t(planeDataConfig),

		op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
	{
		reset();
	}



	void reset()
	{
		simVars.timecontrol.current_timestep_nr = 0;
#if 0
		prog_h.set_all(simVars.sim.h0);

		if (std::isinf(simVars.bogus.var[0]) != 0)
		{
			prog_u.set_all(0);
			prog_v.set_all(0);
		}
		else
		{
			prog_u.set_all(simVars.bogus.var[0]);
			prog_v.set_all(simVars.bogus.var[1]);
		}
#endif


		prog_h.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double x = (((double)i)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
				double y = (((double)j)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

				io_data = SWEPlaneBenchmarks::return_h(simVars, x, y);
			}
		);

//		prog_u += 4.0;
//		prog_v += 4.0;

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

				io_data = ((double)i)*((double)simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0]);
			}
		);
		posy_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				//int i = idx % planeDataConfig->physical_data_size[0];
				int j = idx / planeDataConfig->physical_data_size[0];

				io_data = ((double)j)*((double)simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1]);
			}
		);

		sampler2D.setup(simVars.sim.domain_size, planeDataConfig);

		//PXT- This just calls sampler2D.setup, so any reason for having it?
		semiLagrangian.setup(simVars.sim.domain_size, planeDataConfig);
	}



	void run_timestep()
	{
		double dt = std::min(
				std::abs((double)planeDataConfig->physical_res[0]/param_velocity_u),
				std::abs((double)planeDataConfig->physical_res[1]/param_velocity_v)
			)*simVars.sim.CFL;

		if ((param_velocity_u == 0) && (param_velocity_v == 0))
			dt = 0;

		simVars.timecontrol.current_timestep_size = dt;

		// update velocities
		prog_u = param_velocity_u;
		prog_v = param_velocity_v;

		// velocities at t
		PlaneData* vel[2] = {&prog_u, &prog_v};
		// velocities at t-1
		PlaneData* vel_prev[2] = {&prog_u_prev, &prog_v_prev};

		// position of departure points at t
		ScalarDataArray posx_d(planeDataConfig->physical_array_data_number_of_elements);
		ScalarDataArray posy_d(planeDataConfig->physical_array_data_number_of_elements);
		ScalarDataArray* output_pos_departure[2] = {&posx_d, &posy_d};

#if 0
		*output_pos_departure[0] = *input_pos_arrival[0];
		*output_pos_departure[1] = *input_pos_arrival[1];
#else
		semiLagrangian.compute_departure_points_settls(
				vel_prev,
				vel,
				input_pos_arrival,
				dt,
				output_pos_departure
		);
#endif

		prog_u_prev = prog_u;
		prog_v_prev = prog_v;

		PlaneData new_prog_h(planeDataConfig);
		sampler2D.bicubic_scalar(
				prog_h,
				posx_d,
				posy_d,
				new_prog_h
		);

		prog_h = new_prog_h;

#if 0
		PlaneData *x[2] = {&data_x, &data_y};

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position
		for (int j = 0; j < simVars.disc.res_physical[1]; j++)
		{
			for (int i = 0; i < simVars.disc.res_physical[0]; i++)
			{
				x[0]->set(j, i, ((double)i)*(simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0])+simVars.timecontrol.current_timestep_nr);
				x[1]->set(j, i, ((double)j)*(simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1])+simVars.timecontrol.current_timestep_nr*2);
			}
		}

		sampler2D.bicubic_scalar(
				prog_h,	///< input scalar field
				x,		///< sampling positions
				prog_u	///< output field
		);

		// interpolate data from x/y coordinates
//		interpol_values = op_sample.bilinear(prog_h, x, y);
#endif

		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;
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
	}


	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data
	)
	{
		switch (simVars.misc.vis_id)
		{
		case 0:
			*o_dataArray = &prog_h;
			break;

		case 1:
			*o_dataArray = &prog_u;
			break;

		case 2:
			*o_dataArray = &prog_v;
			break;
		}
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
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

	if (std::isinf(simVars.bogus.var[0]) != 0 || std::isinf(simVars.bogus.var[1]) != 0)
	{
		std::cout << "Both velocities have to be set, see parameters --velocity-u, --velocity-v" << std::endl;
		return -1;
	}

	param_velocity_u = simVars.bogus.var[0];
	param_velocity_v = simVars.bogus.var[1];

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

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
