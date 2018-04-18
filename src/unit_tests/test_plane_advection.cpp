/*
 * test_plane_advection.cpp
 *
 *  Created on: 4th April 2018
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */


#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#include "../include/sweet/plane/PlaneData.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <benchmarks_plane/SWEBenchmarksCombined.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#include "test_plane_advection/Adv_Plane_TimeSteppers.hpp"



// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

SimulationVariables simVars;

double param_velocity_u;
double param_velocity_v;



class SimulationInstance
{
public:
	PlaneData prog_h;
	PlaneData prog_h0;	// at t0
	PlaneData prog_u, prog_v;

	Adv_Plane_TimeSteppers timeSteppers;

	double center_x = 0.5;
	double center_y = 0.5;
	double exp_fac = 50.0;

	PlaneOperators op;

#if SWEET_GUI
	PlaneData viz_plane_data;

	int render_primitive_id = 0;
#endif


	double max_error_h0 = -1;
	double rms_error_h0 = -1;


public:
	SimulationInstance()	:
		prog_h(planeDataConfig),
		prog_h0(planeDataConfig),

		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		op(planeDataConfig, simVars.sim.domain_size)

#if SWEET_GUI
		,
		viz_plane_data(planeDataConfig)
#endif
	{
		reset();
	}


	~SimulationInstance()
	{
	}



	static
	double gaussianValue(
			double i_center_x, double i_center_y,
			double i_x, double i_y,
			double i_exp_fac
	)
	{
		double sx = simVars.sim.domain_size[0];
		double sy = simVars.sim.domain_size[1];

		// Gaussian
		double dx = i_x-i_center_x*sx;
		double dy = i_y-i_center_y*sy;

		if (dx > 0.5*simVars.sim.domain_size[0])
			dx -= simVars.sim.domain_size[0];
		else if (dx < -0.5*simVars.sim.domain_size[0])
			dx += simVars.sim.domain_size[0];

		if (dy > 0.5*simVars.sim.domain_size[1])
			dy -= simVars.sim.domain_size[1];
		else if (dy < -0.5*simVars.sim.domain_size[1])
			dy += simVars.sim.domain_size[1];

		dx /= sx*simVars.setup.radius_scale;
		dy /= sy*simVars.setup.radius_scale;

		return std::exp(-i_exp_fac*(dx*dx + dy*dy));
	}


	void reset()
	{
		simVars.reset();

		prog_h.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				double x = (double)i*(simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0]);
				double y = (double)j*(simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1]);

				io_data = gaussianValue(center_x, center_y, x, y, exp_fac);
			}
		);

		prog_h.spectral_zeroAliasingModes();
		prog_h0 = prog_h;

		prog_u = param_velocity_u;
		prog_v = param_velocity_v;

		// setup planeDataconfig instance again
		planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

		timeSteppers.setup(simVars.disc.timestepping_method, op, simVars);

		if (simVars.misc.verbosity > 2)
			simVars.outputConfig();
	}



	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_h, prog_u, prog_v,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);

		double dt = simVars.timecontrol.current_timestep_size;

		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.misc.verbosity > 2)
			std::cout << simVars.timecontrol.current_timestep_nr << ": " << simVars.timecontrol.current_simulation_time/(60*60*24.0) << std::endl;

		max_error_h0 = (prog_h-prog_h0).reduce_maxAbs();
		rms_error_h0 = (prog_h-prog_h0).reduce_rms();
	}



	void compute_error()
	{
#if 0
		double t = simVars.timecontrol.current_simulation_time;

		PlaneData prog_testh(planeDataConfig);
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

				io_data = SWEPlaneBenchmarks::return_h(simVars, x, y);
			}
		);

		std::cout << "Lmax Error: " << (prog_h-prog_testh).reduce_maxAbs() << std::endl;
#endif
	}



	bool should_quit()
	{
		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		double diff = std::abs(simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time);

		if (	simVars.timecontrol.max_simulation_time != -1 &&
				(
						simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time	||
						diff/simVars.timecontrol.max_simulation_time < 1e-11	// avoid numerical issues in time stepping if current time step is 1e-14 smaller than max time step
				)
			)
			return true;

		return false;
	}

#if SWEET_GUI
	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
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
//		*o_bogus_data = planeDataConfig;

		int id = simVars.misc.vis_id % 3;
		switch (id)
		{
		case 0:
			viz_plane_data = prog_h;
			break;

		case 1:
			viz_plane_data = prog_u;
			break;

		case 2:
			viz_plane_data = prog_v;
			break;
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 1;
	}



	const char* vis_get_status_string()
	{
		const char* description = "";
		int id = simVars.misc.vis_id % 3;

		switch (id)
		{
		default:
		case 0:
			description = "H";
			break;

		case 1:
			description = "u";
			break;

		case 2:
			description = "v";
			break;
		}


		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string,
#if SWEET_MPI
				"Rank %i - "
#endif
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
#if SWEET_MPI
				mpi_rank,
#endif
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy,
				viz_plane_data.reduce_max(),
				viz_plane_data.reduce_min()
		);

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
#endif
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
		FatalError("Timestep size not set");


	int initial_spectral_modes = simVars.disc.res_spectral[0];

	if (simVars.timecontrol.current_timestep_size < 0)
		FatalError("Timestep size not set");

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= 256; i *= 2)
	{
		simVars.timecontrol.current_timestep_size *= 0.5;

		if (simVars.disc.timestepping_method == "na_sl")
		{
			simVars.disc.res_spectral[0] = i;
			simVars.disc.res_spectral[1] = i;

//			simVars.disc.res_physical[0] = 2*i;
//			simVars.disc.res_physical[1] = i;
			simVars.disc.res_physical[0] = 0;
			simVars.disc.res_physical[1] = 0;
		}
		else
		{
			simVars.disc.res_spectral[0] = initial_spectral_modes;
			simVars.disc.res_spectral[1] = initial_spectral_modes;

			simVars.disc.res_physical[0] = 0;
			simVars.disc.res_physical[1] = 0;
		}



		planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

		std::cout << "Testing with " << planeDataConfigInstance.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << simVars.timecontrol.current_timestep_size << std::endl;

		SimulationInstance simulation;


	#if SWEET_GUI
		if (simVars.misc.gui_enabled)
		{
			planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical);
			VisSweet<SimulationInstance> visSweet(&simulation);
			return 0;
		}
		else
	#endif
		{
			while (!simulation.should_quit())
				simulation.run_timestep();

			std::cout << "Error compared to initial condition" << std::endl;
			std::cout << "Lmax error: " << simulation.max_error_h0 << std::endl;
			std::cout << "RMS error: " << simulation.rms_error_h0 << std::endl;

			if (prev_max_error >= 0)
			{
				//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
				double conv = prev_max_error / simulation.max_error_h0;
				std::cout << "Convergence: " << conv << std::endl;

				if (conv*1.1 < 4.0)
				{
					std::cerr << "Convergence not given!" << std::endl;
					exit(1);
				}
			}
			prev_max_error = simulation.max_error_h0;
			std::cout << std::endl;
		}
	}

	return 0;
}
