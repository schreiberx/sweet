/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pdeAdvectionPlane/
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#include "../include/sweet/plane/PlaneData_Spectral.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include "../programs/swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#include "test_plane_advection/Adv_Plane_TimeSteppers.hpp"



// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

SimulationVariables simVars;



class SimulationInstance
{
public:
	PlaneData_Spectral prog_h;
	PlaneData_Spectral prog_h0;	// at t0
	PlaneData_Spectral prog_u, prog_v;

	Adv_Plane_TimeSteppers timeSteppers;

	double center_x = 0.5;
	double center_y = 0.5;
	double exp_fac = 50.0;

	PlaneOperators op;

#if SWEET_GUI
	PlaneData_Physical viz_plane_data;

	int render_primitive_id = 0;
#endif

	SWEPlaneBenchmarksCombined planeBenchmarkCombined;

	double max_error_h0 = -1;
	double rms_error_h0 = -1;


public:
	SimulationInstance()	:
		prog_h(planeDataConfig),
		prog_h0(planeDataConfig),

		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		op(planeDataConfig, simVars.sim.plane_domain_size)

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


	void reset()
	{
		simVars.reset();

		planeBenchmarkCombined.setupInitialConditions(
				prog_h, prog_u, prog_v,
				simVars, op
			);

		prog_h0 = prog_h;


		// setup planeDataconfig instance again
		planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

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

		if (simVars.misc.verbosity >= 10)
			std::cout << simVars.timecontrol.current_timestep_nr << ": " << simVars.timecontrol.current_simulation_time/(60*60*24.0) << std::endl;

		max_error_h0 = (prog_h-prog_h0).toPhys().physical_reduce_max_abs();
		rms_error_h0 = (prog_h-prog_h0).toPhys().physical_reduce_rms();
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
			const PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		*o_render_primitive_id = render_primitive_id;
//		*o_bogus_data = planeDataConfig;

		int id = simVars.misc.vis_id % 3;
		switch (id)
		{
		case 0:
			viz_plane_data = prog_h.toPhys();
			break;

		case 1:
			viz_plane_data = prog_u.toPhys();
			break;

		case 2:
			viz_plane_data = prog_v.toPhys();
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
				viz_plane_data.physical_reduce_max(),
				viz_plane_data.physical_reduce_min()
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

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.sim.advection_velocity[0] == 0 && simVars.sim.advection_velocity[1] == 0)
	{
		std::cout << "At least one velocity has to be set, see parameters --advection-velocity=..." << std::endl;
		return -1;
	}

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");


	int initial_spectral_modes = simVars.disc.space_res_spectral[0];
	if (initial_spectral_modes <= 0)
	{
		SWEETError("Please specify the number of MODES");
	}

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");

	int max_modes = 256;

	// double max modes if only 1st order
	if (simVars.disc.timestepping_order == 1)
		max_modes *= 2;

	double dt = simVars.timecontrol.current_timestep_size;

	int c = 0;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		simVars.timecontrol.current_timestep_size = dt/std::pow(2.0, c);
		simVars.timecontrol.setup_timestep_size = simVars.timecontrol.current_timestep_size;

		c++;

		if (simVars.disc.timestepping_method == "na_sl")
		{
			simVars.disc.space_res_spectral[0] = i;
			simVars.disc.space_res_spectral[1] = i;

			simVars.disc.space_res_physical[0] = 0;
			simVars.disc.space_res_physical[1] = 0;
		}
		else
		{
			simVars.disc.space_res_spectral[0] = initial_spectral_modes;
			simVars.disc.space_res_spectral[1] = initial_spectral_modes;

			simVars.disc.space_res_physical[0] = 0;
			simVars.disc.space_res_physical[1] = 0;
		}



		planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

		std::cout << "Testing with " << planeDataConfigInstance.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << simVars.timecontrol.current_timestep_size << std::endl;

		SimulationInstance simulation;

	#if SWEET_GUI
		if (simVars.misc.gui_enabled)
		{
			planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);
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

				if (conv*1.1 < std::pow(2.0, (double)simVars.disc.timestepping_order))
				{
					std::cerr << "Convergence too low!" << std::endl;
					exit(1);
				}


				if (conv*0.9 > std::pow(2.0, (double)(simVars.disc.timestepping_order+1)))
				{
					std::cout<< "Convergence too high, might be still correct!" << std::endl;
					//exit(1);
				}
			}
			prev_max_error = simulation.max_error_h0;
			std::cout << std::endl;
		}
	}

	return 0;
}
