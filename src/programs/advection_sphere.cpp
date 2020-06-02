/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/advection_sphere_timeintegrators
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/advection_sphere_benchmarks
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif

#include <sweet/sphere/SphereData_Spectral.hpp>
#include "advection_sphere_benchmarks/BenchmarksSphereAdvection.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/Convert_SphereDataSpectral_To_PlaneData.hpp>
#include <sweet/Convert_SphereDataPhysical_To_PlaneData.hpp>

#include <sweet/sphere/SphereData_DebugContainer.hpp>

#include <programs/advection_sphere_timeintegrators/SphereAdvectionTimeSteppers.hpp>



// Sphere data config
SphereData_Config sphereDataConfigInstance;
SphereData_Config *sphereDataConfig = &sphereDataConfigInstance;

#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif

SimulationVariables simVars;


class SimulationInstance
{
public:
	std::vector<SphereData_Spectral*> prognostic_variables;
	std::vector<SphereData_Spectral*> prognostic_variables_t0;

	//SphereData_Spectral prog_vort, prog_div;
	SphereData_Physical velocity_field_u, velocity_field_v;

	SphereAdvectionTimeSteppers timeSteppers;


	SphereOperators_SphereData op;

	bool time_varying_fields;

#if SWEET_GUI
	PlaneData viz_plane_data;

	int render_primitive_id = 1;
#endif

	BenchmarksSphereAdvection sphereBenchmarks;


public:
	SimulationInstance()	:
		velocity_field_u(sphereDataConfig),
		velocity_field_v(sphereDataConfig),

		op(sphereDataConfig, &(simVars.sim)),
		time_varying_fields(false)

#if SWEET_GUI
		,
		viz_plane_data(planeDataConfig)
#endif
	{
		reset();
	}


	~SimulationInstance()
	{
		free_prognostic_variables();
	}


	void alloc_prognostic_variables(std::size_t i_size)
	{
		prognostic_variables.resize(i_size);
		for (std::size_t i = 0; i < prognostic_variables.size(); i++)
			prognostic_variables[i] = new SphereData_Spectral(sphereDataConfig);

		prognostic_variables_t0.resize(i_size);
		for (std::size_t i = 0; i < prognostic_variables_t0.size(); i++)
			prognostic_variables_t0[i] = new SphereData_Spectral(sphereDataConfig);
	}


	void free_prognostic_variables()
	{
		for (std::size_t i = 0; i < prognostic_variables.size(); i++)
			delete prognostic_variables[i];
		prognostic_variables.clear();

		for (std::size_t i = 0; i < prognostic_variables_t0.size(); i++)
			delete prognostic_variables_t0[i];
		prognostic_variables_t0.clear();
	}


	void reset()
	{
		simVars.reset();

		SphereData_Spectral tmp_vort(sphereDataConfig);
		SphereData_Spectral tmp_div(sphereDataConfig);

		sphereBenchmarks.setup(simVars, op);

		int num_field_variables = sphereBenchmarks.master->get_num_prognostic_fields();

		free_prognostic_variables();
		alloc_prognostic_variables(num_field_variables);

		sphereBenchmarks.master->get_initial_state(
				prognostic_variables, velocity_field_u, velocity_field_v
			);

		for (std::size_t i = 0; i < prognostic_variables.size(); i++)
			*prognostic_variables_t0[i] = *prognostic_variables[i];

		// has this benchmark time-varying fields?
		time_varying_fields = sphereBenchmarks.master->has_time_varying_state();

		// setup sphereDataconfig instance again
		sphereDataConfigInstance.setupAuto(
				simVars.disc.space_res_physical,
				simVars.disc.space_res_spectral,
				simVars.misc.reuse_spectral_transformation_plans
			);

		timeSteppers.setup(simVars.disc.timestepping_method, op, simVars);

		simVars.outputConfig();
	}



	void run_timestep()
	{
		SphereData_DebugContainer::clear();

		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;


		timeSteppers.master->run_timestep(
				prognostic_variables, velocity_field_u, velocity_field_v,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time,
				(time_varying_fields ? &sphereBenchmarks : nullptr)
			);

		if (time_varying_fields)
		{
			/*
			 * Update velocities just for sake of the correction visualization
			 */
			sphereBenchmarks.master->get_varying_velocities(
					velocity_field_u,
					velocity_field_v,
					simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size
				);
		}

		double dt = simVars.timecontrol.current_timestep_size;

		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.misc.verbosity > 2)
			std::cout << simVars.timecontrol.current_timestep_nr << ": " << simVars.timecontrol.current_simulation_time/(60*60*24.0) << std::endl;
	}



	bool should_quit()
	{
		if (simVars.misc.gui_enabled)
			return false;

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
	}


	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		if (simVars.misc.vis_id < 0)
		{
			int id = -simVars.misc.vis_id-1;

#if 0
			if (id <  (int)SphereData_DebugContainer().size())
			{
				SphereData_DebugContainer::DataContainer &d = SphereData_DebugContainer().container_data()[id];
				if (d.is_spectral)
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneData::physical_convert(d.data_spectral, planeDataConfig);
				else
					viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(d.data_physical, planeDataConfig);

				*o_dataArray = &viz_plane_data;
				*o_aspect_ratio = 0.5;
				return;
			}
#else
			if (id <  (int)prognostic_variables.size())
			{
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneData::physical_convert(
							*prognostic_variables[id] - *prognostic_variables_t0[id],
							planeDataConfig
						);

				*o_dataArray = &viz_plane_data;
				*o_aspect_ratio = 0.5;
				return;
			}
#endif
		}

		std::size_t id = simVars.misc.vis_id;

		if (id >= 0 && id < prognostic_variables.size())
		{
			viz_plane_data = Convert_SphereDataSpectral_To_PlaneData::physical_convert(*prognostic_variables[id], planeDataConfig);
		}
		else if (id >= prognostic_variables.size() && id < prognostic_variables.size() + 2)
		{
			switch (id - prognostic_variables.size())
			{
			case 0:
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(velocity_field_u, planeDataConfig);
				break;

			case 1:
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(velocity_field_v, planeDataConfig);
				break;
			}
		}
		else if (id >= prognostic_variables.size() + 2 && id < prognostic_variables.size() + 4)
		{
			SphereData_Physical u, v;
			op.vrtdiv_to_uv(*prognostic_variables[0], *prognostic_variables[1], u, v);

			switch (id - prognostic_variables.size() - 2)
			{
			case 0:
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(u, planeDataConfig);
				break;

			case 1:
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(v, planeDataConfig);
				break;
			}
		}
		else
		{
			viz_plane_data.physical_set_zero();
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	const char* vis_get_status_string()
	{
		std::string description = "";

		bool found = false;
		if (simVars.misc.vis_id < 0)
		{
			int id = -simVars.misc.vis_id-1;

#if 0
			if (id < (int)SphereData_DebugContainer().size())
			{
				description = std::string("DEBUG_")+SphereData_DebugContainer().container_data()[id].description;
				found = true;
			}
#else
			if (id <  (int)prognostic_variables.size())
			{
				std::ostringstream msg;
				msg << "DIFF prog. field " << id;
				description = msg.str();
			}
#endif
		}

		if (!found)
		{
			std::size_t id = simVars.misc.vis_id;

			if (id >= 0 && id < prognostic_variables.size())
			{
				std::ostringstream msg;
				msg << "Prog. field " << id;
				description = msg.str();
			}
			else if (id >= prognostic_variables.size() && id < prognostic_variables.size() + 2)
			{
				switch (id - prognostic_variables.size())
				{
				case 0:
					description = "u velocity";
					break;

				case 1:
					description = "v velocity";
					break;
				}
			}
			else if (id >= prognostic_variables.size() + 2 && id < prognostic_variables.size() + 4)
			{
				switch (id - prognostic_variables.size() - 2)
				{
				case 0:
					description = "prognostic field: u velocity";
					break;

				case 1:
					description = "prognostic field: v velocity";
					break;
				}
			}
			else
			{
				description = "field doesn't exist";
			}
		}


		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string,
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e "
				","
				"Colorscale: lowest [Blue... green ... red] highest",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description.c_str(),
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
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
	{
		std::cout << std::endl;
		return -1;
	}

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");

	sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

	SimulationInstance *simulation = new SimulationInstance;

#if SWEET_GUI
	if (simVars.misc.gui_enabled)
	{
		planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);

		VisSweet<SimulationInstance> visSweet(simulation);
	}
	else
#endif
	{
		simulation->reset();
		while (!simulation->should_quit())
		{
			simulation->run_timestep();

			if (simVars.timecontrol.max_simulation_time != -1)
				if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
					break;

			if (simVars.timecontrol.max_timesteps_nr != -1)
				if (simVars.timecontrol.current_timestep_nr > simVars.timecontrol.max_timesteps_nr)
					break;
		}
	}

	delete simulation;

	return 0;
}
