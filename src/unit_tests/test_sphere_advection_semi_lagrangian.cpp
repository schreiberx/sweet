/*
 * test_sphere_advection.cpp
 *
 *  Created on: 3 Apr 2018
 *      Author: martin
 */


#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#include <sweet/sphere/SphereData_Spectral.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <benchmarks_sphere/SWESphereBenchmarksCombined.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/Convert_SphereDataSpectral_To_PlaneData.hpp>
#include <sweet/Convert_SphereDataPhysical_To_PlaneData.hpp>

#include "test_sphere_advection_semi_lagrangian/Adv_Sphere_TimeSteppers.hpp"



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
	SphereData_Spectral prog_phi_pert;
	SphereData_Spectral prog_phi_pert_t0;	// at t0
	SphereData_Spectral prog_vort, prog_div;

	Adv_Sphere_TimeSteppers timeSteppers;

	SphereOperators_SphereData op;

	bool time_varying_fields;

	/*
	 * LMax error to h0
	 */
	double max_error_h0;
	/*
	 * RMS error to h0
	 */
	double rms_error_h0;


#if SWEET_GUI
	PlaneData viz_plane_data;

	int render_primitive_id = 1;
#endif

	SWESphereBenchmarksCombined sphereBenchmarks;



public:
	SimulationInstance()	:
		prog_phi_pert(sphereDataConfig),
		prog_phi_pert_t0(sphereDataConfig),

		prog_vort(sphereDataConfig),
		prog_div(sphereDataConfig),

		op(sphereDataConfig, &(simVars.sim)),
		time_varying_fields(false)

#if SWEET_GUI
		,
		viz_plane_data(planeDataConfig)
#endif
	{
		reset();
	}



	void reset()
	{
		simVars.reset();

		sphereBenchmarks.setup(simVars, op);
		SphereData_Physical o_phi_pert_phys(sphereDataConfig);
		sphereBenchmarks.compute_initial_condition_pert(
				prog_phi_pert, prog_vort, prog_div,
				&o_phi_pert_phys,
				&time_varying_fields
			);

		prog_phi_pert_t0 = prog_phi_pert;

		// setup sphereDataconfig instance again
		sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

		timeSteppers.setup(simVars.disc.timestepping_method, op, simVars);
	}



	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		/*
		 * Update time varying fields
		 */
		if (time_varying_fields)
			sphereBenchmarks.update_time_varying_fields_pert(prog_phi_pert, prog_vort, prog_div, simVars.timecontrol.current_simulation_time);

		SphereData_Physical asdf(sphereDataConfig);
		timeSteppers.master->run_timestep(
				prog_phi_pert, prog_vort, prog_div,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time,
				&sphereBenchmarks,
				asdf
			);

		double dt = simVars.timecontrol.current_timestep_size;

		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.misc.verbosity >= 10)
			std::cout << simVars.timecontrol.current_timestep_nr << ": " << simVars.timecontrol.current_simulation_time/(60*60*24.0) << std::endl;

		max_error_h0 = (prog_phi_pert_t0-prog_phi_pert).getSphereDataPhysical().physical_reduce_max_abs();
		rms_error_h0 = (prog_phi_pert_t0-prog_phi_pert).getSphereDataPhysical().physical_reduce_rms();
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
	}


	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max
	)
	{
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		int id = simVars.misc.vis_id % 5;
		switch (id)
		{
		case 0:
			viz_plane_data = Convert_SphereDataSpectral_To_PlaneData::physical_convert(prog_phi_pert, planeDataConfig);
			break;

		case 1:
			viz_plane_data = Convert_SphereDataSpectral_To_PlaneData::physical_convert(prog_vort, planeDataConfig);
			break;

		case 2:
			viz_plane_data = Convert_SphereDataSpectral_To_PlaneData::physical_convert(prog_div, planeDataConfig);
			break;

		case 3:
			{
				SphereData_Physical u(sphereDataConfig);
				SphereData_Physical v(sphereDataConfig);

//				op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);
				op.vortdiv_to_uv(prog_vort, prog_div, u, v);
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(u, planeDataConfig);
			}
			break;

		case 4:
			{
				SphereData_Physical u(sphereDataConfig);
				SphereData_Physical v(sphereDataConfig);

//				op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);
				op.vortdiv_to_uv(prog_vort, prog_div, u, v);
				viz_plane_data = Convert_SphereDataPhysical_To_PlaneData::physical_convert(v, planeDataConfig);
			}
			break;

		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	const char* vis_get_status_string()
	{
		const char* description = "";
		int id = simVars.misc.vis_id % 5;

		switch (id)
		{
		default:
		case 0:
			description = "H";
			break;

		case 1:
			description = "vort";
			break;

		case 2:
			description = "div";
			break;

		case 3:
			description = "u";
			break;

		case 4:
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
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
	{
		std::cout << std::endl;
		return -1;
	}

	simVars.outputConfig();

	int initial_spectral_modes = simVars.disc.space_res_spectral[0];

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");

	int max_modes = 256;

	if (simVars.disc.timestepping_order == 1)
		max_modes = 512;
	else if (simVars.disc.timestepping_order == 2)
		max_modes = 256;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		simVars.timecontrol.current_timestep_size *= 0.5;
		simVars.timecontrol.setup_timestep_size = simVars.timecontrol.current_timestep_size;

		if (simVars.disc.timestepping_method == "na_sl")
		{
			simVars.disc.space_res_spectral[0] = i;
			simVars.disc.space_res_spectral[1] = i;

			simVars.disc.space_res_physical[0] = 2*i;
			simVars.disc.space_res_physical[1] = i;
		}
		else
		{
			simVars.disc.space_res_spectral[0] = initial_spectral_modes;
			simVars.disc.space_res_spectral[1] = initial_spectral_modes;

			simVars.disc.space_res_physical[0] = 0;
			simVars.disc.space_res_physical[1] = 0;
		}



		sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

		if (1)
		{
			std::cout << "Checking for right velocity" << std::endl;
			SphereTimestepping_SemiLagrangian sl(simVars, sphereDataConfig);


			/*
			 * Convert to Cartesian velocity space
			 */
			ScalarDataArray V_lon_D(sphereDataConfig->physical_array_data_number_of_elements);
			ScalarDataArray V_lat_D(sphereDataConfig->physical_array_data_number_of_elements);
			V_lon_D.set_all(1.0);
			V_lat_D.set_all(2.0);

			ScalarDataArray V_x_D, V_y_D, V_z_D;
			SWEETMath::latlon_velocity_to_cartesian_velocity(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_lon_D,
					V_lat_D,
					V_x_D,
					V_y_D,
					V_z_D
				);

			ScalarDataArray V_lon_tmp, V_lat_tmp;
			SWEETMath::cartesian_velocity_to_latlon_velocity(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_x_D,
					V_y_D,
					V_z_D,
					V_lon_tmp, V_lat_tmp
			);

			double err_lon = (V_lon_D - V_lon_tmp).reduce_maxAbs();
			double err_lat = (V_lat_D - V_lat_tmp).reduce_maxAbs();

			if (err_lon > 1e-10)
			{
				std::cerr << "Error: " << err_lon << std::endl;
				SWEETError("Error lon too high!");
			}

			if (err_lat > 1e-10)
			{
				std::cerr << "Error: " << err_lat << std::endl;
				SWEETError("Error lat too high!");
			}
		}

		std::cout << "Testing with " << sphereDataConfigInstance.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << simVars.timecontrol.current_timestep_size << std::endl;

		SimulationInstance simulation;


#if SWEET_GUI
		if (simVars.misc.gui_enabled)
		{
			planeDataConfigInstance.setupAutoSpectralSpace(
					simVars.disc.space_res_physical,
					simVars.misc.reuse_spectral_transformation_plans
				);
			VisSweet<SimulationInstance> visSweet(&simulation);
			return 0;
		}
		else
#endif
		{
//			simulation.reset();
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
					SWEETError("Convergence not given!");
			}

			if (simulation.max_error_h0  > 1e10)
				SWEETError("Lmax error exceeded threshold!");

			prev_max_error = simulation.max_error_h0;

			std::cout << "*********************************************" << std::endl;
		}
	}

	return 0;
}
