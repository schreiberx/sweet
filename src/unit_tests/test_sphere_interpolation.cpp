/*
 * test_sphere_interpolation.cpp
 *
 *  Created on: 2nd April 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SWEET_GUI
#define SWEET_GUI 1
#endif

#include "../include/sweet/sphere/SphereDataSpectral.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereDataSampler.hpp>
#include <sweet/Convert_SphereDataSpectral_To_PlaneData.hpp>
#include <sweet/Convert_SphereDataPhysical_To_PlaneData.hpp>



// Sphere data config
SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

SphereDataConfig sphereDataConfigOversamplingInstance;
SphereDataConfig *sphereDataConfigOversampling = &sphereDataConfigOversamplingInstance;


#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif

SimulationVariables simVars;


class SimulationInstance
{
public:
	SphereDataSpectral prog_h;

	SphereOperators op;

#if SWEET_GUI
	PlaneData viz_plane_data;

	int render_primitive_id = 1;
#endif

	int interpolation_order = 3;

	ScalarDataArray posx_a, posy_a;

	double center_lon = 0;
	double center_lat = 0;
	double exp_fac = 10.0;

	SphereDataSampler sphereDataSampler;

	double max_error;

public:
	SimulationInstance()	:
		prog_h(sphereDataConfig),

		op(sphereDataConfig, simVars.sim.sphere_radius)

#if SWEET_GUI
		,
		viz_plane_data(planeDataConfig)
#endif
	{
		reset();
	}



	static
	double gaussianValue(
			double i_center_lon, double i_center_lat,
			double i_lon, double i_lat,
			double i_exp_fac
	)
	{
		double center_lon = i_center_lon;
		double center_lat = i_center_lat;

		double mu = std::sin(i_lat);
		double phi1 = asin(mu);
		double phi2 = center_lat;
		double lambda1 = i_lon;
		double lambda2 = center_lon;

		double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		return std::exp(-d*d*i_exp_fac)*0.1*simVars.sim.h0 + simVars.sim.h0;
	}



	void reset()
	{
		simVars.reset();

		SphereDataSpectral tmp_vort(sphereDataConfig);
		SphereDataSpectral tmp_div(sphereDataConfig);

		SphereDataPhysical prog_h_phys(sphereDataConfig);
		prog_h_phys.physical_update_lambda(
			[&](double i_lon, double i_lat, double &o_data)
			{
				o_data = gaussianValue(center_lon, center_lat, i_lon, i_lat, exp_fac);
			}
		);
		prog_h.loadSphereDataPhysical(prog_h_phys);

		posx_a.setup(sphereDataConfigOversampling->physical_array_data_number_of_elements);
		posy_a.setup(sphereDataConfigOversampling->physical_array_data_number_of_elements);

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position

		posx_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
			{
				int i = idx % sphereDataConfigOversampling->physical_num_lon;

				io_data = 2.0*M_PI*(double)i/(double)sphereDataConfigOversampling->physical_num_lon;
				assert(io_data >= 0);
				assert(io_data < 2.0*M_PI);
			}
		);
		posy_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				//int i = idx % sphereDataConfig->physical_data_size[0];
				int j = idx / sphereDataConfigOversampling->physical_num_lon;

				io_data = sphereDataConfigOversampling->lat[j];

				assert(io_data >= -M_PI*0.5);
				assert(io_data <= M_PI*0.5);
			}
		);

		sphereDataSampler.setup(sphereDataConfig);
	}



	void run_timestep()
	{
		ScalarDataArray out_data;
		out_data.setup(posx_a.number_of_elements);

		if (interpolation_order == 2)
		{
			sphereDataSampler.bilinear_scalar(
					prog_h,
					posx_a,
					posy_a,
					out_data,
					false
			);
		}
		else if (interpolation_order == 3)
		{
			sphereDataSampler.bicubic_scalar(
					prog_h,
					posx_a,
					posy_a,
					out_data,
					false
			);
		}
		else
		{
			FatalError("Interpolation order not available");
		}

		max_error = 0;
		for (std::size_t i = 0; i < posx_a.number_of_elements; i++)
		{
			double value = gaussianValue(center_lon, center_lat, posx_a.scalar_data[i], posy_a.scalar_data[i], exp_fac);
			max_error = std::max(max_error, std::abs(value - out_data.scalar_data[i]));
		}
	}


	bool should_quit()
	{
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
			{
				run_timestep();
				std::cout << max_error << std::endl;
			}
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

		int id = simVars.misc.vis_id % 1;
		switch (id)
		{
		case 0:
			viz_plane_data = Convert_SphereDataSpectral_To_PlaneData::physical_convert(prog_h, planeDataConfig);
			break;
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	const char* vis_get_status_string()
	{
		const char* description = "";
		int id = simVars.misc.vis_id % 1;

		switch (id)
		{
		default:
		case 0:
			description = "H";
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

	int initial_spectral_modes = simVars.disc.space_res_spectral[0];

	double gaussian_center_array[6][2] = {
			{0.0, M_PI*0.5},	// equator
			{0.0, 0.0},			// equator
			{0.0, -M_PI*0.5},	// bottom
			{0.0, M_PI*0.5},	// top
			{1.0, M_PI*0.4},	// misc
			{2.3, -M_PI*0.34},	// misc
	};

	for (int i_gaussians = 0; i_gaussians < 6; i_gaussians++)
	{
		double center_lon = gaussian_center_array[i_gaussians][0];
		double center_lat = gaussian_center_array[i_gaussians][1];

		std::cout << "*********************************************************" << std::endl;
		std::cout << "* Running studies for Gaussian at " << center_lon << ", " << center_lat << std::endl;
		std::cout << "*********************************************************" << std::endl;

		for (int interpolation_order = 2; interpolation_order <= 3; interpolation_order++)
		{
//interpolation_order = 3;
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* Running studies for interpolation of order " << interpolation_order << std::endl;
			std::cout << "*********************************************************" << std::endl;

			int oversampling = 13;
			std::cout << "Using oversampling of " << oversampling << std::endl;

			double prev_max_error = -1;
			for (int i = initial_spectral_modes; i <= 256; i *= 2)
			{
				simVars.disc.space_res_physical[0] = 2*i;
				simVars.disc.space_res_physical[1] = i;

				simVars.disc.space_res_spectral[0] = i;
				simVars.disc.space_res_spectral[1] = i;

				sphereDataConfigInstance.setupAuto(
						simVars.disc.space_res_physical,
						simVars.disc.space_res_spectral,
						simVars.misc.reuse_spectral_transformation_plans
					);

				std::cout << "Testing with " << sphereDataConfigInstance.getUniqueIDString() << std::endl;

				int res_physical_overs[2] = {simVars.disc.space_res_physical[0]*oversampling, simVars.disc.space_res_physical[1]*oversampling};
				int res_spectral_overs[2] = {simVars.disc.space_res_spectral[0]*oversampling, simVars.disc.space_res_spectral[1]*oversampling};

				sphereDataConfigOversamplingInstance.setupAuto(
						res_physical_overs,
						res_spectral_overs,
						simVars.misc.reuse_spectral_transformation_plans
					);

				{
					SimulationInstance simulation;

					// Update interpolation order
					simulation.interpolation_order = interpolation_order;

					// center of Gaussian bump
					simulation.center_lon = center_lon;
					simulation.center_lat = center_lat;

					simulation.reset();

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
						simulation.reset();
						simulation.run_timestep();
						std::cout << "Lmax error: " << simulation.max_error << std::endl;

						if (prev_max_error >= 0)
						{
							//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
							double conv = prev_max_error / simulation.max_error;
							std::cout << "Convergence: " << conv << std::endl;

							if (conv*1.1 < std::pow(2.0, interpolation_order))
								FatalError("Convergence not given!");
						}
						prev_max_error = simulation.max_error;
					}
				}
			}
		}
	}

	return 0;
}
