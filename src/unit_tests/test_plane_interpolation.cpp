/*
 * test_plane_interpolation.cpp
 *
 *  Created on: 2nd April 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SWEET_GUI
#define SWEET_GUI 1
#endif

#include "../include/sweet/plane/PlaneData.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>



// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

PlaneDataConfig planeDataConfigOversamplingInstance;
PlaneDataConfig *planeDataConfigOversampling = &planeDataConfigOversamplingInstance;


SimulationVariables simVars;


class SimulationInstance
{
public:
	PlaneData prog_h;

	PlaneOperators op;

#if SWEET_GUI
	PlaneData viz_plane_data;

	int render_primitive_id = 0;
#endif

	int interpolation_order = 3;

	ScalarDataArray posx_a, posy_a;

	double center_x = 0;
	double center_y = 0;
	double exp_fac = 50.0;

	PlaneDataSampler planeDataSampler;

	double max_error;

public:
	SimulationInstance()	:
		prog_h(planeDataConfig),

		op(planeDataConfig, simVars.sim.domain_size)

#if SWEET_GUI
		,
		viz_plane_data(planeDataConfig)
#endif
	{
		reset();
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

		posx_a.setup(planeDataConfigOversampling->physical_array_data_number_of_elements);
		posy_a.setup(planeDataConfigOversampling->physical_array_data_number_of_elements);

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position

		posx_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
			{
				int i = idx % planeDataConfigOversampling->physical_res[0];
				//int j = idx / planeDataConfig->physical_data_size[0];

				io_data = simVars.sim.domain_size[0]*(double)i/(double)planeDataConfigOversampling->physical_res[0];
				assert(io_data >= 0);
				assert(io_data < simVars.sim.domain_size[0]);
			}
		);
		posy_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				//int i = idx % planeDataConfig->physical_data_size[0];
				int j = idx / planeDataConfigOversampling->physical_res[0];

				io_data = simVars.sim.domain_size[1]*(double)j/(double)planeDataConfigOversampling->physical_res[1];

				assert(io_data >= -M_PI*0.5);
				assert(io_data < simVars.sim.domain_size[1]);
			}
		);

		planeDataSampler.setup(simVars.sim.domain_size, planeDataConfig);
	}



	void run_timestep()
	{
		ScalarDataArray out_data;
		out_data.setup(posx_a.number_of_elements);

		if (interpolation_order == 2)
		{
			planeDataSampler.bilinear_scalar(
					prog_h,
					posx_a,
					posy_a,
					out_data
			);
		}
		else if (interpolation_order == 3)
		{
			planeDataSampler.bicubic_scalar(
					prog_h,
					posx_a,
					posy_a,
					out_data
			);
		}
		else
		{
			FatalError("Interpolation order not available");
		}

		max_error = 0;
		for (std::size_t i = 0; i < posx_a.number_of_elements; i++)
		{
			double value = gaussianValue(center_x, center_y, posx_a.scalar_data[i], posy_a.scalar_data[i], exp_fac);
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
//		*o_bogus_data = planeDataConfig;

		int id = simVars.misc.vis_id % 1;
		switch (id)
		{
		case 0:
			viz_plane_data = prog_h;
			break;
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 1.0;
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

	int initial_spectral_modes = simVars.disc.res_spectral[0];

	double gaussian_center_array[6][2] = {
			{0.5, 0.5},
			{0.9, 0.4},
			{0.4, 0.9},
			{0.3, 0.8},
			{0.7, 0.3},
			{0.3, 0.7},
	};

	for (int i_gaussians = 2; i_gaussians < 6; i_gaussians++)
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
				simVars.disc.res_physical[0] = 0;
				simVars.disc.res_physical[1] = 0;

				simVars.disc.res_spectral[0] = i;
				simVars.disc.res_spectral[1] = i;

				planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

				std::cout << "Testing with " << planeDataConfigInstance.getUniqueIDString() << std::endl;

				int res_physical_overs[2] = {simVars.disc.res_physical[0]*oversampling, simVars.disc.res_physical[1]*oversampling};
				int res_spectral_overs[2] = {simVars.disc.res_spectral[0]*oversampling, simVars.disc.res_spectral[1]*oversampling};

				planeDataConfigOversamplingInstance.setupAuto(res_physical_overs, res_spectral_overs);

				{
					SimulationInstance simulation;

					// Update interpolation order
					simulation.interpolation_order = interpolation_order;

					// center of Gaussian bump
					simulation.center_x = center_lon;
					simulation.center_y = center_lat;

					simulation.reset();

					if (simVars.misc.verbosity > 2)
						simVars.outputConfig();

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
