/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif

#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <unistd.h>
#include <stdio.h>

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;



SimulationVariables simVars;

double velx = 0;
double vely = 0;

class SimulationSWE
{
public:
	PlaneData_Physical h;
	PlaneData_Physical u;
	PlaneData_Physical v;

	PlaneData_Physical hu;
	PlaneData_Physical hv;

	PlaneData_Physical h_t;

	PlaneOperators op;

public:
	SimulationSWE()	:
		h(planeDataConfig),
		u(planeDataConfig),
		v(planeDataConfig),

		hu(planeDataConfig),
		hv(planeDataConfig),

		h_t(planeDataConfig),

		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)
	{
		reset();
	}



	void reset()
	{
		//double cell_size_x = simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
		//double cell_size_y = simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];

		simVars.timecontrol.current_timestep_nr = 0;

		h.physical_set_all_value(simVars.sim.h0);
		u.physical_set_all_value(velx);
		v.physical_set_all_value(vely);

		double center_x = 0.7;
		double center_y = 0.6;

		if (simVars.benchmark.benchmark_name == "advection_benchmark_0")
		{
			/*
			 * radial dam break
			 */
			double radius = 0.2;

			h.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.space_res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.space_res_physical[1];

					double dx = x-center_x;
					double dy = y-center_y;

					if (radius*radius >= dx*dx+dy*dy)
						io_data = simVars.sim.h0+1.0;
				}
			);
		}

		if (simVars.benchmark.benchmark_name == "advection_benchmark_1")
		{
			/*
			 * fun with Gaussian
			 */

			h.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.space_res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.space_res_physical[1];

					double dx = x-center_x;
					double dy = y-center_y;

					io_data = simVars.sim.h0+std::exp(-50.0*(dx*dx + dy*dy));
				}
			);
		}
	}



	void run_timestep()
	{
		//double cell_size_x = simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
		//double cell_size_y = simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];

		double dt = simVars.timecontrol.current_timestep_size;//*std::min(cell_size_x/u.reduce_maxAbs(), cell_size_y/v.reduce_maxAbs());

		if (std::isinf(dt) != 0)
			dt = simVars.timecontrol.current_timestep_size;//*cell_size_x/0.000001;

		simVars.timecontrol.current_timestep_size = dt;

		// 0: staggered
		// 1: non-staggered
		// 2: up/downwinding

#define GRID_LAYOUT_AND_ADVECTION	2

#if ADVECTION_METHOD == 0

        // staggered
		h -= dt*(
				op.diff_b_x(op.avg_f_x(h)*u) +
				op.diff_b_y(op.avg_f_y(h)*v)
			).toPhys();
#endif

#if ADVECTION_METHOD == 1
		// non-staggered
		h = h - dt*(
				op.diff_c_x(h*u) +
				op.diff_c_y(h*v)
			).toPhys();
#endif

#if ADVECTION_METHOD == 2

		h += dt*
			(
				(
					// u is positive
					op.shift_right(h)*u.return_value_if_positive()	// inflow
					-h*op.shift_left(u.return_value_if_positive())	// outflow

					// u is negative
					+(h*u.return_value_if_negative())				// outflow
					-op.shift_left(h*u.return_value_if_negative())	// inflow
				)*(1.0/simVars.disc.cell_size_x)				// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(h)*v.return_value_if_positive()		// inflow
					-h*op.shift_down(v.return_value_if_positive())	// outflow

					// v is negative
					+(h*v.return_value_if_negative())				// outflow
					-op.shift_down(h*v.return_value_if_negative())	// inflow
				)*(1.0/simVars.disc.cell_size_y)
			);
#endif
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
			const PlaneData_Physical **o_dataArray,
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
			*o_dataArray = &h;
			break;

		case 1:
			*o_dataArray = &u;
			break;

		case 2:
			*o_dataArray = &v;
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
	const char *user_defined_parameters[] = {
			"velocity-u",
			"velocity-v",
			nullptr
	};

	if (!simVars.setupFromMainParameters(i_argc, i_argv, user_defined_parameters))
	{
		std::cout << std::endl;
		std::cout << "Program-specific options:" << std::endl;
		std::cout << "	--velocity-u [advection velocity u]" << std::endl;
		std::cout << "	--velocity-v [advection velocity v]" << std::endl;
		return -1;
	}

	if (std::isinf(atof(simVars.user_defined.var[0].c_str()) != 0))
		velx = atof(simVars.user_defined.var[0].c_str());

	if (std::isinf(atof(simVars.user_defined.var[1].c_str()) != 0))
		vely = atof(simVars.user_defined.var[1].c_str());

	if (velx == 0 && vely == 0)
	{
		std::cout << "At least one velocity has to be set, see parameters --velocity-u, --velocity-v" << std::endl;
		return -1;
	}

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

	SimulationSWE *simulationSWE = new SimulationSWE;

#if SWEET_GUI
	VisSweet<SimulationSWE> visSweet(simulationSWE);
#else
	simulationSWE->reset();
	while (!simulationSWE->should_quit())
	{
		simulationSWE->run_timestep();

		if (simVars.misc.verbosity > 2)
			std::cout << simVars.timecontrol.current_simulation_time << std::endl;

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			break;
	}
#endif

	delete simulationSWE;

	return 0;
}
