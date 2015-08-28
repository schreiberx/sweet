
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include "sweet/SimulationParameters.hpp"
#include "sweet/Operators2D.hpp"
#include <unistd.h>
#include <stdio.h>

SimulationParameters parameters;


class SimulationSWE
{
public:
	DataArray<2> h;
	DataArray<2> u;
	DataArray<2> v;

	DataArray<2> hu;
	DataArray<2> hv;

	DataArray<2> h_t;

	Operators2D op;

public:
	SimulationSWE()	:
		h(parameters.res),
		u(parameters.res),
		v(parameters.res),

		hu(parameters.res),
		hv(parameters.res),

		h_t(parameters.res),

		op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs)
	{
		reset();
	}



	void reset()
	{
		parameters.status_timestep_nr = 0;

		h.setAll(parameters.setup_h0);

		if (std::isinf(parameters.bogus_var0))
		{
			u.setAll(0);
			v.setAll(0);
		}
		else
		{
			u.setAll(parameters.bogus_var0);
			v.setAll(parameters.bogus_var1);
		}

		double center_x = 0.7;
		double center_y = 0.6;

		if (parameters.setup_scenario == 0)
		{
			/*
			 * radial dam break
			 */
			double radius = 0.2;
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

					double dx = x-center_x;
					double dy = y-center_y;

					if (radius*radius >= dx*dx+dy*dy)
						h.set(j,i, parameters.setup_h0+1.0);
				}
			}
		}

		if (parameters.setup_scenario == 1)
		{
			/*
			 * fun with Gaussian
			 */
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

					double dx = x-center_x;
					double dy = y-center_y;

					h.set(j,i, parameters.setup_h0+std::exp(-50.0*(dx*dx + dy*dy)));
				}
			}
		}
	}



	void run_timestep()
	{
        double dt = parameters.sim_CFL*std::min(parameters.sim_cell_size[0]/u.reduce_maxAbs(), parameters.sim_cell_size[1]/v.reduce_maxAbs());
        if (std::isinf(dt))
        	dt = parameters.sim_CFL*parameters.sim_cell_size[0]/0.000001;

        parameters.timestepping_timestep_size = dt;

        // 0: staggered
        // 1: non-staggered
        // 2: up/downwinding

#define GRID_LAYOUT_AND_ADVECTION	2

#if ADVECTION_METHOD == 0

        // staggered
		h -= dt*(
				op.diff_b_x(op.avg_f_x(h)*u) +
				op.diff_b_y(op.avg_f_y(h)*v)
			);
#endif

#if ADVECTION_METHOD == 1
		// non-staggered
		h = h - dt*(
				op.diff_c_x(h*u) +
				op.diff_c_y(h*v)
			);
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
				)*(1.0/parameters.sim_cell_size[0])				// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(h)*v.return_value_if_positive()		// inflow
					-h*op.shift_down(v.return_value_if_positive())	// outflow

					// v is negative
					+(h*v.return_value_if_negative())				// outflow
					-op.shift_down(h*v.return_value_if_negative())	// inflow
				)*(1.0/parameters.sim_cell_size[1])
			);
#endif
		parameters.status_timestep_nr++;
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
		if (parameters.run_simulation)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}


	void vis_get_vis_data_array(
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		switch (parameters.vis_id)
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
		*o_aspect_ratio = parameters.sim_domain_size[1] / parameters.sim_domain_size[0];
	}



	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Timestep: %i, timestep size: %e, Mass: %e, Energy: %e, Potential Entrophy: %e", parameters.status_timestep_nr, parameters.timestepping_timestep_size, parameters.diagnostics_mass, parameters.diagnostics_energy, parameters.diagnostics_potential_entrophy);
		return title_string;
	}

	void vis_pause()
	{
		parameters.run_simulation = !parameters.run_simulation;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			parameters.vis_id++;
			break;

		case 'V':
			parameters.vis_id--;
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

	if (!parameters.setup(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Program-specific options:" << std::endl;
		std::cout << "	--velocity-u [advection velocity u]" << std::endl;
		std::cout << "	--velocity-v [advection velocity v]" << std::endl;
		return -1;
	}

	if (std::isinf(parameters.bogus_var0) || std::isinf(parameters.bogus_var1))
	{
		std::cout << "Both velocities have to be set, see parameters --velocity-u, --velocity-v" << std::endl;
		return -1;
	}


	SimulationSWE *simulationSWE = new SimulationSWE;

#if SWEET_GUI
	VisSweet<SimulationSWE> visSweet(simulationSWE);
#else
	simulationSWE->reset();
	while (!simulationSWE->should_quit())
	{
		simulationSWE->run_timestep();

		if (parameters.verbosity > 2)
			std::cout << parameters.status_simulation_time << std::endl;

		if (parameters.status_simulation_time > parameters.max_simulation_time)
			break;
	}
#endif

	delete simulationSWE;

	return 0;
}
