
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include "sweet/Operators2D.hpp"
#include <unistd.h>
#include <stdio.h>

SimulationVariables simVars;


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
		h(simVars.disc.res),
		u(simVars.disc.res),
		v(simVars.disc.res),

		hu(simVars.disc.res),
		hv(simVars.disc.res),

		h_t(simVars.disc.res),

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_diffs)
	{
		reset();
	}



	void reset()
	{
		simVars.timecontrol.current_timestep_nr = 0;

		h.setAll(simVars.setup.h0);

		if (std::isinf(simVars.bogus.var[0]))
		{
			u.setAll(0);
			v.setAll(0);
		}
		else
		{
			u.setAll(simVars.bogus.var[0]);
			v.setAll(simVars.bogus.var[1]);
		}

		double center_x = 0.7;
		double center_y = 0.6;

		if (simVars.setup.scenario == 0)
		{
			/*
			 * radial dam break
			 */
			double radius = 0.2;
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.res[0];
					double y = ((double)j+0.5)/(double)simVars.disc.res[1];

					double dx = x-center_x;
					double dy = y-center_y;

					if (radius*radius >= dx*dx+dy*dy)
						h.set(j,i, simVars.setup.h0+1.0);
				}
			}
		}

		if (simVars.setup.scenario == 1)
		{
			/*
			 * fun with Gaussian
			 */
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.res[0];
					double y = ((double)j+0.5)/(double)simVars.disc.res[1];

					double dx = x-center_x;
					double dy = y-center_y;

					h.set(j,i, simVars.setup.h0+std::exp(-50.0*(dx*dx + dy*dy)));
				}
			}
		}
	}



	void run_timestep()
	{
        double dt = simVars.sim.CFL*std::min(simVars.disc.cell_size[0]/u.reduce_maxAbs(), simVars.disc.cell_size[1]/v.reduce_maxAbs());
        if (std::isinf(dt))
        	dt = simVars.sim.CFL*simVars.disc.cell_size[0]/0.000001;

        simVars.disc.timestepping_timestep_size = dt;

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
				)*(1.0/simVars.disc.cell_size[0])				// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(h)*v.return_value_if_positive()		// inflow
					-h*op.shift_down(v.return_value_if_positive())	// outflow

					// v is negative
					+(h*v.return_value_if_negative())				// outflow
					-op.shift_down(h*v.return_value_if_negative())	// inflow
				)*(1.0/simVars.disc.cell_size[1])
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
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
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
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
	}



	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Timestep: %i, timestep size: %e", simVars.timecontrol.current_timestep_nr, simVars.disc.timestepping_timestep_size);
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

	if (std::isinf(simVars.bogus.var[0]) || std::isinf(simVars.bogus.var[1]))
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

		if (simVars.misc.verbosity > 2)
			std::cout << simVars.timecontrol.current_simulation_time << std::endl;

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			break;
	}
#endif

	delete simulationSWE;

	return 0;
}
