/*
 * lagrangian_test.cpp
 *
 *  Created on: 3 Dec 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */


#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Sampler2D.hpp>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <array>

SimulationVariables simVars;


class SimulationSWE
{
public:
	DataArray<2> prog_h;
	DataArray<2> prog_u;
	DataArray<2> prog_v;

	DataArray<2> interpol_values;
//	DataArray<2> hu;
//	DataArray<2> hv;

	DataArray<2> h_t;

	Operators2D op;
	Sampler2D sampler2D;


public:
	SimulationSWE()	:
		prog_h(simVars.disc.res),
		prog_u(simVars.disc.res),
		prog_v(simVars.disc.res),

		interpol_values(simVars.disc.res),

//		hu(simVars.disc.res),
//		hv(simVars.disc.res),

		h_t(simVars.disc.res),

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_diffs)
	{
		reset();
	}



	void reset()
	{
		simVars.timecontrol.current_timestep_nr = 0;
#if 0
		prog_h.set_all(simVars.setup.h0);

		if (std::isinf(simVars.bogus.var[0]))
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


		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
				double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

				prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
				prog_u.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
				prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
			}
		}

		prog_u += 4.0;
		prog_v += 4.0;

		sampler2D.setup(simVars.sim.domain_size, simVars.disc.res);

	}



	void run_timestep()
	{
		double dt = simVars.sim.domain_size[0] / simVars.disc.res[0];
		simVars.timecontrol.current_timestep_size = dt;


		// setup x/y coordinates
		DataArray<2> data_x(prog_h.resolution);
		DataArray<2> data_y(prog_h.resolution);

		DataArray<2> *x[2] = {&data_x, &data_y};

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position
		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
				x[0]->set(j, i, ((double)i)*(simVars.sim.domain_size[0]/(double)simVars.disc.res[0])+simVars.timecontrol.current_timestep_nr);
				x[1]->set(j, i, ((double)j)*(simVars.sim.domain_size[1]/(double)simVars.disc.res[1])+simVars.timecontrol.current_timestep_nr*2);
			}
		}

		sampler2D.bicubic_scalar(
				prog_h,	///< input scalar field
				x,		///< sampling positions
				prog_u	///< output field
		);

		// interpolate data from x/y coordinates
//		interpol_values = op_sample.bilinear(prog_h, x, y);

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
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
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
#if 1
	const char *bogus_var_names[] = {
			"velocity-u",
			"velocity-v",
			nullptr
	};

	if (!simVars.setupFromMainParameters(i_argc, i_argv/*, bogus_var_names*/))
	{
		std::cout << std::endl;
		std::cout << "Program-specific options:" << std::endl;
		std::cout << "	--velocity-u [advection velocity u]" << std::endl;
		std::cout << "	--velocity-v [advection velocity v]" << std::endl;
		return -1;
	}
/*
	if (std::isinf(simVars.bogus.var[0]) || std::isinf(simVars.bogus.var[1]))
	{
		std::cout << "Both velocities have to be set, see parameters --velocity-u, --velocity-v" << std::endl;
		return -1;
	}
*/
#endif

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
