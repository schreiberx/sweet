
#include "../include/sweet/plane/PlaneData.hpp"
#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/Stopwatch.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;



SimulationVariables simVars;

int time_test_function_order = 0;
int timestepping_runge_kutta_order = 0;

class SimulationTestRK
{
public:
	PlaneData prog_h;
	PlaneData prog_u;
	PlaneData prog_v;

	PlaneDataTimesteppingExplicitRK timestepping;


	/**
	 * Function of 4th order.
	 */
public:
	double test_function(
			int i_order,
			double z
	)
	{
		switch(i_order)
		{
		case 0:
			return	2.0;

		case 1:
			return	z
					+2.0;

		case 2:
			return	+(5./12.)*z*z
					+z
					+2.0;

		case 3:
			return	-(1./2.)*z*z*z
					+(5./12.)*z*z
					+z
					+2.0;

		case 4:
			return	(1./12.)*z*z*z*z
					-(1./2.)*z*z*z
					+(5./12.)*z*z
					+z
					+2.0;
		}

		FatalError("Not supported (i_order)");
	}



public:
	double test_function_diff(int i_order, double z)
	{
		switch(i_order)
		{
		case 0:
			return	0.0;

		case 1:
			return	1.0;

		case 2:
			return	+(10./12.)*z
					+1.0;

		case 3:
			return	-(3./2.)*z*z
					+(10./12.)*z
					+1.0;

		case 4:
			return	(4./12.)*z*z*z
					-(3./2.)*z*z
					+(10./12.)*z
					+1.0;
		}
		FatalError("Not supported test_function_diff (i_order)");
	}


public:
	SimulationTestRK()	:
		prog_h(planeDataConfig),
		prog_u(planeDataConfig),
		prog_v(planeDataConfig)
	{
		reset();
	}


	void reset()
	{
		simVars.timecontrol.current_timestep_nr = 0;

		prog_h.physical_set_all(test_function(time_test_function_order, 0));
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);
	}


	void p_run_euler_timestep_update(
			const PlaneData &i_h,	///< prognostic variables
			const PlaneData &i_u,	///< prognostic variables
			const PlaneData &i_v,	///< prognostic variables

			PlaneData &o_h_t,	///< time updates
			PlaneData &o_u_t,	///< time updates
			PlaneData &o_v_t,	///< time updates

			double i_current_timestamp = -1
	)
	{
		o_h_t.physical_set_all(test_function_diff(time_test_function_order, i_current_timestamp));
		o_u_t.physical_set_all(0);
		o_v_t.physical_set_all(0);

		simVars.timecontrol.current_timestep_nr++;
	}



	void run()
	{
	}



	void run_timestep()
	{
		// either set time step size to 0 for autodetection or to
		// a positive value to use a fixed time step size
		//simVars.timecontrol.current_timestep_size = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 0);
		assert(simVars.timecontrol.current_timestep_size > 0);

		timestepping.run_timestep(
				this,
				&SimulationTestRK::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_h, prog_u, prog_v,
				simVars.timecontrol.current_timestep_size,
				timestepping_runge_kutta_order,
				simVars.timecontrol.current_simulation_time
			);

		// provide information to parameters
		simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
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
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data
	)
	{
		*o_dataArray = &prog_h;
		*o_aspect_ratio = simVars.sim.plane_domain_size[1] / simVars.sim.plane_domain_size[0];
	}

	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				simVars.diag.total_mass, simVars.diag.total_energy, simVars.diag.total_potential_enstrophy);
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

	bool instability_detected()
	{
		return !(	prog_h.reduce_boolean_all_finite() &&
					prog_u.reduce_boolean_all_finite() &&
					prog_v.reduce_boolean_all_finite()
				);
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{
	const char *bogus_var_names[] = {
			nullptr
	};

	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names, false))
		return -1;

	if (simVars.timecontrol.current_timestep_size <= 0)
		FatalError("Timestep size <= 0!");

	if (simVars.disc.space_res_physical[0] <= 0)
		simVars.disc.space_res_physical[0] = 16;
	if (simVars.disc.space_res_physical[1] <= 0)
		simVars.disc.space_res_physical[1] = 16;

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

#if SWEET_GUI
	if (simVars.misc.gui_enabled)
	{
		SimulationTestRK *simulationTestRK = new SimulationTestRK;
		VisSweet<SimulationTestRK> visSweet(simulationTestRK);
		delete simulationTestRK;
		return 0;
	}
#endif

	if (simVars.timecontrol.max_simulation_time == -1)
	{
		std::cerr << "Max. simulation time is unlimited, please specify e.g. -t 6" << std::endl;
		return -1;
	}

	for (int fun_order = 0; fun_order <= 4; fun_order++)
	{
		time_test_function_order = fun_order;

		for (int ts_order = 1; ts_order <= 4; ts_order++)
		{
			timestepping_runge_kutta_order = ts_order;

			/*
			 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
			 */
			simVars.reset();

			SimulationTestRK *simulationTestRK = new SimulationTestRK;

			while(true)
			{
				if (simVars.misc.verbosity >= 10)
					std::cout << simVars.timecontrol.current_simulation_time << ": " << simulationTestRK->prog_h.p_physical_get(0,0) << std::endl;

				simulationTestRK->run_timestep();

				if (simulationTestRK->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}

				if (simVars.timecontrol.max_simulation_time < simVars.timecontrol.current_simulation_time)
				{
					PlaneData benchmark_h(planeDataConfig);

					benchmark_h.physical_set_all(simulationTestRK->test_function(time_test_function_order, simVars.timecontrol.current_simulation_time));

					double error = (simulationTestRK->prog_h-benchmark_h).reduce_rms_quad();
					std::cout << "with function order " << fun_order << " with RK timestepping of order " << ts_order << " resulted in RMS error " << error << "\t\t" << std::flush;

					if (fun_order <= ts_order)
					{
						if (error > 1e-8)
						{
							std::cout << "ERROR threshold exceeded!" << std::endl;
							return 1;
						}

						std::cout << "OK" << std::endl;
					}
					else
					{
						std::cout << "OK (errors expected)" << std::endl;

						if (error < 1e-8)
						{
							std::cout << "WARNING: Error expected to be larger, however relatively small error found" << std::endl;
							return 1;
						}

					}

					break;
				}
			}

			delete simulationTestRK;
		}
	}


	return 0;
}
