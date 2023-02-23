/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pdeAdvectionPlane/
 * MULE_COMPILE_FILES_AND_DIRS: src/include/sweet/plane/
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <sweet/defaultPrecompilerValues.hpp>

#include <sweet/ErrorBase.hpp>
#include <sweet/ProgramArguments.hpp>

#include <sweet/plane/Plane.hpp>

#include <sweet/shacks/ShackDictionary.hpp>
//#include <sweet/shacksShared/ShackDiagnostics.hpp>

#include <sweet/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/shacksShared/ShackMisc.hpp>

#include "pdeAdvectionPlane/ShackPDEAdvectionPlane.hpp"
#include "pdeAdvectionPlaneBenchmarks/ShackPDEAdvectionPlaneBenchmarks.hpp"

// Benchmarks
#include "pdeAdvectionPlaneBenchmarks/PDEAdvectionPlaneBenchmarksCombined.hpp"

// Time steppers
#include "pdeAdvectionPlane/PDEAdvPlaneTimeSteppers.hpp"


#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif



class SimulationInstance
{
public:
	sweet::ErrorBase error;

	sweet::ProgramArguments programArguments;

	PlaneDataConfig planeDataConfig;

	PlaneOperators ops;

	class SimPlaneData
	{
	public:
		PlaneData_Spectral prog_h;
		PlaneData_Spectral prog_h_t0;	// at t0
		PlaneData_Spectral prog_u;
		PlaneData_Spectral prog_v;


		void setup(PlaneDataConfig &i_planeDataConfig)
		{
			prog_h.setup(i_planeDataConfig);
			prog_h_t0.setup(i_planeDataConfig);
			prog_u.setup(i_planeDataConfig);
			prog_v.setup(i_planeDataConfig);
		}

		void clear()
		{
			prog_h.clear();
			prog_h_t0.clear();
			prog_u.clear();
			prog_v.clear();
		}
	};

	SimPlaneData simPlaneData;

	PDEAdvPlaneTimeSteppers timeSteppers;


#if SWEET_GUI
	PlaneData_Physical viz_plane_data;

	int render_primitive_id = 0;

	bool gui_active = false;
#endif

	PDEAdvectionPlaneBenchmarksCombined planeBenchmarksCombined;

	double max_error_h0 = -1;

	ShackPlaneDataOps *shackPlaneDataOps;

	ShackPDEAdvectionPlane *pdeAdvectionPlane;
	ShackPDEAdvectionPlaneBenchmarks *pdeAdvectionPlaneBenchmarks;

//	ShackDiagnostics *diagnostics;
	Misc *misc;
	ShackTimestepControl *timestepControl;

	sweet::ShackDictionary shackDict;



public:
	SimulationInstance()	:
		shackPlaneDataOps(nullptr),
		pdeAdvectionPlane(nullptr),
//		diagnostics(nullptr),
		misc(nullptr),
		timestepControl(nullptr),
		prog_argc(0),
		prog_argv(nullptr)
	{
	}


	void clear()
	{
		programArguments.clear();
		shackDict.clear();

		simPlaneData.clear();
		ops.clear();

#if SWEET_GUI
		viz_plane_data.clear();
#endif
	}

	/*
	 * Chekc if help should be printed a do so
	 */
	void checkAndPrintHelp()
	{
		/*
		 * First, check for --help or -h
		 */
		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			shackDict.printProgramArguments();
		}
	}

	/*
	 * Setup main
	 */
	int prog_argc;
	const char *const * prog_argv;

	bool setup(int i_argc, const char *const * i_argv)
	{
		prog_argc = i_argc;
		prog_argv = i_argv;

		/*
		 * Parse program arguments
		 */
		if (!programArguments.setup(prog_argc, prog_argv))
			return error.forwardWithPositiveReturn(programArguments.error);

		/*
		 * SHACK: Register classes which we require
		 */
		pdeAdvectionPlane = shackDict.getAutoRegistration<ShackPDEAdvectionPlane>();
		shackPlaneDataOps = shackDict.getAutoRegistration<ShackPlaneDataOps>();
		//diagnostics = shackDict.getAutoRegistration<ShackDiagnostics>();
		timestepControl = shackDict.getAutoRegistration<ShackTimestepControl>();
		misc = shackDict.getAutoRegistration<Misc>();

		if (shackDict.error.exists())
			return error.forwardWithPositiveReturn(shackDict.error);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		if (!planeBenchmarksCombined.shackRegistration(shackDict))
			return error.forwardWithPositiveReturn(planeBenchmarksCombined.error);

		if (!timeSteppers.shackRegistration(shackDict))
			return error.forwardWithPositiveReturn(timeSteppers.error);


		/*
		 * First, check for --help or -h
		 */
		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			shackDict.printProgramArguments();
			return false;
		}

		/*
		 * SHACK: Process arguments
		 */
		if (!shackDict.processProgramArguments(programArguments))
			return error.forwardWithPositiveReturn(shackDict.error);

		shackDict.printShackData();

#if SWEET_GUI
		/*
		 * Forward information about GUI whether it's enabled
		 */
		gui_active = misc->gui_enabled;
#endif

		/*
		 * Do some validation
		 */
		if (!timestepControl->validateTimestepSize())
			return error.forwardWithPositiveReturn(timestepControl->error);

		simPlaneData.prog_h_t0 = simPlaneData.prog_h;

		if (!timeSteppers.setup(ops, shackDict))
			return error.forwardWithPositiveReturn(timeSteppers.error);

		/*
		 * Setup Plane Data Config & Operators
		 */
		if (!planeDataConfig.setupAuto(
				shackPlaneDataOps->space_res_physical,
				shackPlaneDataOps->space_res_spectral,
				misc->reuse_spectral_transformation_plans
			))
			return error.forwardWithPositiveReturn(planeDataConfig.error);

		if (!ops.setup(
				&planeDataConfig,
				shackPlaneDataOps
			))
			return error.forwardWithPositiveReturn(ops.error);

		simPlaneData.setup(planeDataConfig);

#if SWEET_GUI
		viz_plane_data.setup(planeDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackDict.printShackData();

		if (!planeBenchmarksCombined.setupInitialConditions(
				simPlaneData.prog_h,
				simPlaneData.prog_u,
				simPlaneData.prog_v,
				ops,
				shackDict,
				programArguments
			))
			return error.forwardWithPositiveReturn(planeBenchmarksCombined.error);

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackDict.registrationFinished();
		shackDict.getFinished();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		if (!programArguments.checkAllArgumentsProcessed())
			return error.forwardWithPositiveReturn(programArguments.error);

		return true;
	}

	bool reset()
	{
		clear();

		return setup(prog_argc, prog_argv);
	}

	void printErrors()
	{
		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Lmax error: " << (simPlaneData.prog_h_t0-simPlaneData.prog_h).toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << "RMS error: " << (simPlaneData.prog_h_t0-simPlaneData.prog_h).toPhys().physical_reduce_rms() << std::endl;
	}

	~SimulationInstance()
	{
	}

	void run_timestep()
	{
		if (timestepControl->current_simulation_time + timestepControl->current_timestep_size > timestepControl->max_simulation_time)
			timestepControl->current_timestep_size = timestepControl->max_simulation_time - timestepControl->current_simulation_time;

		timeSteppers.master->run_timestep(
				simPlaneData.prog_h, simPlaneData.prog_u, simPlaneData.prog_v,
				timestepControl->current_timestep_size,
				timestepControl->current_simulation_time
			);

		double dt = timestepControl->current_timestep_size;

		// advance in time
		timestepControl->current_simulation_time += dt;
		timestepControl->current_timestep_nr++;

		if (misc->verbosity > 2)
			std::cout << "timestep: " << timestepControl->current_timestep_nr << ": t=" << timestepControl->current_simulation_time << std::endl;

		max_error_h0 = (simPlaneData.prog_h-simPlaneData.prog_h_t0).toPhys().physical_reduce_max_abs();
	}



	void compute_error()
	{
	}



	bool should_quit()
	{
		if (timestepControl->max_timesteps_nr != -1 && timestepControl->max_timesteps_nr <= timestepControl->current_timestep_nr)
			return true;

		double diff = std::abs(timestepControl->max_simulation_time - timestepControl->current_simulation_time);

		if (	timestepControl->max_simulation_time != -1 &&
				(
						timestepControl->max_simulation_time <= timestepControl->current_simulation_time	||
						diff/timestepControl->max_simulation_time < 1e-11	// avoid numerical issues in time stepping if current time step is 1e-14 smaller than max time step
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
		if (timestepControl->run_simulation_timesteps)
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

		int id = misc->vis_id % 3;
		switch (id)
		{
		case 0:
			sweet::Convert_PlaneDataSpectral_To_PlaneDataPhysical::convert(simPlaneData.prog_h, viz_plane_data);
			break;

		case 1:
			sweet::Convert_PlaneDataSpectral_To_PlaneDataPhysical::convert(simPlaneData.prog_u, viz_plane_data);
			break;

		case 2:
			sweet::Convert_PlaneDataSpectral_To_PlaneDataPhysical::convert(simPlaneData.prog_v, viz_plane_data);
			break;
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 1;
	}


	const char* vis_get_status_string()
	{
		const char* description = "";
		int id = misc->vis_id % 3;

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

		sprintf(title_string,
#if SWEET_MPI
				"Rank %i - "
#endif
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, MaxVal: %.6e, MinVal: %.6e ",
#if SWEET_MPI
				-1,	// TODO: mpi_rank,
#endif
				timestepControl->current_simulation_time,
				timestepControl->current_simulation_time/(60.0*60.0*24.0),
				timestepControl->current_timestep_nr,
				timestepControl->current_timestep_size,
				description,
				viz_plane_data.physical_reduce_max(),
				viz_plane_data.physical_reduce_min()
		);

		return title_string;
	}


	void vis_pause()
	{
		timestepControl->run_simulation_timesteps = !timestepControl->run_simulation_timesteps;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			misc->vis_id++;
			break;

		case 'V':
			misc->vis_id--;
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
	SimulationInstance simulation;
	if (simulation.error.exists())
	{
		std::cerr << "ERROR: " << simulation.error.get() << std::endl;
		return EXIT_FAILURE;
	}

	if (!simulation.setup(i_argc, i_argv))
	{
		if (simulation.error.exists())
			std::cerr << "ERROR: " << simulation.error.get() << std::endl;
		return EXIT_FAILURE;
	}


#if SWEET_GUI
	if (simulation.gui_active)
	{
		std::cout << "GUI" << std::endl;

		VisSweet<SimulationInstance> visSweet(&simulation);
		std::cout << "Max error h0: "<< simulation.max_error_h0 << std::endl;
	}
	else
#endif
	{
		if (!simulation.timestepControl->validateMaxSimulationTime())
		{
			simulation.timestepControl->error.print();
			exit(1);
		}

		while (!simulation.should_quit())
			simulation.run_timestep();

		simulation.printErrors();
	}

	std::cout << "FIN" << std::endl;
	return 0;
}
