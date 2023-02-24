/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pdeAdvectionPlane/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pdeAdvectionPlane/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pdeAdvectionPlane/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/ErrorBase.hpp>

// Parse program arguments
#include <sweet/ProgramArguments.hpp>

// Include everything we need for simulations on the plane
#include <sweet/plane/Plane.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/shacks/ShackDictionary.hpp>

// Different shacks we need in this file
#include <sweet/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/shacksShared/ShackIOData.hpp>

// Benchmarks
#include "pdeAdvectionPlane/benchmarks/PDEAdvectionPlaneBenchmarksCombined.hpp"

// Time steppers
#include "pdeAdvectionPlane/time/PDEAdvPlaneTimeSteppers.hpp"

#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif


class SimulationPDEAdvectionPlane
#if SWEET_GUI
		:	public SimulationGUICallbacks
#endif
{
public:
	sweet::ErrorBase error;
	sweet::ProgramArguments programArguments;

	/*
	 * Just a class to store simulation data all together
	 */
	class SimPlaneData
	{
	public:
		sweet::ErrorBase error;

		sweet::PlaneDataConfig planeDataConfig;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral prog_h;
		sweet::PlaneData_Spectral prog_h_t0;	// at t0
		sweet::PlaneData_Spectral prog_u;
		sweet::PlaneData_Spectral prog_v;


		bool setup(sweet::ShackPlaneDataOps *i_shackPlaneDataOps)
		{
			/*
			 * Setup Plane Data Config & Operators
			 */
			if (!planeDataConfig.setupAuto(*i_shackPlaneDataOps))
				return error.forwardWithPositiveReturn(planeDataConfig.error);

			if (!ops.setup(
					planeDataConfig,
					*i_shackPlaneDataOps
				))
				return error.forwardWithPositiveReturn(ops.error);

			prog_h.setup(planeDataConfig);
			prog_h_t0.setup(planeDataConfig);
			prog_u.setup(planeDataConfig);
			prog_v.setup(planeDataConfig);

			return true;
		}

		void clear()
		{
			prog_h.clear();
			prog_h_t0.clear();
			prog_u.clear();
			prog_v.clear();

			ops.clear();
			planeDataConfig.clear();
		}
	};

	// Simulation data
	SimPlaneData simPlaneData;

	// time integrators
	PDEAdvPlaneTimeSteppers timeSteppers;


	// Handler to all benchmarks
	PDEAdvectionPlaneBenchmarksCombined planeBenchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackDictionary shackDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;


#if SWEET_GUI
	// Data to visualize is stored to this variable
	sweet::PlaneData_Physical viz_plane_data;

	// Which primitive to use for rendering
	int viz_render_type_of_primitive_id = 0;

	// Which primitive to use for rendering
	int viz_data_id = 0;

#endif

public:
	SimulationPDEAdvectionPlane()	:
		shackPlaneDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		prog_argc(0),
		prog_argv(nullptr)
	{
	}


	/*
	 * Clear all data
	 */
	void clear()
	{
		programArguments.clear();
		shackDict.clear();

		simPlaneData.clear();

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
		programArguments.setup(prog_argc, prog_argv);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(programArguments);

		/*
		 * SHACK: Register classes which we require
		 */
		shackPlaneDataOps = shackDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackTimestepControl = shackDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackIOData = shackDict.getAutoRegistration<sweet::ShackIOData>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		planeBenchmarksCombined.shackRegistration(shackDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(planeBenchmarksCombined);

		timeSteppers.shackRegistration(shackDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(planeBenchmarksCombined);

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
		shackDict.processProgramArguments(programArguments);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackDict);


		/*
		 * Do some validation of program arguments
		 */
		shackTimestepControl->validateTimestepSize();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*shackTimestepControl);

		/*
		 * Setup Plane Data Config & Operators
		 */
		simPlaneData.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(simPlaneData);


		timeSteppers.setup(shackDict, simPlaneData.ops);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(timeSteppers);


#if SWEET_GUI
		viz_plane_data.setup(simPlaneData.planeDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackDict.printShackData();

		planeBenchmarksCombined.setupInitialConditions(
				simPlaneData.prog_h,
				simPlaneData.prog_u,
				simPlaneData.prog_v,
				simPlaneData.ops,
				shackDict,
				programArguments
			);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(planeBenchmarksCombined);

		simPlaneData.prog_h_t0 = simPlaneData.prog_h;

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackDict.closeRegistration();
		shackDict.closeGet();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		programArguments.checkAllArgumentsProcessed();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(programArguments);

		return true;
	}

	bool reset()
	{
		clear();

		setup(prog_argc, prog_argv);

		return !error.exists();
	}

	void printSimulationErrors()
	{
		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Lmax error: " << (simPlaneData.prog_h_t0-simPlaneData.prog_h).toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << "RMS error: " << (simPlaneData.prog_h_t0-simPlaneData.prog_h).toPhys().physical_reduce_rms() << std::endl;
	}

	~SimulationPDEAdvectionPlane()
	{
		clear();
	}

	bool run_timestep()
	{
		shackTimestepControl->timestepHelperStart();

		timeSteppers.master->run_timestep(
				simPlaneData.prog_h, simPlaneData.prog_u, simPlaneData.prog_v,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);

		shackTimestepControl->timestepHelperEnd();

		if (shackIOData->verbosity > 2)
			std::cout << "timestep: " << shackTimestepControl->current_timestep_nr << ": dt=" << shackTimestepControl->current_timestep_size << ": t=" << shackTimestepControl->current_simulation_time << std::endl;

		return true;
	}


	bool should_quit()
	{
		return shackTimestepControl->isFinalTimestepReached();
	}


#if SWEET_GUI
	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (shackTimestepControl->run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				run_timestep();
	}

	void vis_get_vis_data_array(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		*o_render_primitive_id = viz_render_type_of_primitive_id;

		int id = viz_data_id % 3;
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
		int id = viz_data_id % 3;

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
				shackTimestepControl->current_simulation_time,
				shackTimestepControl->current_simulation_time/(60.0*60.0*24.0),
				shackTimestepControl->current_timestep_nr,
				shackTimestepControl->current_timestep_size,
				description,
				viz_plane_data.physical_reduce_max(),
				viz_plane_data.physical_reduce_min()
		);

		return title_string;
	}


	void vis_pause()
	{
		shackTimestepControl->run_simulation_timesteps = !shackTimestepControl->run_simulation_timesteps;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			viz_data_id++;
			break;

		case 'V':
			viz_data_id--;
			break;

		case 'b':
			viz_render_type_of_primitive_id = (viz_render_type_of_primitive_id + 1) % 2;
			break;
		}
	}
#endif
};




int main(int i_argc, char *i_argv[])
{
	SimulationPDEAdvectionPlane simulation;
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
	if (simulation.shackIOData->gui_enabled)
	{
		VisSweet visSweet(simulation);
	}
	else
#endif
	{
		if (!(	simulation.shackTimestepControl->validateMaxSimulationTime() ||
				simulation.shackTimestepControl->validateMaxTimestepNr()
		))
		{
			simulation.shackTimestepControl->error.print();
			return EXIT_FAILURE;
		}

		while (!simulation.should_quit())
			simulation.run_timestep();
	}

	simulation.printSimulationErrors();

	std::cout << "FIN" << std::endl;
	return 0;
}
