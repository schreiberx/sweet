/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_PROGRAMODESCALAR_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_PROGRAMODESCALAR_HPP_


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackIOData.hpp>

#include <sweet/core/shacksShared/ShackParallelization.hpp>

// Benchmarks
#include "ODEScalarBenchmarksCombined.hpp"

// Time steppers
#include "ODEScalarTimeSteppers.hpp"

#include "ODEScalar_FileOutput.hpp"

class ProgramODEScalar
#if SWEET_GUI
		:	public SimulationGUICallbacks
#endif
{
public:
	sweet::ErrorBase error;

	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::ErrorBase error;

		double prog_u;
		double prog_u_t0;

		bool setup()
		{
			return true;
		}

		void clear()
		{
		}
	};

	// Simulation data
	Data data;

	// time integrators
	ODEScalarTimeSteppers timeSteppers;

	// Handler to all benchmarks
	ODEScalarBenchmarksCombined scalarBenchmarksCombined;

	ODEScalar_FileOutput fileOutput;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackParallelization *shackParallelization;
	ShackODEScalar *shackODEScalar;
	ShackODEScalarTimeDiscretization *shackTimeDisc;
	ShackODEScalarBenchmarks *shackBenchmarks;

	int timestep_nr_last_output_simtime = -1;

public:
	ProgramODEScalar(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackParallelization(nullptr),
		shackODEScalar(nullptr),
		shackTimeDisc(nullptr),
		shackBenchmarks(nullptr)
	{
		ERROR_CHECK_COND_RETURN(shackProgArgDict);
	}



	bool setup_1_shackRegistration()
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::ShackParallelization>();
		shackODEScalar = shackProgArgDict.getAutoRegistration<ShackODEScalar>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackODEScalarTimeDiscretization>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackODEScalarBenchmarks>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		scalarBenchmarksCombined.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(scalarBenchmarksCombined);

		timeSteppers.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackODEScalar = nullptr;
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackParallelization = nullptr;
		shackTimeDisc = nullptr;
		shackBenchmarks = nullptr;

		scalarBenchmarksCombined.clear();
		timeSteppers.clear();
		shackProgArgDict.clear();
	}

	bool setup_2_processArguments()
	{
		shackProgArgDict.setup();

		/*
		 * SHACK: Process arguments
		 */
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Do some validation of program arguments
		 */
		shackTimestepControl->validateTimestepSize();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);

		return true;
	}

	void clear_2_process_arguments()
	{
		shackProgArgDict.clear();
	}

	bool setup_3_data()
	{

		/*
		 * Setup the time steppers and their buffers
		 */
		timeSteppers.setup(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		scalarBenchmarksCombined.setupInitialConditions(
				data.prog_u,
				shackProgArgDict
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(scalarBenchmarksCombined);

		data.prog_u_t0 = data.prog_u;

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackProgArgDict.closeRegistration();
		shackProgArgDict.closeGet();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		shackProgArgDict.checkAllArgumentsProcessed();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		fileOutput.setup(shackIOData, shackTimestepControl, shackODEScalar);

		return true;
	}
	void clear_3_data()
	{
#if SWEET_GUI
		vis_plane_data.clear();
#endif

		fileOutput.clear();

		timeSteppers.clear();

		data.clear();
	}

	bool setup()
	{
		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_data())
			return false;

		/*
		 * Output data for the first time step as well if output of datafiels is requested
		 */
		if (shackIOData->output_each_sim_seconds >= 0)
			_timestepDoOutput();

		std::cout << "SETUP FINISHED" << std::endl;
		return true;
	}
	void clear()
	{
		clear_3_data();
		clear_2_process_arguments();
		clear_1_shackRegistration();
	}

	bool reset()
	{
		clear();

		if (!setup())
		{
			error.print();
			return false;
		}

		return !error.exists();
	}

	void printSimulationErrors()
	{
		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Error: " << std::abs(data.prog_u_t0-data.prog_u) << std::endl;
	}

	~ProgramODEScalar()
	{
		clear();
	}

	bool runTimestep()
	{
		shackTimestepControl->timestepHelperStart();

		timeSteppers.master->runTimestep(
				data.prog_u,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);

		shackTimestepControl->timestepHelperEnd();

		if (shackIOData->verbosity > 2)
		{
			double error = std::abs(data.prog_u_t0-data.prog_u);
			std::cout << "timestep: " << shackTimestepControl->current_timestep_nr << ": dt=" << shackTimestepControl->current_timestep_size << ": t=" << shackTimestepControl->current_simulation_time << std::endl;
			std::cout << "error:" << error << std::endl;
		}
		return true;
	}


	void _timestepDoOutput()
	{

		if (shackParallelization->isMPIRoot)
		{
			fileOutput.write_file_output(
					data.prog_u
			);
		}

		if (shackIOData->output_each_sim_seconds > 0)
			while (shackIOData->output_next_sim_seconds <= shackTimestepControl->current_simulation_time)
				shackIOData->output_next_sim_seconds += shackIOData->output_each_sim_seconds;
	}

public:
	bool timestepHandleOutput()
	{
		if (shackIOData->output_each_sim_seconds < 0)
			return false;

		if (shackTimestepControl->current_simulation_time == timestep_nr_last_output_simtime)
			return false;

		timestep_nr_last_output_simtime = shackTimestepControl->current_simulation_time;

		if (shackTimestepControl->current_simulation_time < shackTimestepControl->max_simulation_time - shackIOData->output_each_sim_seconds*1e-10)
		{
			if (shackIOData->output_next_sim_seconds > shackTimestepControl->current_simulation_time)
				return false;
		}

		if (shackParallelization->isMPIRoot)
			if (shackIOData->verbosity > 0)
				std::cout << std::endl;

		_timestepDoOutput();

		return true;
	}

	bool should_quit()
	{
		return shackTimestepControl->isFinalTimestepReached();
	}

};




#endif
