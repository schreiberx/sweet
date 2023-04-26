/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_XBRAID_ODE_SCALAR_PROGRAMXBRAIDODESCALAR_HPP_
#define SRC_PROGRAMS_XBRAID_ODE_SCALAR_PROGRAMXBRAIDODESCALAR_HPP_


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackIOData.hpp>

#include <sweet/xbraid/ShackXBraid.hpp>
#include <sweet/xbraid/XBraid_sweet_lib.hpp>

// Benchmarks
#include "ODEScalarBenchmarksCombined.hpp"

// Time steppers
#include "ODEScalarTimeSteppers.hpp"

class ProgramXBraidODEScalar
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

	////// time integrators
	////ODEScalarTimeSteppers timeSteppers;

	////// Handler to all benchmarks
	////ODEScalarBenchmarksCombined scalarBenchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackODEScalarTimeDiscretization *shackTimeDisc;
	sweet::ShackXBraid *shackXBraid;

	// XBraid
	sweet_BraidApp* xbraid_app = nullptr;
	BraidCore* xbraid_core = nullptr;

	// MPI
	MPI_Comm mpi_comm;
	int mpi_rank;

public:
	ProgramXBraidODEScalar(
			int i_argc,
			char *const * const i_argv,
			MPI_Comm i_mpi_comm,
			int i_mpi_rank
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		mpi_comm(i_mpi_comm),
		mpi_rank(i_mpi_rank)
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
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackODEScalarTimeDiscretization>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::ShackXBraid>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		////////scalarBenchmarksCombined.shackRegistration(shackProgArgDict);
		////////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(scalarBenchmarksCombined);

		////////timeSteppers.shackRegistration(shackProgArgDict);
		////////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;
		shackXBraid = nullptr;

		/////scalarBenchmarksCombined.clear();
		/////timeSteppers.clear();
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

		/////////*
		//////// * Setup the time steppers and their buffers
		//////// */
		////////timeSteppers.setup(shackProgArgDict);
		////////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		//////////scalarBenchmarksCombined.setupInitialConditions(
		//////////		data.prog_u,
		//////////		shackProgArgDict
		//////////	);
		//////////ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(scalarBenchmarksCombined);

		//////////////data.prog_u_t0 = data.prog_u;

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

		return true;
	}

	void clear_3_data()
	{

		/////////timeSteppers.clear();

		data.clear();
	}

	bool setup_4_xbraid()
	{

		// get the number of timesteps in the finest level
		int nt = (int) (shackTimestepControl->max_simulation_time / shackTimestepControl->current_timestep_size);
		if (nt * shackTimestepControl->current_timestep_size < shackTimestepControl->max_simulation_time - 1e-10)
			nt++;

		// XBraid app (user-defined)
		this->xbraid_app = new sweet_BraidApp(this->mpi_comm, this->mpi_rank, 0., shackTimestepControl->max_simulation_time, nt, &shackProgArgDict);

		// XBraid core
		this->xbraid_core = new BraidCore(this->mpi_comm, this->xbraid_app);

		this->xbraid_app->setup(*this->xbraid_core);
	}

	void clear_4_xbraid()
	{

		if (this->xbraid_core)
		{
			delete this->xbraid_core;
			this->xbraid_core = nullptr;
		}

		if (this->xbraid_app)
		{
			delete this->xbraid_app;
			this->xbraid_app = nullptr;
		}

	}

	bool setup()
	{
		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_data())
			return false;

		if (!setup_4_xbraid())
			return false;

		std::cout << "SETUP FINISHED" << std::endl;
		return true;
	}
	void clear()
	{
		clear_4_xbraid();
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

	~ProgramXBraidODEScalar()
	{
		clear();
	}


	bool runXBraid()
	{

		shackTimestepControl->timestepHelperStart();

		// Run Simulation
		this->xbraid_core->Drive();

		shackTimestepControl->timestepHelperEnd();

		return true;

	}

	bool should_quit()
	{
		////////////return shackTimestepControl->isFinalTimestepReached();
	}

};




#endif
