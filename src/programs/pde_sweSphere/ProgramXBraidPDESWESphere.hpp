/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_XBRAID_PDE_SWESPHERE_PROGRAMXBRAIDPDESWESPHERE_HPP_
#define SRC_PROGRAMS_XBRAID_PDE_SWESPHERE_PROGRAMXBRAIDPDESWESPHERE_HPP_


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Include everything we need for simulations on the plane
#include <sweet/core/sphere/Sphere.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackSphereDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include "benchmarks/ShackPDESWESphereBenchmarks.hpp"

#include <sweet/xbraid/ShackXBraid.hpp>
#include <sweet/xbraid/XBraid_sweet_lib.hpp>

// Benchmarks
#include "PDESWESphere_BenchmarksCombined.hpp"

// Time steppers
#include "time/PDESWESphere_TimeSteppers.hpp"

#include<vector>

class ProgramXBraidPDESWESphere
{
public:
	sweet::ErrorBase error;

	/*
	 * Just a class to store simulation data all together
	 */
	class DataConfigOps
	{
	public:
		sweet::ErrorBase error;

		sweet::SphereData_Config sphereDataConfig;
		sweet::SphereOperators ops;

		sweet::SphereData_Spectral prog_phi_pert;
		sweet::SphereData_Spectral prog_div;
		sweet::SphereData_Spectral prog_vrt;


		sweet::SphereData_Spectral t0_prog_phi_pert;
		sweet::SphereData_Spectral t0_prog_div;
		sweet::SphereData_Spectral t0_prog_vrt;

		bool setup(
				sweet::ShackSphereDataOps *i_shackSphereDataOps,
				bool i_setup_spectral_transforms = true		// for reset()
		)
		{
			/*
			 * Setup Sphere Data Config & Operators
			 */
			if (i_setup_spectral_transforms)
			{
				sphereDataConfig.setupAuto(i_shackSphereDataOps);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereDataConfig);
			}

			ops.setup(&sphereDataConfig, i_shackSphereDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			prog_phi_pert.setup(sphereDataConfig);
			prog_div.setup(sphereDataConfig);
			prog_vrt.setup(sphereDataConfig);

			return true;
		}

		void clear(bool i_clear_spectral_transforms = true)
		{
			prog_phi_pert.clear();
			prog_div.clear();
			prog_vrt.clear();

			t0_prog_phi_pert.clear();
			t0_prog_div.clear();
			t0_prog_vrt.clear();

			ops.clear();

			if (i_clear_spectral_transforms)
				sphereDataConfig.clear();
		}
	};


	// Simulation data
	DataConfigOps dataConfigOps;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWESphereTimeDiscretization *shackTimeDisc;
	sweet::ShackParallelization *shackParallelization;
	ShackPDESWESphere *shackPDESWESphere;
	ShackPDESWESphereBenchmarks *shackBenchmarks;
	sweet::ShackXBraid *shackXBraid;

	// XBraid
	sweet_BraidApp* xbraid_app = nullptr;
	BraidCore* xbraid_core = nullptr;

	// MPI
	MPI_Comm mpi_comm;
	int mpi_rank;

public:
	ProgramXBraidPDESWESphere(
			int i_argc,
			char *const * const i_argv,
			MPI_Comm i_mpi_comm,
			int i_mpi_rank
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphereDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackParallelization(nullptr),
		shackPDESWESphere(nullptr),
		shackBenchmarks(nullptr),
		shackXBraid(nullptr),
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
		shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::ShackParallelization>();
		shackPDESWESphere = shackProgArgDict.getAutoRegistration<ShackPDESWESphere>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDESWESphereBenchmarks>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::ShackXBraid>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackSphereDataOps = nullptr;
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackTimeDisc = nullptr;
		shackParallelization = nullptr;
		shackPDESWESphere = nullptr;
		shackBenchmarks = nullptr;
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

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		/*
		 * Setup Sphere Data Config & Operators
		 */
		dataConfigOps.setup(shackSphereDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

		// get the number of timesteps in the finest level
		int nt = (int) (shackTimestepControl->max_simulation_time / shackTimestepControl->current_timestep_size);
		if (nt * shackTimestepControl->current_timestep_size < shackTimestepControl->max_simulation_time - 1e-10)
			nt++;

		// XBraid app (user-defined)
		///this->xbraid_app = new sweet_BraidApp(this->mpi_comm, this->mpi_rank, 0., shackTimestepControl->max_simulation_time, nt);/////, planeDataConfigs, ops);//, &shackProgArgDict);
		this->xbraid_app = new sweet_BraidApp(this->mpi_comm, this->mpi_rank, 0., shackTimestepControl->max_simulation_time, nt, &dataConfigOps.sphereDataConfig, &dataConfigOps.ops);
		this->xbraid_app->shackRegistration(shackProgArgDict);

		// XBraid core
		if (shackXBraid->xbraid_run_wrapper_tests)
			this->xbraid_app->setup();
		else
		{
			this->xbraid_core = new BraidCore(this->mpi_comm, this->xbraid_app);
			this->xbraid_app->setup(*this->xbraid_core);
		}

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

		dataConfigOps.clear();
	}

	bool setup()
	{
		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_data())
			return false;

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

	////void printSimulationErrors()
	////{
	////	std::cout << "Error compared to initial condition" << std::endl;
	////	std::cout << "Error: " << std::abs(data.prog_u_t0-data.prog_u) << std::endl;
	////}

	~ProgramXBraidPDESWESphere()
	{
		clear();
	}


	bool runXBraid()
	{

		shackTimestepControl->timestepHelperStart();

		// Run wrapper tests
		if (shackXBraid->xbraid_run_wrapper_tests)
		{
			BraidUtil braid_util;
			int test = braid_util.TestAll(this->xbraid_app, this->mpi_comm, stdout, 0., shackTimestepControl->current_timestep_size, shackTimestepControl->current_timestep_size * 2);
			if (test == 0)
				SWEETError("Tests failed!");
			else
				std::cout << "Tests successful!" << std::endl;
		}
		else
		{
			// Run Simulation
			this->xbraid_core->Drive();
		}

		shackTimestepControl->timestepHelperEnd();

		return true;

	}

	bool should_quit()
	{
		////////////return shackTimestepControl->isFinalTimestepReached();
	}

};




#endif
