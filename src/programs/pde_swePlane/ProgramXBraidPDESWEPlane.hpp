/*
 * 		Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_XBRAID_PDE_SWEPLANE_PROGRAMXBRAIDPDESWEPLANE_HPP_
#define SRC_PROGRAMS_XBRAID_PDE_SWEPLANE_PROGRAMXBRAIDPDESWEPLANE_HPP_


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Include everything we need for simulations on the plane
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/plane/PlaneData_Config.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/plane/PlaneDataGridMapping.hpp>
#include "ShackPDESWEPlane_Diagnostics.hpp"
#include "benchmarks/ShackPDESWEPlaneBenchmarks.hpp"

#include <sweet/xbraid/ShackXBraid.hpp>
#include <sweet/xbraid/XBraid_sweet_lib.hpp>

// Benchmarks
#include "PDESWEPlane_BenchmarksCombined.hpp"

// Time steppers
#include "PDESWEPlane_TimeSteppers.hpp"

#include<vector>

class ProgramXBraidPDESWEPlane
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

		sweet::PlaneData_Config planeDataConfig;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral prog_h_pert;
		sweet::PlaneData_Spectral prog_u;
		sweet::PlaneData_Spectral prog_v;

		// TODO: Get rid of me right here
		// Initial values for comparison with analytical solution
		sweet::PlaneData_Spectral t0_prog_h_pert;
		sweet::PlaneData_Spectral t0_prog_u;
		sweet::PlaneData_Spectral t0_prog_v;

		// Mapping between grids
		sweet::PlaneDataGridMapping gridMapping;
		
		// Diagnostics measures
		int last_timestep_nr_update_diagnostics = -1;
			
		bool setup(sweet::ShackPlaneDataOps *i_shackPlaneDataOps)
		{
			/*
			 * Setup Plane Data Config & Operators
			 */
			planeDataConfig.setupAuto(*i_shackPlaneDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeDataConfig);

			ops.setup(planeDataConfig, *i_shackPlaneDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			prog_h_pert.setup(planeDataConfig);
			prog_u.setup(planeDataConfig);
			prog_v.setup(planeDataConfig);

			t0_prog_h_pert.setup(planeDataConfig);
			t0_prog_u.setup(planeDataConfig);
			t0_prog_v.setup(planeDataConfig);

			last_timestep_nr_update_diagnostics = -1;
			
			if (i_shackPlaneDataOps->space_grid_use_c_staggering)
				gridMapping.setup(i_shackPlaneDataOps, &planeDataConfig);

			return true;
		}

		void clear()
		{
			prog_h_pert.clear();
			prog_u.clear();
			prog_v.clear();
			
			t0_prog_h_pert.clear();
			t0_prog_u.clear();
			t0_prog_v.clear();

			ops.clear();
			planeDataConfig.clear();
		}

	};

	// Simulation data
	Data dataAndOps;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWEPlaneTimeDiscretization *shackTimeDisc;
	sweet::ShackParallelization *shackParallelization;
	ShackPDESWEPlane *shackPDESWEPlane;
	ShackPDESWEPlaneBenchmarks *shackBenchmarks;
	sweet::ShackXBraid *shackXBraid;

	// XBraid
	sweet_BraidApp* xbraid_app = nullptr;
	BraidCore* xbraid_core = nullptr;

	// MPI
	MPI_Comm mpi_comm;
	int mpi_rank;

	//DataConfig and Ops for each level
	std::vector<sweet::PlaneData_Config*> planeDataConfigs = {};
	std::vector<sweet::PlaneOperators*> ops = {};

public:
	ProgramXBraidPDESWEPlane(
			int i_argc,
			char *const * const i_argv,
			MPI_Comm i_mpi_comm,
			int i_mpi_rank
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackPlaneDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackParallelization(nullptr),
		shackPDESWEPlane(nullptr),
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
		shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneTimeDiscretization>();
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::ShackParallelization>();
		shackPDESWEPlane = shackProgArgDict.getAutoRegistration<ShackPDESWEPlane>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneBenchmarks>();
		shackXBraid = shackProgArgDict.getAutoRegistration<sweet::ShackXBraid>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackPlaneDataOps = nullptr;
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackTimeDisc = nullptr;
		shackParallelization = nullptr;
		shackPDESWEPlane = nullptr;
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
		 * Setup Plane Data Config & Operators
		 */
		dataAndOps.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataAndOps);

		/////////////////////////////
		///////////// SETUP XBRAID //
		/////////////////////////////

		///////////// Set planeDataConfig and planeOperators for each level
		//////////////sweet::PlaneOperators op(dataAndOps.planeDataConfig, this->shackPlaneDataOps->plane_domain_size, this->shackPlaneDataOps->space_use_spectral_basis_diffs);
		///////////for (int i = 0; i < this->shackXBraid->xbraid_max_levels; i++)
		///////////{
		///////////	if (this->shackXBraid->xbraid_spatial_coarsening)
		///////////	{
		///////////		int N_physical[2] = {-1, -1};
		///////////		int N_spectral[2];
		///////////		for (int j = 0; j < 2; j++)
		///////////		{
		///////////			// proportional to time step
		///////////			if (this->shackXBraid->xbraid_spatial_coarsening == 1)
		///////////				N_spectral[j] = std::max(4, int(this->shackPlaneDataOps->space_res_spectral[j] / std::pow(this->shackXBraid->xbraid_cfactor, i)));
		///////////			else if (this->shackXBraid->xbraid_spatial_coarsening > 1)
		///////////			{
		///////////				if (i == 0)
		///////////					N_spectral[j] = std::max(4, this->shackPlaneDataOps->space_res_spectral[j]);
		///////////				else
		///////////					N_spectral[j] = std::max(4, this->shackXBraid->xbraid_spatial_coarsening);
		///////////			}
		///////////			else
		///////////				SWEETError("Invalid parameter xbraid_spatial_coarsening");
		///////////		}
		///////////		planeDataConfigs.push_back(new sweet::PlaneData_Config);
		///////////		planeDataConfigs.back()->setupAuto(N_physical, N_spectral, this->shackPlaneDataOps->reuse_spectral_transformation_plans);

		///////////		//PlaneOperators op_level(planeDataConfigs.back(), simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);
		///////////		///ops.push_back(new sweet::PlaneOperators(planeDataConfigs.back(), this->shackPlaneDataOps->plane_domain_size, this->shackPlaneDataOps->space_use_spectral_basis_diffs));
		///////////		ops.push_back(new sweet::PlaneOperators(planeDataConfigs.back(), this->shacksPlaneDataOps_levels[i]));

		////////////////////PlaneData_Config *i_planeDataConfig,		///< data config setup for spectral transformations
		////////////////////ShackPlaneDataOps *i_planeDataOps

		///////////		std::cout << "Spectral resolution at level " << i << " : " << N_spectral[0] << " " << N_spectral[1] << std::endl;
		///////////	}
		///////////	else
		///////////	{
		///////////		planeDataConfigs.push_back(&dataAndOps.planeDataConfig);
		///////////		ops.push_back(&dataAndOps.ops);
		///////////	}
		///////////}


		// get the number of timesteps in the finest level
		int nt = (int) (shackTimestepControl->max_simulation_time / shackTimestepControl->current_timestep_size);
		if (nt * shackTimestepControl->current_timestep_size < shackTimestepControl->max_simulation_time - 1e-10)
			nt++;

		// XBraid app (user-defined)
		///this->xbraid_app = new sweet_BraidApp(this->mpi_comm, this->mpi_rank, 0., shackTimestepControl->max_simulation_time, nt);/////, planeDataConfigs, ops);//, &shackProgArgDict);
		this->xbraid_app = new sweet_BraidApp(this->mpi_comm, this->mpi_rank, 0., shackTimestepControl->max_simulation_time, nt, &dataAndOps.planeDataConfig, &dataAndOps.ops);
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

		if (this->shackXBraid->xbraid_spatial_coarsening)
			for (int i = 0; i < this->shackXBraid->xbraid_max_levels; i++)
			{
				delete planeDataConfigs[i];
				delete ops[i];
				planeDataConfigs[i] = nullptr;
				ops[i] = nullptr;
			}

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

		dataAndOps.clear();
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

	~ProgramXBraidPDESWEPlane()
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
