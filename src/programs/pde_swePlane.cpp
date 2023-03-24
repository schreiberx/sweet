/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */



#include "pde_swePlane/ProgramPDESWEPlane.hpp"



#if SWEET_MPI
int mpi_rank;
#endif

bool isMPIRoot()
{
#if SWEET_MPI
	return mpi_rank == 0;
#else
	return true;
#endif
}



#if 0



int main(int i_argc, char *i_argv[])
{
#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

#if SWEET_MPI

	#if SWEET_THREADING_SPACE
		int provided;
		MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE)
			SWEETError("MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye...");
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

	//input parameter names (specific ones for this program)
	const char *user_defined_prog_params[] = {
			"initial-freq-x-mul",		/// frequency multipliers for special scenario setup
			"initial-freq-y-mul",
			nullptr
	};

	// default values for specific input (for general input see shacks/ShackDictionary.hpp)
	//shackDict.bogus.var[0];  //frequency in x for waves test case
	//shackDict.bogus.var[1];  //frequency in y for waves test case

	// Help menu
	if (!shackDict.setupFromMainParameters(i_argc, i_argv, user_defined_prog_params))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--initial-freq-x-mul [float]" << std::endl;
		std::cout << "	--initial-freq-y-mul [float]" << std::endl;
		std::cout << "" << std::endl;


#if SWEET_PARAREAL
		shackDict.parareal.printProgramArguments();
#endif
		return -1;
	}

	if (shackDict.misc.verbosity > 5)
		std::cout << " + Setting up FFT plans..." << std::flush;

	planeDataConfigInstance.setupAuto(shackDict.disc.space_res_physical, shackDict.disc.space_res_spectral, shackDict.misc.reuse_spectral_transformation_plans);

	if (shackDict.misc.verbosity > 5)
		std::cout << " done" << std::endl;

	planeDataConfig->printInformation();

	// Print header
	std::cout << std::endl;
	shackDict.outputConfig();
	std::cout << "Computing error: " << shackDict.misc.computeErrors << std::endl;
	std::cout << std::endl;

	std::ostringstream buf;
	buf << std::setprecision(14);


#if SWEET_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	#if (SWEET_PARAREAL != 2) && (!SWEET_XBRAID)
	// only start simulation and time stepping for first rank
	if (isMPIRoot())
	#endif
#endif
	{
#if SWEET_MPI
		std::cout << "MPI_RANK: " << mpi_rank << std::endl;
#endif

		StopwatchBox::getInstance().main.start();

#if SWEET_PARAREAL
		if (shackDict.parareal.enabled)
		{

			PlaneOperators op(planeDataConfig, shackDict.sim.plane_domain_size, shackDict.disc.space_use_spectral_basis_diffs);

			// Set planeDataConfig and planeOperators for each level
			std::vector<PlaneDataConfig*> planeDataConfigs;
			std::vector<PlaneOperators*> ops;

			// fine
			planeDataConfigs.push_back(planeDataConfig);
			ops.push_back(&op);

			// coarse
			if (shackDict.parareal.spatial_coarsening)
			{
				///for (int j = 0; j < 2; j++)
				///	assert(shackDict.disc.space_res_physical[j] == -1);
				int N_physical[2] = {-1, -1};
				int N_spectral[2];
				double frac;
				if ( shackDict.parareal.coarse_timestep_size > 0)
					frac = shackDict.timecontrol.current_timestep_size / shackDict.parareal.coarse_timestep_size;
				else
					frac = shackDict.timecontrol.current_timestep_size / (shackDict.timecontrol.max_simulation_time / shackDict.parareal.coarse_slices );
				for (int j = 0; j < 2; j++)
					N_spectral[j] = std::max(4, int(shackDict.disc.space_res_spectral[j] * frac));
				planeDataConfigs.push_back(new PlaneDataConfig);
				planeDataConfigs.back()->setupAuto(N_physical, N_spectral, shackDict.misc.reuse_spectral_transformation_plans);

				ops.push_back(new PlaneOperators(planeDataConfigs.back(), shackDict.sim.plane_domain_size, shackDict.disc.space_use_spectral_basis_diffs));
			}
			else
			{
				planeDataConfigs.push_back(planeDataConfig);
				ops.push_back(&op);
			}


			SWE_Plane_TimeSteppers* timeSteppersFine = new SWE_Plane_TimeSteppers;
			SWE_Plane_TimeSteppers* timeSteppersCoarse = new SWE_Plane_TimeSteppers;

			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller<SWE_Plane_TimeSteppers, 3> parareal_Controller(	&shackDict,
												planeDataConfigs,
												ops,
												timeSteppersFine,
												timeSteppersCoarse);

			// setup controller. This initializes several simulation instances
			parareal_Controller.setup();

			// execute the simulation
			parareal_Controller.run();

			delete timeSteppersFine;
			delete timeSteppersCoarse;

			if (shackDict.parareal.spatial_coarsening)
			{
				delete planeDataConfigs[1];
				delete ops[1];
			}
		}
		else
#endif


#if SWEET_XBRAID

		if (shackDict.xbraid.xbraid_enabled)
		{

			PlaneOperators op(planeDataConfig, shackDict.sim.plane_domain_size, shackDict.disc.space_use_spectral_basis_diffs);

			// Set planeDataConfig and planeOperators for each level
			std::vector<PlaneDataConfig*> planeDataConfigs;
			std::vector<PlaneOperators*> ops;
			for (int i = 0; i < shackDict.xbraid.xbraid_max_levels; i++)
			{
				if (shackDict.xbraid.xbraid_spatial_coarsening)
				{
					int N_physical[2] = {-1, -1};
					int N_spectral[2];
					for (int j = 0; j < 2; j++)
					{
						// proportional to time step
						if (shackDict.xbraid.xbraid_spatial_coarsening == 1)
							N_spectral[j] = std::max(4, int(shackDict.disc.space_res_spectral[j] / std::pow(shackDict.xbraid.xbraid_cfactor, i)));
						else if (shackDict.xbraid.xbraid_spatial_coarsening > 1)
						{
							if (i == 0)
								N_spectral[j] = std::max(4, shackDict.disc.space_res_spectral[j]);
							else
								N_spectral[j] = std::max(4, shackDict.xbraid.xbraid_spatial_coarsening);
						}
						else
							SWEETError("Invalid parameter xbraid_spatial_coarsening");
					}
					planeDataConfigs.push_back(new PlaneDataConfig);
					planeDataConfigs.back()->setupAuto(N_physical, N_spectral, shackDict.misc.reuse_spectral_transformation_plans);

					//PlaneOperators op_level(planeDataConfigs.back(), shackDict.sim.plane_domain_size, shackDict.disc.space_use_spectral_basis_diffs);
					ops.push_back(new PlaneOperators(planeDataConfigs.back(), shackDict.sim.plane_domain_size, shackDict.disc.space_use_spectral_basis_diffs));

					std::cout << "Spectral resolution at level " << i << " : " << N_spectral[0] << " " << N_spectral[1] << std::endl;
				}
				else
				{
					planeDataConfigs.push_back(planeDataConfig);
					ops.push_back(&op);
				}
			}

			MPI_Comm comm = MPI_COMM_WORLD;
			MPI_Comm comm_x, comm_t;

			//////braid_Core core;
			///sweet_App* app = (sweet_App *) malloc(sizeof(sweet_App))
			int nt = (int) (shackDict.timecontrol.max_simulation_time / shackDict.timecontrol.current_timestep_size);
                        if (nt * shackDict.timecontrol.current_timestep_size < shackDict.timecontrol.max_simulation_time - 1e-10)
				nt++;
			///sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., shackDict.timecontrol.max_simulation_time, nt, &shackDict, planeDataConfig, &op);
			sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., shackDict.timecontrol.max_simulation_time, nt, &shackDict, planeDataConfigs, ops);


			if ( shackDict.xbraid.xbraid_run_wrapper_tests)
			{
				app.setup();

				BraidUtil braid_util;
				int test = braid_util.TestAll(&app, comm, stdout, 0., shackDict.timecontrol.current_timestep_size, shackDict.timecontrol.current_timestep_size * 2);
				////int test = braid_util.TestBuf(app, comm, stdout, 0.);
				if (test == 0)
					SWEETError("Tests failed!");
				else
					std::cout << "Tests successful!" << std::endl;

			}
			else
			{
				BraidCore core(MPI_COMM_WORLD, &app);
				app.setup(core);
				// Run Simulation
				core.Drive();
			}


			if (shackDict.xbraid.xbraid_spatial_coarsening)
				for (int i = 0; i < shackDict.xbraid.xbraid_max_levels; i++)
				{
					delete planeDataConfigs[i];
					delete ops[i];
					planeDataConfigs[i] = nullptr;
					ops[i] = nullptr;
				}

		}
		else
#endif

}

#endif




int main_mpi(int i_argc, char *i_argv[])
{
	StopwatchBox::getInstance().main.start();

	ProgramPDESWEPlane simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

#if SWEET_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

#if SWEET_GUI
	if (simulation.shackIOData->gui_enabled)
	{
		VisSweet visSweet(simulation);
	}
	else
#endif
	{

#if SWEET_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(*(simulation.shackTimestepControl));

		if (simulation.shackPDESWEPlane->normal_mode_analysis_generation > 0)
		{
			simulation.normal_mode_analysis();
		}
		else
		{
			StopwatchBox::getInstance().main_timestepping.start();

			while (!simulation.should_quit())
			{
				simulation.runTimestep();

				// Instability
				if (simulation.shackPDESWEPlane->instability_checks)
				{
					if (simulation.instability_detected())
						SWEETError("INSTABILITY DETECTED");
				}
			}

			StopwatchBox::getInstance().main_timestepping.stop();
		}
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);
	}

	if (isMPIRoot())
	{
		if (simulation.shackIOData->output_file_name.size() > 0)
			std::cout << "[MULE] reference_filenames: " << simulation.output_filenames << std::endl;

		// End of run output results
		std::cout << "***************************************************" << std::endl;
		std::cout << "Number of time steps: " << simulation.shackTimestepControl->current_timestep_nr << std::endl;
		std::cout << "Time per time step: " << StopwatchBox::getInstance().main_timestepping()/(double)simulation.shackTimestepControl->current_timestep_nr << " sec/ts" << std::endl;
		std::cout << "Last time step size: " << simulation.shackTimestepControl->current_timestep_size << std::endl;

		simulation.computeErrors();
		simulation.printErrors();

		std::cout << "[MULE] simulation_successfully_finished: 1" << std::endl;
	}

	StopwatchBox::getInstance().main.stop();

#if SWEET_MPI
	if (mpi_rank == 0)
#endif
	{
		std::cout << std::endl;
		StopwatchBox::getInstance().output();
	}

	if (isMPIRoot())
	{
		std::cout << "FIN" << std::endl;
	}
	return 0;
}



int main(int i_argc, char *i_argv[])
{
	return main_mpi(i_argc, i_argv);
}
