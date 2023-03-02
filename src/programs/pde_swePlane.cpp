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

	// default values for specific input (for general input see SimulationVariables.hpp)
	//simVars.bogus.var[0];  //frequency in x for waves test case
	//simVars.bogus.var[1];  //frequency in y for waves test case

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, user_defined_prog_params))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--initial-freq-x-mul [float]" << std::endl;
		std::cout << "	--initial-freq-y-mul [float]" << std::endl;
		std::cout << "" << std::endl;


#if SWEET_PARAREAL
		simVars.parareal.printProgramArguments();
#endif
		return -1;
	}

	if (simVars.misc.verbosity > 5)
		std::cout << " + Setting up FFT plans..." << std::flush;

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

	if (simVars.misc.verbosity > 5)
		std::cout << " done" << std::endl;

	planeDataConfig->printInformation();

	// Print header
	std::cout << std::endl;
	simVars.outputConfig();
	std::cout << "Computing error: " << simVars.misc.compute_errors << std::endl;
	std::cout << std::endl;

	std::ostringstream buf;
	buf << std::setprecision(14);


#if SWEET_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	#if (SWEET_PARAREAL != 2) && (!SWEET_XBRAID)
	// only start simulation and time stepping for first rank
	if (mpi_rank == 0)
	#endif
#endif
	{
#if SWEET_MPI
		std::cout << "MPI_RANK: " << mpi_rank << std::endl;
#endif

		SimulationBenchmarkTimings::getInstance().main.start();

#if SWEET_PARAREAL
		if (simVars.parareal.enabled)
		{

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			// Set planeDataConfig and planeOperators for each level
			std::vector<PlaneDataConfig*> planeDataConfigs;
			std::vector<PlaneOperators*> ops;

			// fine
			planeDataConfigs.push_back(planeDataConfig);
			ops.push_back(&op);

			// coarse
			if (simVars.parareal.spatial_coarsening)
			{
				///for (int j = 0; j < 2; j++)
				///	assert(simVars.disc.space_res_physical[j] == -1);
				int N_physical[2] = {-1, -1};
				int N_spectral[2];
				double frac;
				if ( simVars.parareal.coarse_timestep_size > 0)
					frac = simVars.timecontrol.current_timestep_size / simVars.parareal.coarse_timestep_size;
				else
					frac = simVars.timecontrol.current_timestep_size / (simVars.timecontrol.max_simulation_time / simVars.parareal.coarse_slices );
				for (int j = 0; j < 2; j++)
					N_spectral[j] = std::max(4, int(simVars.disc.space_res_spectral[j] * frac));
				planeDataConfigs.push_back(new PlaneDataConfig);
				planeDataConfigs.back()->setupAuto(N_physical, N_spectral, simVars.misc.reuse_spectral_transformation_plans);

				ops.push_back(new PlaneOperators(planeDataConfigs.back(), simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs));
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
			Parareal_Controller<SWE_Plane_TimeSteppers, 3> parareal_Controller(	&simVars,
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

			if (simVars.parareal.spatial_coarsening)
			{
				delete planeDataConfigs[1];
				delete ops[1];
			}
		}
		else
#endif


#if SWEET_XBRAID

		if (simVars.xbraid.xbraid_enabled)
		{

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			// Set planeDataConfig and planeOperators for each level
			std::vector<PlaneDataConfig*> planeDataConfigs;
			std::vector<PlaneOperators*> ops;
			for (int i = 0; i < simVars.xbraid.xbraid_max_levels; i++)
			{
				if (simVars.xbraid.xbraid_spatial_coarsening)
				{
					int N_physical[2] = {-1, -1};
					int N_spectral[2];
					for (int j = 0; j < 2; j++)
					{
						// proportional to time step
						if (simVars.xbraid.xbraid_spatial_coarsening == 1)
							N_spectral[j] = std::max(4, int(simVars.disc.space_res_spectral[j] / std::pow(simVars.xbraid.xbraid_cfactor, i)));
						else if (simVars.xbraid.xbraid_spatial_coarsening > 1)
						{
							if (i == 0)
								N_spectral[j] = std::max(4, simVars.disc.space_res_spectral[j]);
							else
								N_spectral[j] = std::max(4, simVars.xbraid.xbraid_spatial_coarsening);
						}
						else
							SWEETError("Invalid parameter xbraid_spatial_coarsening");
					}
					planeDataConfigs.push_back(new PlaneDataConfig);
					planeDataConfigs.back()->setupAuto(N_physical, N_spectral, simVars.misc.reuse_spectral_transformation_plans);

					//PlaneOperators op_level(planeDataConfigs.back(), simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);
					ops.push_back(new PlaneOperators(planeDataConfigs.back(), simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs));

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
			int nt = (int) (simVars.timecontrol.max_simulation_time / simVars.timecontrol.current_timestep_size);
                        if (nt * simVars.timecontrol.current_timestep_size < simVars.timecontrol.max_simulation_time - 1e-10)
				nt++;
			///sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., simVars.timecontrol.max_simulation_time, nt, &simVars, planeDataConfig, &op);
			sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., simVars.timecontrol.max_simulation_time, nt, &simVars, planeDataConfigs, ops);


			if ( simVars.xbraid.xbraid_run_wrapper_tests)
			{
				app.setup();

				BraidUtil braid_util;
				int test = braid_util.TestAll(&app, comm, stdout, 0., simVars.timecontrol.current_timestep_size, simVars.timecontrol.current_timestep_size * 2);
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


			if (simVars.xbraid.xbraid_spatial_coarsening)
				for (int i = 0; i < simVars.xbraid.xbraid_max_levels; i++)
				{
					delete planeDataConfigs[i];
					delete ops[i];
					planeDataConfigs[i] = nullptr;
					ops[i] = nullptr;
				}

		}
		else
#endif


#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			ProgramPDESWEPlane *simulationSWE = new ProgramPDESWEPlane;

			SimulationBenchmarkTimings::getInstance().main_timestepping.start();

			VisSweet<ProgramPDESWEPlane> visSweet(simulationSWE);

			SimulationBenchmarkTimings::getInstance().main_timestepping.stop();

			delete simulationSWE;
		}
		else
#endif
		{
			ProgramPDESWEPlane *simulationSWE = new ProgramPDESWEPlane;
			//Setting initial conditions and workspace - in case there is no GUI

			// also initializes diagnostics
			// already called in constructor
			//simulationSWE->reset();

#if SWEET_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif


			if (simVars.misc.normal_mode_analysis_generation > 0)
			{
				simulationSWE->normal_mode_analysis();
			}
			else
			{
				SimulationBenchmarkTimings::getInstance().main_timestepping.start();

				// Main time loop
				while (true)
				{
					// Stop simulation if requested
					if (simulationSWE->should_quit())
						break;

					// Main call for timestep run
					simulationSWE->run_timestep();

					// Instability
					if (simVars.misc.instability_checks)
					{
						if (simulationSWE->instability_detected())
							SWEETError("INSTABILITY DETECTED");
					}
				}

				SimulationBenchmarkTimings::getInstance().main_timestepping.stop();
			}


			if (simVars.iodata.output_file_name.size() > 0)
				std::cout << "[MULE] reference_filenames: " << simulationSWE->output_filenames << std::endl;

			// End of run output results
			std::cout << "***************************************************" << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << SimulationBenchmarkTimings::getInstance().main_timestepping()/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;

			simulationSWE->compute_errors();


#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				if (simVars.misc.verbosity > 0)
				{
					std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((simVars.diag.total_energy-simulationSWE->diagnostics_energy_start)/simulationSWE->diagnostics_energy_start) << std::endl;
					std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((simVars.diag.total_mass-simulationSWE->diagnostics_mass_start)/simulationSWE->diagnostics_mass_start) << std::endl;
					std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((simVars.diag.total_potential_enstrophy-simulationSWE->diagnostics_potential_enstrophy_start)/simulationSWE->diagnostics_potential_enstrophy_start) << std::endl;

					if (simVars.misc.compute_errors)
					{
						std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark.t0_error_max_abs_h_pert << std::endl;
						std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark.t0_error_max_abs_u << std::endl;
						std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark.t0_error_max_abs_v << std::endl;
					}

					std::cout << "[MULE] error_end_linf_h_pert: " << simulationSWE->benchmark.t0_error_max_abs_h_pert << std::endl;
					std::cout << "[MULE] error_end_linf_u: " << simulationSWE->benchmark.t0_error_max_abs_u << std::endl;
					std::cout << "[MULE] error_end_linf_v: " << simulationSWE->benchmark.t0_error_max_abs_v << std::endl;
					std::cout << std::endl;
				}

				if (simulationSWE->compute_error_to_analytical_solution)
				{
					std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << simulationSWE->benchmark.analytical_error_rms_h << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << simulationSWE->benchmark.analytical_error_rms_u << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << simulationSWE->benchmark.analytical_error_rms_v << std::endl;

					std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << simulationSWE->benchmark.analytical_error_maxabs_h << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << simulationSWE->benchmark.analytical_error_maxabs_u << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << simulationSWE->benchmark.analytical_error_maxabs_v << std::endl;
				}
			}

			std::cout << "[MULE] simulation_successfully_finished: 1" << std::endl;


			delete simulationSWE;
		} // end of gui not enabled

		SimulationBenchmarkTimings::getInstance().main.stop();

#if SWEET_MPI
		if (mpi_rank == 0)
#endif
		{
			std::cout << std::endl;
			SimulationBenchmarkTimings::getInstance().output();
		}
	}
#if SWEET_MPI && (SWEET_PARAREAL != 2) && (!SWEET_XBRAID)
	else
///#if SWEET_MPI
///	#if (SWEET_PARAREAL != 2) && (!SWEET_XBRAID)
///	else	// mpi_rank != 0
///	#else
///	if (mpi_rank != 0)
///	#endif
	{

		if (simVars.disc.timestepping_method.find("rexi") != std::string::npos)
		{
			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			SWE_Plane_TS_l_rexi rexiSWE(simVars, op);

			/*
			 * Setup our little dog REXI
			 */
			rexiSWE.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);

			sweet::PlaneData_Spectral prog_h_pert(planeDataConfig);
			sweet::PlaneData_Spectral prog_u(planeDataConfig);
			sweet::PlaneData_Spectral prog_v(planeDataConfig);

			MPI_Barrier(MPI_COMM_WORLD);

			do
			{
				// REXI time stepping
				rexiSWE.run_timestep(
						prog_h_pert,
						prog_u,
						prog_v,
						simVars.timecontrol.current_timestep_size,
						simVars.timecontrol.current_simulation_time

				);
			}
			while(!rexiSWE.final_timestep);
		}
	}

	if (simVars.disc.timestepping_method.find("rexi") != std::string::npos)
	{
		// synchronize REXI
		if (mpi_rank == 0)
			SWE_Plane_TS_l_rexi::MPI_quitWorkers(planeDataConfig);
	}
#endif

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}

#endif

int main(int i_argc, char *i_argv[])
{
	ProgramPDESWEPlane simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

#if SWEET_GUI
	if (simulation.shackIOData->gui_enabled)
	{
		VisSweet visSweet(simulation);
	}
	else
#endif
	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(*(simulation.shackTimestepControl));

		while (!simulation.should_quit())
			simulation.run_timestep();
	}

	//simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	std::cout << "FIN" << std::endl;
	return 0;
}
