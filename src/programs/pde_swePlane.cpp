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


int main_mpi(int i_argc, char *i_argv[])
{
	StopwatchBox::getInstance().main.start();

	ProgramPDESWEPlane simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

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
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimestepControl));

		if (simulation.shackPDESWEPlane->normal_mode_analysis_generation > 0)
		{
			simulation.normal_mode_analysis();
		}
		else
		{
			simulation.timestep_do_output();

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
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);
	}

	if (isMPIRoot())
	{
		if (simulation.shackIOData->output_file_name.size() > 0)
			std::cout << "[MULE] reference_filenames: " << simulation.output_filenames << std::endl;

		// End of run output results
		std::cout << "***************************************************" << std::endl;
		std::cout << "Number of time steps: " << simulation.shackTimestepControl->current_timestep_nr << std::endl;
		std::cout << "Time per time step: " << StopwatchBox::getInstance().main_timestepping()/(double)simulation.shackTimestepControl->current_timestep_nr << " sec/ts" << std::endl;
		std::cout << "Last time step size: " << simulation.shackTimestepControl->current_timestepSize << std::endl;

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

	int retval = main_mpi(i_argc, i_argv);

#if SWEET_MPI
	MPI_Finalize();
#endif

	return retval;
}
