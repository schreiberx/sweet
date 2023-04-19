/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane_MoriZwanzig/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane_MoriZwanzig/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane_MoriZwanzig/benchmarks/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_swePlane/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */



#include "pde_swePlane_MoriZwanzig/ProgramPDESWEPlaneMoriZwanzig.hpp"



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

	ProgramPDESWEPlaneMoriZwanzig simulation(i_argc, i_argv);
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

		////if (simulation.shackPDESWEPlane->normal_mode_analysis_generation > 0)
		////{
		////	simulation.normal_mode_analysis();
		////}
		////else
		////{
			StopwatchBox::getInstance().main_timestepping.start();

			simulation.timestep_do_output();
			std::cout << "CHECKING STABILITY" << std::endl;
			simulation.instability_detected();
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
		////}
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
