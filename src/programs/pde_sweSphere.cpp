/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/time
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/benchmarks
 *
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 * MULE_SCONS_OPTIONS: --fortran-source=enable
 * MULE_SCONS_OPTIONS: --lapack=enable
 */

#include <sweet/core/defaultPrecompilerValues.hpp>

#if SWEET_MPI
#	include <mpi.h>
#endif

#include "pde_sweSphere/ProgramPDESWESphere.hpp"


int main_mpi(int i_argc, char *i_argv[])
{
#if SWEET_MPI
	MPI_Init(&i_argc, &i_argv);

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	ProgramPDESWESphere simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);


#if SWEET_GUI
	if (simulation.shackIOData->gui_enabled)
	{
		VisSweet visSweet(simulation);
	}
	else
#endif
	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(*(simulation.shackTimestepControl));

		if (simulation.shackPDESWESphere->normal_mode_analysis_generation > 0)
		{
			simulation.normalmode_analysis();
		}
		else
		{
			simulation.timestepHandleOutput();

			while (!simulation.should_quit())
			{
				simulation.runTimestep();

				if (simulation.shackPDESWESphere->instability_checks)
				{
#if SWEET_MPI
					if (mpi_rank == 0)
#endif
					{
						if (simulation.detect_instability())
						{
							std::cerr << "INSTABILITY DETECTED" << std::endl;
							exit(1);
							break;
						}
					}
				}

#if SWEET_MPI
				if (mpi_rank == 0)
#endif
				{
					simulation.timestepHandleOutput();
				}
			}
		}

		if (simulation.fileOutput.output_reference_filenames.size() > 0)
			std::cout << "[MULE] reference_filenames: " << simulation.fileOutput.output_reference_filenames << std::endl;
	}

	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

	std::cout << "FIN" << std::endl;

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}


int main(int i_argc, char *i_argv[])
{
	return main_mpi(i_argc, i_argv);
}

