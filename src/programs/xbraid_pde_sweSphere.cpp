/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/benchmarks/
 *
 */

#include "pde_sweSphere/ProgramXBraidPDESWESphere.hpp"

int main(int i_argc, char *i_argv[])
{

#if SWEET_MPI
	int mpi_rank;
	MPI_Init(&i_argc, &i_argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	ProgramXBraidPDESWESphere simulation(
						i_argc, i_argv
#if SWEET_MPI
						,
						MPI_COMM_WORLD,
						mpi_rank
#endif
			);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(*(simulation.shackTimestepControl));

		simulation.runXBraid();

	}

	///simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

	std::cout << "FIN" << std::endl;

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}
