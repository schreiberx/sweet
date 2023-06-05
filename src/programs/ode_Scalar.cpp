/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ode_Scalar/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ode_Scalar/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ode_Scalar/benchmarks/
 *
 */

#include "ode_Scalar/ProgramODEScalar.hpp"

int main(int i_argc, char *i_argv[])
{

#if SWEET_MPI
	MPI_Init(&i_argc, &i_argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	ProgramODEScalar simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimestepControl));

		simulation.timestepHandleOutput();

		while (!simulation.should_quit())
		{
			simulation.runTimestep();
			simulation.timestepHandleOutput();
		}
	}

	simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	std::cout << "FIN" << std::endl;

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}
