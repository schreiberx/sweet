/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ode_Scalar/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ode_Scalar/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/ode_Scalar/benchmarks/
 *
 */

#if ! SWEET_MPI
#error "MPI should be enabled"
#endif

#include "ode_Scalar/ProgramODEScalar.hpp"

int main(int i_argc, char *i_argv[])
{


	MPI_Init(&i_argc, &i_argv);

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


	if (shackDict.xbraid.xbraid_enabled)
	{

		MPI_Comm comm = MPI_COMM_WORLD;
		MPI_Comm comm_x, comm_t;

		int nt = (int) (shackDict.timecontrol.max_simulation_time / shackDict.timecontrol.current_timestepSize);
                if (nt * shackDict.timecontrol.current_timestepSize < shackDict.timecontrol.max_simulation_time - 1e-10)
			nt++;
		sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., shackDict.timecontrol.max_simulation_time, nt, &shackDict);

		if ( shackDict.xbraid.xbraid_run_wrapper_tests)
		{

			app.setup();

			BraidUtil braid_util;
			int test = braid_util.TestAll(&app, comm, stdout, 0., shackDict.timecontrol.current_timestepSize, shackDict.timecontrol.current_timestepSize * 2);
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

	}



	ProgramODEScalar simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*(simulation.shackTimestepControl));

		while (!simulation.should_quit())
			simulation.runTimestep();
	}

	simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

	std::cout << "FIN" << std::endl;
	return 0;
}
