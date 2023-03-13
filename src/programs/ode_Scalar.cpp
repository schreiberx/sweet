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
	ProgramODEScalar simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	{
		simulation.shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(*(simulation.shackTimestepControl));

		while (!simulation.should_quit())
			simulation.runTimestep();
	}

	simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	std::cout << "FIN" << std::endl;
	return 0;
}
