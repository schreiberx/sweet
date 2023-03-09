/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/time
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/benchmarks
 *
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#include "pde_sweSphere/ProgramPDESWESphere.hpp"


int main(int i_argc, char *i_argv[])
{
	ProgramPDESWESphere simulation(i_argc, i_argv);
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
			simulation.runTimestep();
	}

	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

#if 0
	if (simulation.shackPDESWESphere->compute_errors)
	{
		simulation.printSimulationErrors();
		ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);
	}
#endif

	std::cout << "FIN" << std::endl;
	return 0;
}
