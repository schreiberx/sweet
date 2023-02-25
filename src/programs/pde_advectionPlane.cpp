/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionPlane/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionPlane/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionPlane/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */


#include "pde_advectionPlane/ProgramPDEAdvectionPlane.hpp"


int main(int i_argc, char *i_argv[])
{
	ProgramPDEAdvectionPlane simulation(i_argc, i_argv);
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

	simulation.printSimulationErrors();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	std::cout << "FIN" << std::endl;
	return 0;
}
