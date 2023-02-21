/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */

#include <iostream>

#include "../include/sweet/shacks/ShackDictionary.hpp"
#include "swe_sphere_variables/PDESWESphereParameters.hpp"
#include "swe_sphere_variables/IODataParameters.hpp"

int main(int i_argc, char *i_argv[])
{
	/*
	 * We start by setting up the class to parse the program arguments
	 */
	std::cout << " + ProgramArguments()" << std::endl;
	sweet::ProgramArguments pa;
	if (!pa.setup(i_argc, i_argv))
	{
		std::cout << "Error: " << pa.error.get() << std::endl;
		return 1;
	}

	/*
	 * We instantiate a new class PDESWEParametersSphere
	 */
	PDESWEParametersSphere sweParametersSphere;

	/*
	 * After registering all classes, we can check whether we should output the help information
	 */
	if (pa.argumentWithKeyExists("-h") || pa.argumentWithKeyExists("--help"))
	{
		sweParametersSphere.printProgramArguments();
		return EXIT_FAILURE;
	}

	sweParametersSphere.processProgramArguments(pa);

	std::cout << " + sweParametersSphere->printClass()" << std::endl;
	sweParametersSphere.printClass("    ");

#if 0
	// TODO: Activate this if SWEET migrated entirely to this new interface
	if (!pa.checkAllArgumentsProcessed())
	{
		std::cerr << pa.error.get() << std::endl;
		return EXIT_FAILURE;
	}
#endif

	return 0;
}

