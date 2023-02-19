/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include <iostream>
#include <sweet/class_dict/ClassInstanceDictionary.hpp>
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
		std::cout << "Error: " << pa.error.errorGet() << std::endl;
		return 1;
	}

	/*
	 * We instantiate a new class PDESWEParametersSphere
	 */
	PDESWEParametersSphere sweParametersSphere;

	sweParametersSphere.processProgramArguments(pa);

	std::cout << " + sweParametersSphere->outputVariables()" << std::endl;
	sweParametersSphere.outputVariables("    ");

#if 0
	// TODO: Activate this if SWEET migrated entirely to this new interface
	if (!pa.checkAllArgumentsProcessed())
	{
		std::cerr << pa.error.errorGet() << std::endl;
		return EXIT_FAILURE;
	}
#endif

	return 0;
}

