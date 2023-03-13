/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <iostream>

#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>

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
	sweet::ShackIOData sweParametersSphere;

	/*
	 * After registering all classes, we can check whether we should output the help information
	 */
	if (pa.argumentWithKeyExists("-h") || pa.argumentWithKeyExists("--help"))
	{
		sweParametersSphere.printProgramArguments();
		return EXIT_FAILURE;
	}

	sweParametersSphere.processProgramArguments(pa);

	std::cout << " + sweParametersSphere->printShack()" << std::endl;
	sweParametersSphere.printShack("    ");

	/*
	 * If you activate this, this should indeed trigger an error
	 */
	{
		bool dummy;
		if (!pa.getArgumentValueByKey("-doesntexist", dummy, true))
		{
			std::cout << "Key not found, but this is on purpose" << std::endl;
			pa.error.print();
		}
		else
		{
			std::cerr << "This should have triggered an error, but it didn't. Stopping here." << std::endl;
			return EXIT_FAILURE;
		}
	}


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

