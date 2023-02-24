/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <iostream>

#include <sweet/shacks/ShackDictionary.hpp>
#include "../include/sweet/shacksShared/ShackIOData.hpp"
#include "../include/sweet/shacksShared/ShackParallelization.hpp"


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


	{
		/*
		 * Now we warmup the dictionary which can be stored to store arbitrary
		 * classes.
		 *
		 * This is specialized for processing our program arguments and using it
		 * later on in the SWEET programs.
		 */
		std::cout << " + VariablesClassDictionary()" << std::endl;
		sweet::ShackDictionary varClassDict;


		/*
		 * Register new classes
		 */
		std::cout << "   + registerParameterClass<ShackIOData>()" << std::endl;
		varClassDict.registerFirstTime<ShackIOData>();
		varClassDict.registerFirstTime<ShackParallelization>();

		/*
		 * Now we close the registration
		 *
		 * This will avoid performance bugs!
		 */
		varClassDict.closeRegistration();

		{
			/*
			 * If we now try to register a new class, this should raise an error!
			 */
			bool retval = varClassDict.registerFirstTime<ShackParallelization>();

			if (!retval)
			{
				// Just get the error (deleting it) and continue
				varClassDict.error.get();
			}
			else
			{
				std::cerr << "varClassDict.registerFirstTime<> should have raised an error!" << std::endl;
				return EXIT_FAILURE;
			}
		}

		/*
		 * After registering all classes, we can check whether we should output the help information
		 */
		if (pa.argumentWithKeyExists("-h") || pa.argumentWithKeyExists("--help"))
		{
			varClassDict.printProgramArguments();
			return EXIT_FAILURE;
		}

		/*
		 * Now its time to process all program arguments with all registered classes
		 */
		varClassDict.processProgramArguments(pa);

		/*
		 * Get handler to new class ShackIOData
		 */
		ShackIOData *shackIOData = varClassDict.get<ShackIOData>();
		if (shackIOData == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varClassDict.error.get() << std::endl;
			return EXIT_FAILURE;
		}

		/*
		 * Get handler to new class ShackIOData
		 */
		ShackParallelization *shackParallelization = varClassDict.get<ShackParallelization>();
		if (shackParallelization == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varClassDict.error.get() << std::endl;
			return EXIT_FAILURE;
		}

		/*
		 * Now we close getting the parameter class
		 *
		 * This will avoid performance bugs!
		 */
		varClassDict.closeGet();

		{
			/*
			 * If we now access a class instance, this should raise an error!
			 */

			ShackParallelization *ioDataParameters = varClassDict.get<ShackParallelization>();
			if (ioDataParameters == nullptr)
			{
				// Just get the error (deleting it) and continue
				varClassDict.error.get();
			}
			else
			{
				std::cerr << "This should have raised an error!" << std::endl;
				return EXIT_FAILURE;
			}
		}

		/*
		 * Now its time to output all program arguments of all registered classes
		 */
		std::cout << " + varClassDict.printShack()" << std::endl;
		varClassDict.printShackData("    ");

		/*
		 * And we can also print them individually
		 */
		std::cout << " + sweParametersSphere->printShack()" << std::endl;
		shackIOData->printShack("    ");

		std::cout << " + ioDataParameters->printShack()" << std::endl;
		shackParallelization->printShack("    ");

		/*
		 * Of course, we can also access them directly
		 */
		std::cout << " + direct access:" << std::endl;
		std::cout << "   output_file_name: " << shackIOData->output_file_name << std::endl;
	}

	return 0;
}

