/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */

#include <iostream>

#include <sweet/shacks/ShackDictionary.hpp>
#include <sweet/shacksShared/ShackShackIOData.hpp>
#include <sweet/shacksShared/ShackShackParallelization.hpp>


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
		sweet::ClassInstanceDictionary varClassDict;


		/*
		 * Register new classes
		 */
		std::cout << "   + registerParameterClass<ShackIOData>()" << std::endl;
		varClassDict.registerClassInstance<ShackIOData>();
		varClassDict.registerClassInstance<ShackParallelization>();

		/*
		 * Now we close the registration
		 *
		 * This will avoid performance bugs!
		 */
		varClassDict.registrationOfClassInstancesFinished();

		{
			/*
			 * If we now try to register a new class, this should raise an error!
			 */
			bool retval = varClassDict.registerClassInstance<ShackParallelization>();

			if (!retval)
			{
				// Just get the error (deleting it) and continue
				varClassDict.error.get();
			}
			else
			{
				std::cerr << "varClassDict.registerClassInstance<> should have raised an error!" << std::endl;
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
		ShackIOData *shackIOData = varClassDict.getClassInstance<ShackIOData>();
		if (shackIOData == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varClassDict.error.get() << std::endl;
			return EXIT_FAILURE;
		}

		/*
		 * Get handler to new class ShackIOData
		 */
		ShackParallelization *shackParallelization = varClassDict.getClassInstance<ShackParallelization>();
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
		varClassDict.getClassInstancesFinished();

		{
			/*
			 * If we now access a class instance, this should raise an error!
			 */

			ShackParallelization *ioDataParameters = varClassDict.getClassInstance<ShackParallelization>();
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
		std::cout << " + varClassDict.printClass()" << std::endl;
		varClassDict.printClass("    ");

		/*
		 * And we can also print them individually
		 */
		std::cout << " + sweParametersSphere->printClass()" << std::endl;
		shackIOData->printClass("    ");

		std::cout << " + ioDataParameters->printClass()" << std::endl;
		shackParallelization->printClass("    ");

		/*
		 * Of course, we can also access them directly
		 */
		std::cout << " + direct access:" << std::endl;
		std::cout << "   output_file_name: " << shackIOData->output_file_name << std::endl;
	}

	return 0;
}

