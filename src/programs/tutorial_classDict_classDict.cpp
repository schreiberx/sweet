/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include <iostream>

#include <sweet/classDict/ClassInstanceDictionary.hpp>
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


	{
		/*
		 * Now we warmup the dictionary which can be stored to store arbitrary
		 * classes.
		 *
		 * This is specialized for processing our program arguments and using it
		 * later on in the SWEET programs.
		 */
		std::cout << " + VariablesClassDictionary()" << std::endl;
		ClassInstanceDictionary varClassDict;


		/*
		 * Register new classes
		 */
		std::cout << "   + registerParameterClass<PDESWEParametersSphere>()" << std::endl;
		varClassDict.registerClassInstance<PDESWEParametersSphere>();
		varClassDict.registerClassInstance<IODataParameters>();

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
			bool retval = varClassDict.registerClassInstance<IODataParameters>();

			if (!retval)
			{
				// Just get the error (deleting it) and continue
				varClassDict.error.errorGet();
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
			varClassDict.outputProgramArguments();
			return EXIT_FAILURE;
		}

		/*
		 * Now its time to process all program arguments with all registered classes
		 */
		varClassDict.processProgramArguments(pa);

		/*
		 * Get handler to new class PDESWEParametersSphere
		 */
		PDESWEParametersSphere *sweParametersSphere = varClassDict.getClassInstance<PDESWEParametersSphere>();
		if (sweParametersSphere == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varClassDict.error.errorGet() << std::endl;
			return EXIT_FAILURE;
		}

		/*
		 * Get handler to new class PDESWEParametersSphere
		 */
		IODataParameters *ioDataParameters = varClassDict.getClassInstance<IODataParameters>();
		if (ioDataParameters == nullptr)
		{
			std::cerr << "Not a SWEET error: " << varClassDict.error.errorGet() << std::endl;
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

			IODataParameters *ioDataParameters = varClassDict.getClassInstance<IODataParameters>();
			if (ioDataParameters == nullptr)
			{
				// Just get the error (deleting it) and continue
				varClassDict.error.errorGet();
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
		std::cout << " + varClassDict.outputVariables()" << std::endl;
		varClassDict.outputVariables("    ");

		/*
		 * And we can also print them individually
		 */
		std::cout << " + sweParametersSphere->outputVariables()" << std::endl;
		sweParametersSphere->outputVariables("    ");

		std::cout << " + ioDataParameters->outputVariables()" << std::endl;
		ioDataParameters->outputVariables("    ");

		/*
		 * Of course, we can also access them directly
		 */
		std::cout << " + direct access:" << std::endl;
		std::cout << "   h0: " << sweParametersSphere->h0 << std::endl;
	}

	return 0;
}

