/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/core/shacks/ShackDictionary.hpp>
#include <iostream>

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


	}

	return 0;
}

