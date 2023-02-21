/*
 * ShackParallelization.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKPARALLELIZATION_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKPARALLELIZATION_HPP_

#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include "../shacks/ShackInterface.hpp"



/**
 * Parallelization
 */
class Parallelization	:
		public sweet::ClassDictionaryInterface
{
public:
	/// number of threads
	int num_threads_space = -1;

	void outputConfig()
	{
		printClass();
	}

	void setup_longOptionsList(
			struct option *long_options,
			int &next_free_program_option
	)
	{
		long_options[next_free_program_option] = {"num-threads-space", required_argument, 0, 256+next_free_program_option};
		next_free_program_option++;
	}

	void outputProgParams()
	{
		printProgramArguments();
	}

	/*
	 * This method is called to parse a particular
	 * long option related to some ID.
	 *
	 * \return: -1 if the option has been processed
	 */
	int setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{
		switch(i_option_index)
		{
		case 0:
			num_threads_space = atoi(i_value);
			return -1;

		}

		return 0;
	}


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Parallelization options:" << std::endl;
		std::cout << "	--num-threads-space [int]			Specify how many threads to use for spatial parallelization (very useful for nested parallel regions)" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--num-threads-space", num_threads_space);

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "PARALLELIZATION:" << std::endl;
		std::cout << " + num_threads: " << num_threads_space << std::endl;
		std::cout << std::endl;
	}
};




#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKPARALLELIZATION_HPP_ */
