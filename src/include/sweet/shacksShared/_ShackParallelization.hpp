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
#include <sweet/shacks/ShackInterface.hpp>



/**
 * Parallelization
 */
class ShackParallelization	:
		public sweet::ShackInterface
{
public:
	/// number of threads
	int num_threads_space = -1;

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Parallelization options:" << std::endl;
		std::cout << "	--num-threads-space [int]			Specify how many threads to use for spatial parallelization (very useful for nested parallel regions)" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--num-threads-space", num_threads_space);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printShack(
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
