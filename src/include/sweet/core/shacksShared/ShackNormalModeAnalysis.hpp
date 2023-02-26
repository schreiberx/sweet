/*
 *  Created on: Feb 26, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACK_NORMAL_MODE_ANALYSIS_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACK_NORMAL_MODE_ANALYSIS_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>


namespace sweet
{

/**
 * Miscellaneous variables
 */
class ShackNormalModeAnalysis :
		public ShackInterface
{
public:
	/*
	 * Do a normal mode analysis, see
	 * Hillary Weller, John Thuburn, Collin J. Cotter,
	 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
	 */
	int normal_mode_analysis_generation = 0;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Normal mode analysis:" << std::endl;
		std::cout << "	--normal-mode-analysis-generation=[int]			Control generation of normal mode analysis (default: 0)" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--normal-mode-analysis-generation", normal_mode_analysis_generation);

		return error.forwardWithPositiveReturn(i_pa.error);
	}


	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "Normal mode analysis:" << std::endl;
		std::cout << i_prefix << " + normal_mode_analysis_generation: " << normal_mode_analysis_generation << std::endl;
		std::cout << std::endl;
	}
};

}


#endif
