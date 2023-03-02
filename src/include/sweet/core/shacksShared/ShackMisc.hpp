/*
 * ShackMisc.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKMISC_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKMISC_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>


namespace sweet
{

/**
 * Miscellaneous variables
 */
class ShackMisc	:
		public ShackInterface
{
public:
	/// do instability checks for simulation
	int instability_checks = 1;

	/*
	 * Some flexible variable where one can just add options like
	 * --comma-separated-tags=galewsky_analytical_geostrophic_setup
	 */
	std::string comma_separated_tags = "";


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Misc options:" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--instability-checks", instability_checks);
		i_pa.getArgumentValueByKey("--comma-separated-tags", comma_separated_tags);

		return error.forwardWithPositiveReturn(i_pa.error);
	}


	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "MISC:" << std::endl;
		std::cout << " + instability_checks: " << instability_checks << std::endl;
		std::cout << " + comma_separated_tags: " << comma_separated_tags << std::endl;
		std::cout << std::endl;
	}
};

}


#endif
