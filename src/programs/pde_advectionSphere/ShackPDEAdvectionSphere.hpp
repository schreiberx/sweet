/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_ADVECTIONSPHERE_SHACKPDEADVECTIONSPHERE_HPP_
#define SRC_PROGRAMS_PDE_ADVECTIONSPHERE_SHACKPDEADVECTIONSPHERE_HPP_


#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>


/**
 * Coefficients for the PDE Advection equation on the sphere
 */
class ShackPDEAdvectionSphere	:
		public sweet::ShackInterface
{
public:
	/*
	 * Compute errors compared to analytical solution
	 */
	bool compute_errors = false;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "PDE advection sphere:" << std::endl;
		std::cout << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);

		return error.forwardWithPositiveReturn(i_pa.error);
	}


	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << i_prefix << "PDE advection sphere:" << std::endl;
		std::cout << i_prefix << " + compute_errors: " << compute_errors << std::endl;
		std::cout << std::endl;
	}
};


#endif
