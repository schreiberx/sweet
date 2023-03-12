/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SRC_SHACK_SWE_POLVANI_SIMULATION_VARIABLES_HPP_
#define SRC_SHACK_SWE_POLVANI_SIMULATION_VARIABLES_HPP_


#include <getopt.h>
#include <sweet/core/SWEETError.hpp>

#include <sweet/core/shacks/ShackInterface.hpp>



struct ShackPDESWEBenchPolvani	:
	public sweet::ShackInterface
{
	/**
	 * Rossby number
	 */
	double r = 0.01;

	/**
	 * Froude number
	 */
	double f = 0.04;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Polvani benchmark settings (on the plane):" << std::endl;
		std::cout << "	--polvani-rossby [float]	Choose Rossby number, default:0" << std::endl;
		std::cout << "	--polvani-froude [float]	Choose Froude number, default:0" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--polvani-rossby", r);
		i_pa.getArgumentValueByKey("--polvani-froude", f);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "SWEPolvani:" << std::endl;
		std::cout << " + Rossby number: " << r << std::endl;
		std::cout << " + Froude number: " << f << std::endl;
		std::cout << std::endl;
	}
};



#endif
