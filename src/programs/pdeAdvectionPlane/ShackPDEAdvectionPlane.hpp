/*
 *  Created on: Feb 23, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SHACK_PDEADVECTIONPLANECOEFFICIENTS_HPP_
#define SRC_PROGRAMS_SHACK_PDEADVECTIONPLANECOEFFICIENTS_HPP_

#include <sweet/shacks/ShackInterface.hpp>


/**
 * Simulation coefficients for the Advection PDE on the plane
 */
class ShackPDEAdvectionPlane	:
		public sweet::ShackInterface
{
public:
	/**
	 * Velocity and additional parameter for advection test cases
	 */
	//double advection_velocity[3] = {0, 0, 0};
#if 0
	bool validateNonzeroAdvection()
	{
		if (advection_velocity[0] == 0 && advection_velocity[1] == 0)
			return error.set("Both advection velocities are 0, use --advection-velocity=...");

		return true;
	}
#endif
	void printProgramArguments(const std::string& i_prefix = "")
	{
		//std::cout << "	--advection-velocity=[float],[float],[float]	advection velocity components (x, y, rotational)" << std::endl;
		//std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		return true;
		/*
		std::string tmp;
		if (i_pa.getArgumentValueByKey("--advection-velocity", tmp))
		{
			StringSplit::split3double(tmp, &advection_velocity[0], &advection_velocity[1], &advection_velocity[2]);
		}

		return error.forwardWithPositiveReturn(i_pa.error);
		*/
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		/*
		std::cout << "PDE ADVECTION:" << std::endl;
		std::cout << " + advection_velocity (x, y, rotation speed): " << advection_velocity[0] << ", " << advection_velocity[1] << ", " << advection_velocity[2] << std::endl;
		std::cout << std::endl;
		*/
	}
};


#endif
