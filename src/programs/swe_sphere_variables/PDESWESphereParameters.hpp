/*
 * VariablesSWESphere.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_VARIABLES_PDESWESPHEREPARAMETERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_VARIABLES_PDESWESPHEREPARAMETERS_HPP_

#include <sweet/ProgramArguments.hpp>
#include <sweet/ErrorBase.hpp>
#include "../../include/sweet/classDict/ClassDictionaryInterface.hpp"
#include "../swe_common/PDESWEParametersCommon.hpp"


/**
 * SWE PDE on the sphere parameters
 */
class PDESWEParametersSphere	:
	public PDESWEParametersCommon,
	public sweet::ClassDictionaryInterface
{
public:
	sweet::ErrorBase& getError()
	{
		return error;
	}

	/**
	 * Earth radius for simulations on the sphere
	 */
	double sphere_radius = 6.37122e6;

	/**
	 * Simulation on f-sphere? (constant f0 term over entire sphere)
	 */
	bool sphere_use_fsphere = false;

	/**
	 * Coriolis effect
	 * 7.2921 x 10^{-5}
	 */
	double sphere_coriolis_omega = 0.000072921;

	/**
	 * Rotational speed if f-sphere is used (cylindrical)
	 */
	double sphere_fsphere_f0 = 0.00007292*2;


	void printClass(const std::string& i_prefix = "")
	{
		std::cout << std::endl;
		std::cout << i_prefix << "SWE PDE on Sphere coefficients:" << std::endl;
		PDESWEParametersCommon::printClass(i_prefix);
		std::cout << i_prefix << " + radius: " << sphere_radius << std::endl;
		std::cout << i_prefix << " + rotating_coriolis_omega: " << sphere_coriolis_omega << std::endl;
		std::cout << i_prefix << " + use_fsphere: " << sphere_use_fsphere << std::endl;
		std::cout << i_prefix << " + fsphere_f0: " << sphere_fsphere_f0 << std::endl;

		std::cout << std::endl;
	}


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "Simulation parameters:" << std::endl;
		PDESWEParametersCommon::outputProgArguments(i_prefix);
		std::cout << i_prefix << "	-a [float]	earth radius" << std::endl;
		std::cout << i_prefix << "	-f [float]	f-parameter for f-plane or coriolis omega term, default=0" << std::endl;
		std::cout << i_prefix << "	-F [bool]	Simulation on f-sphere, default=0" << std::endl;
		std::cout << i_prefix << "" << std::endl;

	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		if (!i_pa.getArgumentValueBy2Keys("--pde-sphere-radius", "-a", sphere_radius))
		{
			if (error.forwardFrom(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy2Keys("--pde-use-fsphere", "-F", sphere_use_fsphere))
		{
			if (error.forwardFrom(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("--pde-sphere-f0", sphere_fsphere_f0))
		{
			if (error.forwardFrom(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("--pde-coriolis-omega", sphere_coriolis_omega))
		{
			if (error.forwardFrom(i_pa.error))
				return false;
		}

		if (!PDESWEParametersCommon::processProgramArguments(i_pa))
			return false;

		return true;
	}
};



#endif /* SRC_PROGRAMS_SWE_SPHERE_VARIABLES_PDESWESPHEREPARAMETERS_HPP_ */
