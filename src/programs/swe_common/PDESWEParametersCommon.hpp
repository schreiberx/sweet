/*
 * PDESWEParametersCommon.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_COMMON_PDESWEPARAMETERSCOMMON_HPP_
#define SRC_PROGRAMS_SWE_COMMON_PDESWEPARAMETERSCOMMON_HPP_

#include <sweet/ProgramArguments.hpp>
#include <sweet/ErrorBase.hpp>


class PDESWEParametersCommon
{
public:
	sweet::ErrorBase error;

	/**
	 * Average height
	 */
	double h0 = 10000.0;

	/**
	 * For more information on viscosity,
	 * see 13.3.1 "Generic Form of the Explicit Diffusion Mechanism"
	 * in "Numerical Techniques for Global Atmospheric Models"
	 *
	 * viscosity-term on velocities with 2nd order diff operator
	 */

	double viscosity = 0.0;

	/**
	 * Order of viscosity
	 */
	int viscosity_order = 2;

	/**
	 * Gravitational constant
	 */
	double gravitation = 9.80616;


	void outputVariables(std::string i_prefix = "")
	{
		std::cout << i_prefix << " + h0: " << h0 << std::endl;
		std::cout << i_prefix << " + viscosity: " << viscosity << std::endl;
		std::cout << i_prefix << " + viscosity_order: " << viscosity_order << std::endl;
		std::cout << i_prefix << " + gravitation: " << gravitation << std::endl;
	}


	void outputProgArguments(std::string i_prefix = "")
	{
		std::cout << i_prefix << "	-H [float]	Average (initial) height of water" << std::endl;
		std::cout << i_prefix << "	-u [visc]	Viscosity, , default=0" << std::endl;
		std::cout << i_prefix << "	-U [visc]	Viscosity order, default=2" << std::endl;
		std::cout << i_prefix << "	-g [float]	Gravitation" << std::endl;
	}


	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		if (!i_pa.getArgumentValueBy3Keys("pde-h0", "H", "h0", h0))
		{
			if (error.errorForwardFrom(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy3Keys("pde-viscosity", "pde-mu", "mu", viscosity))
		{
			if (error.errorForwardFrom(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("pde-viscosity-order", viscosity_order))
		{
			if (error.errorForwardFrom(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy3Keys("pde-g", "g", "gravitation", gravitation))
		{
			if (error.errorForwardFrom(i_pa.error))
				return false;
		}

		return true;
	}
};



#endif /* SRC_PROGRAMS_SWE_COMMON_PDESWEPARAMETERSCOMMON_HPP_ */
