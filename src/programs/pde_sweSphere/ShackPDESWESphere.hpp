/*
 * PDESWEParametersCommon.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_COMMON_PDESWEPARAMETERSCOMMON_HPP_
#define SRC_PROGRAMS_SWE_COMMON_PDESWEPARAMETERSCOMMON_HPP_

#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/ErrorBase.hpp>


class ShackPDESWESphere
{
public:
	virtual sweet::ErrorBase& getError() = 0;

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


	/**
	 * Simulation on f-sphere? (constant f0 term over entire sphere)
	 */
	bool sphere_use_fsphere = false;

	/**
	 * Coriolis effect
	 * 7.2921 x 10^{-5}
	 */
	double sphere_rotating_coriolis_omega = 0.000072921;

	double sphere_fsphere_f0 = 0.00007292*2; //Sphere

	void printShack(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << " + h0: " << h0 << std::endl;
		std::cout << i_prefix << " + viscosity: " << viscosity << std::endl;
		std::cout << i_prefix << " + viscosity_order: " << viscosity_order << std::endl;
		std::cout << i_prefix << " + gravitation: " << gravitation << std::endl;
		std::cout << i_prefix << " + sphere_rotating_coriolis_omega: " << sphere_rotating_coriolis_omega << std::endl;
		std::cout << i_prefix << " + sphere_use_fsphere: " << sphere_use_fsphere << std::endl;
		std::cout << i_prefix << " + sphere_fsphere_f0: " << sphere_fsphere_f0 << std::endl;
	}


	void outputProgArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "	-H [float]	Average (initial) height of water" << std::endl;
		std::cout << i_prefix << "	-u [visc]	Viscosity, , default=0" << std::endl;
		std::cout << i_prefix << "	-U [visc]	Viscosity order, default=2" << std::endl;
		std::cout << i_prefix << "	-g [float]	Gravitation" << std::endl;
	}


	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		if (!i_pa.getArgumentValueBy3Keys("--pde-h0", "-H", "--h0", h0))
		{
			if (getError().forward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy3Keys("--pde-viscosity", "--pde-mu", "--mu", viscosity))
		{
			if (getError().forward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("--pde-viscosity-order", viscosity_order))
		{
			if (getError().forward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy3Keys("--pde-g", "-g", "--gravitation", gravitation))
		{
			if (getError().forward(i_pa.error))
				return false;
		}

		i_pa.getArgumentValueByKey("-F", sphere_use_fsphere);
		if (i_pa.getArgumentValueByKey("-f", sphere_rotating_coriolis_omega))
			sphere_fsphere_f0 = sphere_rotating_coriolis_omega;

		return true;
	}
};



#endif
