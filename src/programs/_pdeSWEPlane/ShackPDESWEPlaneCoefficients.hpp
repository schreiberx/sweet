/*
 *  Created on: Feb 22, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SHACK_SWE_PLANE_COEFFICIENTS_HPP_
#define SRC_INCLUDE_SHACK_SWE_PLANE_COEFFICIENTS_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>

/**
 * Simulation coefficients for the SWE on the plane
 */
class ShackPDESWEPlaneCoefficients	:
		public sweet::ShackInterface
{
public:
	/**
	 * average height for initialization
	 *
	 * Default value similar to Galewski benchmark
	 */
	double h0 = 10000.0;


	/**
	 * For more information on viscosity,
	 * see 13.3.1 "Generic Form of the Explicit Diffusion Mechanism"
	 * in "Numerical Techniques for Global Atmospheric Models"
	 */

	/**
	 * viscosity-term on velocities with 2nd order diff operator
	 */
	double viscosity = 0.0;

	/**
	 * hyper viscosity-term on velocities with 4th order diff operator
	 */
	int viscosity_order = 2;


	/**
	 * Plane with f-Coriolis rotation
	 */
	double plane_rotating_f0 = 1.0; //Plane


	/**
	 * Gravitational constant
	 */
	double gravitation = 9.80616;


	/**
	 * domain size if running simulation on the plane
	 */
	double plane_domain_size[2] = {1.0, 1.0};


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "Simulation parameters:" << std::endl;
		std::cout << "	-X [length]	length of simulation domain in x direction, default=1" << std::endl;
		std::cout << "	-Y [width]	width of simulation domain in y direction, default=1" << std::endl;
		std::cout << "	-u [visc]	viscosity, default=0" << std::endl;
		std::cout << "	-U [visc]	viscosity order, default=2" << std::endl;
		std::cout << "	-f [float]	f-parameter for f-plane or coriolis omega term, default=0" << std::endl;
		std::cout << "	-g [float]	gravity" << std::endl;
		std::cout << "	-H [float]	average (initial) height of water" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("-u", viscosity);
		i_pa.getArgumentValueByKey("-U", viscosity_order);

		i_pa.getArgumentValueByKey("-X", plane_domain_size[0]);
		i_pa.getArgumentValueByKey("-Y", plane_domain_size[1]);

		i_pa.getArgumentValueByKey("-g", gravitation);
		i_pa.getArgumentValueByKey("-H", h0);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "SIMULATION COEFFICIENTS:" << std::endl;
		std::cout << " + h0: " << h0 << std::endl;
		std::cout << " + viscosity: " << viscosity << std::endl;
		std::cout << " + viscosity_order: " << viscosity_order << std::endl;
		std::cout << " + gravitation: " << gravitation << std::endl;
		std::cout << " + domain_size (2D): " << plane_domain_size[0] << " x " << plane_domain_size[1] << std::endl;
		std::cout << " + plane_rotating_f0: " << plane_rotating_f0 << std::endl;

		std::cout << std::endl;
	}
};


#endif
