/*
 * ShackSimulationCoefficients.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKSIMULATIONCOEFFICIENTS_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKSIMULATIONCOEFFICIENTS_HPP_

#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>



/**
 * simulation coefficients
 */
class SimulationCoefficients	:
		public sweet::ClassDictionaryInterface
{
public:
	/// average height for initialization
	/// use value similar to Galewski benchmark
	double h0 = 10000.0;


	/// For more information on viscosity,
	/// see 13.3.1 "Generic Form of the Explicit Diffusion Mechanism"
	/// in "Numerical Techniques for Global Atmospheric Models"

	/// viscosity-term on velocities with 2nd order diff operator
	double viscosity = 0.0;

	/// hyper viscosity-term on velocities with 4th order diff operator
	int viscosity_order = 2;


#if SWEET_USE_SPHERE_SPECTRAL_SPACE
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
	double sphere_rotating_coriolis_omega = 0.000072921;

	double sphere_fsphere_f0 = 0.00007292*2; //Sphere
#endif


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


	/**
	 * Velocity and additional parameter for advection test cases
	 */
	double advection_velocity[3] = {0, 0, 0};


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "Simulation parameters:" << std::endl;
		std::cout << "	-X [length]	length of simulation domain in x direction, default=1" << std::endl;
		std::cout << "	-Y [width]	width of simulation domain in y direction, default=1" << std::endl;
		std::cout << "	-u [visc]	viscosity, , default=0" << std::endl;
		std::cout << "	-U [visc]	viscosity order, default=2" << std::endl;
		std::cout << "	-f [float]	f-parameter for f-plane or coriolis omega term, default=0" << std::endl;
		std::cout << "	-F [int]	Simulation on f-sphere, default=0" << std::endl;
		std::cout << "	-g [float]	gravity" << std::endl;
		std::cout << "	-a [float]	earth radius" << std::endl;
		std::cout << "	-H [float]	average (initial) height of water" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("-u", viscosity);
		i_pa.getArgumentValueByKey("-U", viscosity_order);

		i_pa.getArgumentValueByKey("-X", plane_domain_size[0]);
		i_pa.getArgumentValueByKey("-Y", plane_domain_size[1]);

		std::string tmp;
		i_pa.getArgumentValueByKey("--advection-velocity", tmp);
		StringSplit::split3double(tmp, &advection_velocity[0], &advection_velocity[1], &advection_velocity[2]);

		if (i_pa.getArgumentValueByKey("-f", plane_rotating_f0))
		{
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
			sphere_rotating_coriolis_omega = plane_rotating_f0;
			sphere_fsphere_f0 = plane_rotating_f0;
#endif
		}

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		i_pa.getArgumentValueByKey("-F", sphere_use_fsphere);
		i_pa.getArgumentValueByKey("-a", sphere_radius);
#endif

		i_pa.getArgumentValueByKey("-g", gravitation);
		i_pa.getArgumentValueByKey("-H", h0);

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
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
		std::cout << " + advection_velocity (x, y, rotation speed): " << advection_velocity[0] << ", " << advection_velocity[1] << ", " << advection_velocity[2] << std::endl;
		std::cout << " + plane_rotating_f0: " << plane_rotating_f0 << std::endl;

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		std::cout << " + sphere_radius: " << sphere_radius << std::endl;
		std::cout << " + sphere_rotating_coriolis_omega: " << sphere_rotating_coriolis_omega << std::endl;
		std::cout << " + sphere_use_fsphere: " << sphere_use_fsphere << std::endl;
		std::cout << " + sphere_fsphere_f0: " << sphere_fsphere_f0 << std::endl;
#endif

		std::cout << std::endl;
	}
};




#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKSIMULATIONCOEFFICIENTS_HPP_ */
