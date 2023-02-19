/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include <sweet/ProgramArguments.hpp>
#include <sweet/ErrorBase.hpp>
#include <iostream>

class CoefficientsPDESWECommon
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
	 * order of viscosity
	 */
	int viscosity_order = 2;

	/**
	 * Gravitational constant
	 */
	double gravitation = 9.80616;



	void outputConfig()
	{
		std::cout << " + h0: " << h0 << std::endl;
		std::cout << " + viscosity: " << viscosity << std::endl;
		std::cout << " + viscosity_order: " << viscosity_order << std::endl;
		std::cout << " + gravitation: " << gravitation << std::endl;
	}

	void outputProgParams()
	{
		std::cout << "Simulation parameters:" << std::endl;
		CoefficientsPDESWECommon::outputProgParams();
		std::cout << "	-H [float]	Average (initial) height of water" << std::endl;
		std::cout << "	-u [visc]	Viscosity, , default=0" << std::endl;
		std::cout << "	-U [visc]	Viscosity order, default=2" << std::endl;
		std::cout << "	-g [float]	Gravitation" << std::endl;
		std::cout << "" << std::endl;
	}


	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		if (!i_pa.getArgumentValueBy3Keys("pde-h0", "H", "h0", h0))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy3Keys("pde-viscosity", "pde-mu", "mu", viscosity))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("pde-viscosity-order", viscosity_order))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy3Keys("pde-g", "g", "gravitation", gravitation))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		return true;
	}
};


/**
 * SWE PDE on the sphere parameters
 */
class CoefficientsPDESWESphere	:
	public CoefficientsPDESWECommon
{
public:
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


	void outputConfig()
	{
		std::cout << std::endl;
		std::cout << "SWE PDE on Sphere coefficients:" << std::endl;
		CoefficientsPDESWECommon::outputConfig();
		std::cout << " + radius: " << sphere_radius << std::endl;
		std::cout << " + rotating_coriolis_omega: " << sphere_coriolis_omega << std::endl;
		std::cout << " + use_fsphere: " << sphere_use_fsphere << std::endl;
		std::cout << " + fsphere_f0: " << sphere_fsphere_f0 << std::endl;

		std::cout << std::endl;
	}


	void outputProgramArguments()
	{
		std::cout << "Simulation parameters:" << std::endl;
		CoefficientsPDESWECommon::outputProgParams();
		std::cout << "	-a [float]	earth radius" << std::endl;
		std::cout << "	-f [float]	f-parameter for f-plane or coriolis omega term, default=0" << std::endl;
		std::cout << "	-F [bool]	Simulation on f-sphere, default=0" << std::endl;
		std::cout << "" << std::endl;

	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		if (!i_pa.getArgumentValueBy2Keys("pde-sphere-radius", "a", sphere_radius))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueBy2Keys("pde-use-fsphere", "F", sphere_use_fsphere))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("pde-sphere-f0", sphere_fsphere_f0))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("pde-coriolis-omega", sphere_coriolis_omega))
		{
			if (error.errorForward(i_pa.error))
				return false;
		}

		if (!CoefficientsPDESWECommon::processProgramArguments(i_pa))
			return false;

		return true;
	}
};


int main(int i_argc, char *i_argv[])
{

	sweet::ProgramArguments pa;
	if (!pa.setup(i_argc, i_argv))
	{
		std::cout << "Error: " << pa.error.errorGet() << std::endl;
		return 1;
	}

	CoefficientsPDESWESphere c;

	if (pa.argumentWithKeyExists("h") || pa.argumentWithKeyExists("help"))
	{
		c.outputProgramArguments();
		return EXIT_FAILURE;
	}

	if (!c.processProgramArguments(pa))
	{
		std::cerr << "Error: " << c.error.errorGet() << std::endl;
		return EXIT_FAILURE;
	}
	c.outputConfig();

	std::cout << pa << std::endl;

	if (!pa.checkAllArgumentsProcessed())
	{
		std::cerr << "Error: " << pa.error.errorGet() << std::endl;
		return EXIT_FAILURE;
	}

	return 0;
}

