/*
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SWE_PLANE_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SWE_PLANE_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>



/**
 * simulation coefficients
 */
class ShackPDESWEPlane	:
		public sweet::ShackInterface
{
public:
	/*
	 * Average height for initialization
	 *
	 * We use a default value similar to the Galewski benchmark
	 */
	double h0 = 10000.0;


	/*
	 * For more information on viscosity,
	 *
	 * see 13.3.1 "Generic Form of the Explicit Diffusion Mechanism"
	 * in "Numerical Techniques for Global Atmospheric Models"
	 */

	/*
	 * viscosity-term on velocities with 2nd order diff operator
	 */
	double viscosity = 0.0;

	/*
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
	 * Avoid nonlinear divergence and only solve linear one
	 */
	bool use_only_linear_divergence = false;


	/**
	 * Diffusion applied only on nonlinear divergence
	 */
	int use_nonlinear_only_visc = 0;

	/*
	 * Do a normal mode analysis, see
	 * Hillary Weller, John Thuburn, Collin J. Cotter,
	 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
	 */
	int normal_mode_analysis_generation = 0;

	/*
	 * Compute errors compared to analytical solution
	 */
	bool compute_errors = false;


	/*
	 * Check for instabilities and stop
	 */
	bool instability_checks = false;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "PDESWEPlane parameters:" << std::endl;
		std::cout << i_prefix << "	-u [visc]	viscosity, , default=0" << std::endl;
		std::cout << i_prefix << "	-U [visc]	viscosity order, default=2" << std::endl;
		std::cout << i_prefix << "	-f [float]	f-parameter for f-plane, default=0" << std::endl;
		std::cout << i_prefix << "	-g [float]	gravity" << std::endl;
		std::cout << i_prefix << "	-H [float]	average (initial) height of water" << std::endl;
		std::cout << i_prefix << "	--use-only-linear-divergence [bool]	Use only linear divergence" << std::endl;
		std::cout << i_prefix << "	--use-nonlinear-only-visc [bool]	Use only viscosity on nonlinear part" << std::endl;
		std::cout << i_prefix << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << i_prefix << "	--normal-mode-analysis-generation=[int]			Control generation of normal mode analysis (default: 0)" << std::endl;
		std::cout << i_prefix << "	--instability-checks=[bool]			Check for instabilities (default: 0)" << std::endl;
		std::cout << i_prefix << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("-u", viscosity);
		i_pa.getArgumentValueByKey("-U", viscosity_order);

		i_pa.getArgumentValueByKey("-f", plane_rotating_f0);

		i_pa.getArgumentValueByKey("-g", gravitation);
		i_pa.getArgumentValueByKey("-H", h0);

		i_pa.getArgumentValueByKey("--use-only-linear-divergence", use_only_linear_divergence);
		i_pa.getArgumentValueByKey("--use-nonlinear-only-visc", use_nonlinear_only_visc);

		i_pa.getArgumentValueByKey("--normal-mode-analysis-generation", normal_mode_analysis_generation);
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);

		i_pa.getArgumentValueByKey("--instability-checks", instability_checks);

		ERROR_CHECK_WITH_RETURN_BOOLEAN(i_pa);
		return true;
	}


	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << i_prefix << "SIMULATION SWE PLANE:" << std::endl;
		std::cout << i_prefix << " + h0: " << h0 << std::endl;
		std::cout << i_prefix << " + viscosity: " << viscosity << std::endl;
		std::cout << i_prefix << " + viscosity_order: " << viscosity_order << std::endl;
		std::cout << i_prefix << " + gravitation: " << gravitation << std::endl;
		std::cout << i_prefix << " + plane_rotating_f0: " << plane_rotating_f0 << std::endl;
		std::cout << i_prefix << " + use_only_linear_divergence: " << use_only_linear_divergence << std::endl;
		std::cout << i_prefix << " + use_nonlinear_only_visc: " << use_nonlinear_only_visc << std::endl;
		std::cout << i_prefix << " + compute_errors: " << compute_errors << std::endl;
		std::cout << i_prefix << " + normal_mode_analysis_generation: " << normal_mode_analysis_generation << std::endl;
		std::cout << i_prefix << " + instability_checks: " << instability_checks << std::endl;
		std::cout << std::endl;
	}
};




#endif
