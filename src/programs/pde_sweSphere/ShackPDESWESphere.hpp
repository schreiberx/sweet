/*
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_COMMON_PDESWEPARAMETERSCOMMON_HPP_
#define SRC_PROGRAMS_SWE_COMMON_PDESWEPARAMETERSCOMMON_HPP_

#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/ErrorBase.hpp>


class ShackPDESWESphere	:
		public sweet::ShackInterface
{
public:
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
	 * Compute diagnostics
	 */
	bool compute_diagnostics = false;


	/*
	 * Check for instabilities and stop
	 */
	bool instability_checks = false;


	void printShack(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "PDESWESphere parameters:" << std::endl;
		std::cout << i_prefix << " + h0: " << h0 << std::endl;
		std::cout << i_prefix << " + viscosity: " << viscosity << std::endl;
		std::cout << i_prefix << " + viscosity_order: " << viscosity_order << std::endl;
		std::cout << i_prefix << " + gravitation: " << gravitation << std::endl;
		std::cout << i_prefix << " + sphere_rotating_coriolis_omega: " << sphere_rotating_coriolis_omega << std::endl;
		std::cout << i_prefix << " + sphere_use_fsphere: " << sphere_use_fsphere << std::endl;
		std::cout << i_prefix << " + sphere_fsphere_f0: " << sphere_fsphere_f0 << std::endl;
		std::cout << i_prefix << " + compute_errors: " << compute_errors << std::endl;
		std::cout << i_prefix << " + compute_diagnostics: " << compute_diagnostics << std::endl;
		std::cout << i_prefix << " + instability_checks: " << instability_checks << std::endl;
	}


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "	-H [float]	Average (initial) height of water" << std::endl;
		std::cout << i_prefix << "	-u [visc]	Viscosity, , default=0" << std::endl;
		std::cout << i_prefix << "	-U [visc]	Viscosity order, default=2" << std::endl;
		std::cout << i_prefix << "	-g [float]	Gravitation" << std::endl;
		std::cout << i_prefix << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << i_prefix << "	--compute-diagnostics [bool]	Compute diagnostics" << std::endl;
		std::cout << i_prefix << "	--normal-mode-analysis-generation=[int]			Control generation of normal mode analysis (default: 0)" << std::endl;
		std::cout << i_prefix << "	--instability-checks=[bool]			Check for instabilities (default: 0)" << std::endl;
	}


	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueBy3Keys("--pde-h0", "-H", "--h0", h0);
		i_pa.getArgumentValueBy2Keys("--pde-viscosity", "-u", viscosity);
		i_pa.getArgumentValueByKey("--pde-viscosity-order", viscosity_order);
		i_pa.getArgumentValueBy3Keys("--pde-g", "-g", "--pde-gravitation", gravitation);
		i_pa.getArgumentValueByKey("-F", sphere_use_fsphere);

		if (i_pa.getArgumentValueByKey("-f", sphere_rotating_coriolis_omega))
			sphere_fsphere_f0 = sphere_rotating_coriolis_omega;

		i_pa.getArgumentValueByKey("--normal-mode-analysis-generation", normal_mode_analysis_generation);
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);
		i_pa.getArgumentValueByKey("--compute-diagnostics", compute_diagnostics);
		i_pa.getArgumentValueByKey("--instability-checks", instability_checks);
		
		ERROR_CHECK_WITH_RETURN_BOOLEAN(i_pa);
		return true;
	}
};



#endif
