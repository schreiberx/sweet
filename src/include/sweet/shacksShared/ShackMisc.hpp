/*
 * ShackMisc.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKMISC_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKMISC_HPP_

#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include "../shacks/ShackInterface.hpp"



/**
 * Miscellaneous variables
 */
class Misc	:
		public sweet::ClassDictionaryInterface
{
public:
	/// set verbosity of simulation
	int verbosity = 0;

	/// compute errors
	int compute_errors = 0;

	/// do instability checks for simulation
	int instability_checks = 1;

	/// activate GUI mode?
	bool gui_enabled = (SWEET_GUI == 0 ? false : true);


	/// id for visualization
	int vis_id = 0;

	/// Diffusion applied only on nonlinear divergence
	int use_nonlinear_only_visc = 0;

	/// Load / Save plans for SHTNS (useful for reproducibility)
	TransformationPlans::TRANSFORMATION_PLAN_CACHE reuse_spectral_transformation_plans = TransformationPlans::QUICK;

	/*
	 * Do a normal mode analysis, see
	 * Hillary Weller, John Thuburn, Collin J. Cotter,
	 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
	 */
	int normal_mode_analysis_generation = 0;

	/*
	 * Some flexible variable where one can just add options like
	 * --comma-separated-tags=galewsky_analytical_geostrophic_setup
	 */
	std::string comma_separated_tags = "";

	void outputConfig()
	{
		printClass();
	}



	void outputProgParams()
	{
		printProgramArguments();
	}

	void setup_longOptionsList(
			struct option *long_options,
			int &next_free_program_option
	)
	{
		// MISC
		long_options[next_free_program_option] = {"compute-errors", required_argument, 0, 256+next_free_program_option};
		next_free_program_option++;

		long_options[next_free_program_option] = {"instability-checks", required_argument, 0, 256+next_free_program_option};
		next_free_program_option++;

		long_options[next_free_program_option] = {"use-nonlinear-only-visc", required_argument, 0, 256+next_free_program_option};
		next_free_program_option++;

		long_options[next_free_program_option] = {"reuse-plans", required_argument, 0, 256+next_free_program_option};
		next_free_program_option++;

		long_options[next_free_program_option] = {"normal-mode-analysis-generation", required_argument, 0, 256+next_free_program_option};
		next_free_program_option++;

		long_options[next_free_program_option] = {"comma-separated-tags", required_argument, 0, 256+next_free_program_option};
		next_free_program_option++;
	}



	/*
	 * This method is called to parse a particular
	 * long option related to some ID.
	 *
	 * \return: -1 if the option has been processed
	 */
	int setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{

		switch(i_option_index)
		{
		case 0:
			compute_errors = atoi(i_value);
			return -1;

		case 1:
			instability_checks = atoi(i_value);
			return -1;

		case 2:
			use_nonlinear_only_visc = atoi(i_value);
			return -1;

		case 3:
			reuse_spectral_transformation_plans = TransformationPlans::getEnumFromString(i_value);
			return -1;

		case 4:
			normal_mode_analysis_generation = atoi(i_value);
			return -1;

		case 5:
			comma_separated_tags = i_value;
			return -1;
		}

		return 6;
	}

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Misc options:" << std::endl;
		std::cout << "	-v [int]			verbosity level" << std::endl;
//		std::cout << "	-V [double]			period of outputConfig" << std::endl;
		std::cout << "	-G [0/1]			graphical user interface" << std::endl;
		std::cout << "	-O [string]			string prefix for filename of output of simulation data (default output_%s_t%020.8f.csv)" << std::endl;
		std::cout << "	-d [int]			accuracy of floating point output" << std::endl;
		std::cout << "	-i [file0][;file1][;file3]...	string with filenames for initial conditions" << std::endl;
		std::cout << "					specify BINARY; as first file name to read files as binary raw data" << std::endl;
		std::cout << "	--compute-errors [int]          Compute errors when possible [1], default=0	" << std::endl;
		std::cout << "	--use-local-visc [0/1]	Viscosity will be applied only on nonlinear divergence, default:0" << std::endl;
		std::cout << "	--reuse-plans [0/1]	Save plans for fftw transformations and SH transformations" << std::endl;
		std::cout << "					-1: use only estimated plans (no wisdom)" << std::endl;
		std::cout << "					0: compute optimized plans (no wisdom)" << std::endl;
		std::cout << "					1: compute optimized plans, use wisdom if available and store wisdom" << std::endl;
		std::cout << "					2: use wisdom if available if not, trigger error if wisdom doesn't exist (not yet working for SHTNS)" << std::endl;
		std::cout << "					default: -1 (quick mode)" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);
		i_pa.getArgumentValueByKey("--instability-checks", instability_checks);
		i_pa.getArgumentValueByKey("--use-nonlinear-only-visc", use_nonlinear_only_visc);
		i_pa.getArgumentValueByKey("--reuse-plans", *(int*)&reuse_spectral_transformation_plans);
		i_pa.getArgumentValueByKey("--normal-mode-analysis-generation", normal_mode_analysis_generation);
		i_pa.getArgumentValueByKey("--comma-separated-tags", comma_separated_tags);
		i_pa.getArgumentValueByKey("-G", gui_enabled);
		i_pa.getArgumentValueByKey("-v", verbosity);

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "MISC:" << std::endl;
		std::cout << " + verbosity: " << verbosity << std::endl;
		std::cout << " + compute_errors " << compute_errors << std::endl;
		std::cout << " + instability_checks: " << instability_checks << std::endl;
		std::cout << " + gui_enabled: " << gui_enabled << std::endl;
		std::cout << " + vis_id: " << vis_id << std::endl;
		std::cout << " + use_nonlinear_only_visc: " << use_nonlinear_only_visc << std::endl;
		std::cout << " + reuse_spectral_transformation_plans: " << TransformationPlans::getStringFromEnum(reuse_spectral_transformation_plans) << std::endl;
		std::cout << std::endl;
		std::cout << " + normal_mode_analysis_generation: " << normal_mode_analysis_generation << std::endl;
		std::cout << " + comma_separated_tags: " << comma_separated_tags << std::endl;
		std::cout << std::endl;
	}
};




#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKMISC_HPP_ */
