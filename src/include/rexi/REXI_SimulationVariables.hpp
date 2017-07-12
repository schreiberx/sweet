/*
 * REXI_SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_REXI_SIMULATION_VARIABLES_HPP_
#define SRC_REXI_SIMULATION_VARIABLES_HPP_

#include <unistd.h>
#include <getopt.h>

/**
 * REXI
 */
struct REXI_SimulationVariables
{
	/**
	 * REXI parameter h
	 */
	double h = 0.15;

	/**
	 * REXI parameter M
	 */
	int M = 128;

	/**
	 * REXI parameter L
	 */
	int L = 0;

	/**
	 * Use only half of the poles for REXI
	 */
	bool use_half_poles = true;

	/**
	 * Extend modes for certain operations
	 */
	int use_extended_modes = 2;

	/**
	 * Normalize REXI for geostrophic balance
	 */
	bool normalization = true;

	/**
	 * Use Next Generation REXI
	 */
	bool use_next_generation = false;

	/**
	 * Use REXI preallocation
	 */
	bool sphere_solver_preallocation = true;

	/**
	 * Use direct solution instead of REXI
	 */
	bool use_direct_solution = false;


	void outputConfig()
	{
		std::cout << std::endl;
		std::cout << "REXI:" << std::endl;
		std::cout << " + rexi_h: " << h << std::endl;
		std::cout << " + rexi_M: " << M << std::endl;
		std::cout << " + rexi_L: " << L << std::endl;
		std::cout << " + rexi_use_half_poles: " << use_half_poles << std::endl;
		std::cout << " + rexi_use_extended_modes: " << use_extended_modes << std::endl;
		std::cout << " + rexi_use_next_generation: " << use_next_generation << std::endl;
		std::cout << " + rexi_normalization: " << normalization << std::endl;
		std::cout << " + rexi_sphere_solver_preallocation: " << sphere_solver_preallocation << std::endl;
		std::cout << " + use_direct_solution: " << use_direct_solution << std::endl;
		std::cout << std::endl;
	}

	void outputProgParams()
	{
		std::cout << "" << std::endl;
		std::cout << "Rexi:" << std::endl;
		std::cout << "	--rexi-h [float]			REXI parameter h" << std::endl;
		std::cout << "	--rexi-m [int]				REXI parameter M" << std::endl;
		std::cout << "	--rexi-l [int]				REXI parameter L" << std::endl;
		std::cout << "	--rexi-half [bool]			Use half REXI poles, default:1" << std::endl;
		std::cout << "	--rexi-normalization [bool]		Use REXI normalization around geostrophic balance, default:1" << std::endl;
		std::cout << "	--rexi-sphere-preallocation [bool]	Use preallocation of SPH-REXI solver coefficients, default:1" << std::endl;
		std::cout << "	--rexi-use-direct-solution [bool]	Use direct solution (analytical) for REXI, default:0" << std::endl;
		std::cout << "	--rexi-use-next-generation [bool]	Use next generation REXI, default:0" << std::endl;
		std::cout << "	--rexi-ext-modes [int]	Use this number of extended modes in spherical harmonics" << std::endl;
		std::cout << "" << std::endl;
	}

	void setup_longOptionList(
			struct option io_long_options[],		///< string and meta information for long options
			int &io_next_free_program_option,	///< number of free options, has to be increased for each new option
			int i_max_options					///< maximum number of options
	)
	{
		io_long_options[io_next_free_program_option] = {"rexi-h", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-m", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-l", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-half", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-normalization", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-sphere-preallocation", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-use-direct-solution", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-use-next-generation", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ext-modes", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;
	}

	/**
	 * Callback method to setup the values for the option with given index.
	 *
	 * \return Number of processed options or 0 in case of processed arguments
	 */
	int setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{
		switch(i_option_index)
		{
			case 0:	h = atof(optarg);	return 0;
			case 1:	M = atoi(optarg);	return 0;
			case 2:	L = atoi(optarg);	return 0;
			case 3:	use_half_poles = atoi(optarg);	return 0;
			case 4:	normalization = atoi(optarg);	return 0;
			case 5:	sphere_solver_preallocation = atoi(optarg);	return 0;
			case 6:	use_direct_solution = atoi(optarg);	return 0;
			case 7:	use_next_generation = atoi(optarg);	return 0;
			case 8:	use_extended_modes = atoi(optarg);	return 0;
		}

		return 9;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
