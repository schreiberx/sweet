/*
 * PararealSimulationVariables.hpp
 *
 *  Created on: 18 Apr 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_

#if !SWEET_PARAREAL
#error "Parareal not activated"
#endif

#include <getopt.h>



/**
 * Simulation variables which are specific to Parareal
 */
class PararealSimulationVariables
{
public:
	/**
	 * Is parareal enabled?
	 */
	bool enabled = false;

	/**
	 * Number of coarse slices for Parareal.
	 *
	 * If set to -1, the number of coarse slices is automatically determined.
	 */
	int coarse_slices = -1;

	/**
	 * Verbosity of Parareal controller
	 */
	int verbosity = 0;

	/**
	 * Convergence error threshold
	 */
	double convergence_error_threshold = -1;

	/**
	 * Maximum simulation time
	 */
	double max_simulation_time = -1;


	/**
	 * setup long options for program arguments
	 */
public:
	void setup_longOptionList(
			struct option io_long_options[],		///< string and meta information for long options
			int &io_next_free_program_option,	///< number of free options, has to be increased for each new option
			int i_max_options					///< maximum number of options
	)
	{
		io_long_options[io_next_free_program_option] = {"parareal-coarse-slices", required_argument, 0, (int)256+'a'+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-convergence-threshold", required_argument, 0, (int)256+'a'+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-verbosity", required_argument, 0, (int)256+'a'+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-enabled", required_argument, 0, (int)256+'a'+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-max-simulation-time", required_argument, 0, (int)256+'a'+io_next_free_program_option};
		io_next_free_program_option++;

		if (io_next_free_program_option > i_max_options)
		{
			std::cerr << "Max number of program options exceeded" << std::endl;
			exit(-1);
		}
	}


	/**
	 * Callback method to setup the values for the option with given index.
	 */
	void setup_printOptions()
	{
		std::cout << std::endl;
		std::cout << "Parareal options:" << std::endl;
		std::cout << "	--parareal-coarse-slices=[int]				Number of coarse time slices" << std::endl;
		std::cout << "	--parareal-convergence-threshold=[float]	Threshold for convergence test" << std::endl;
		std::cout << "	--parareal-verbosity=[int]					Verbosity level" << std::endl;
		std::cout << "	--parareal-enabled=[0/1]					Enable Parareal method" << std::endl;
		std::cout << "	--parareal-max-simulation-time=[float]		Overall simulation time" << std::endl;
		std::cout << std::endl;
	}


	/**
	 * Callback method to setup the values for the option with given index.
	 */
	void setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{
		switch(i_option_index)
		{
		case 0:
			coarse_slices = atoi(i_value);
			break;


		case 1:
			convergence_error_threshold = atof(i_value);
			break;

		case 2:
			verbosity = atoi(i_value);
			break;

		case 3:
			enabled = atoi(i_value);
			break;

		case 4:
			max_simulation_time = atof(i_value);
			break;

		default:
			std::cerr << "Unknown long option with id " << i_option_index << std::endl;
			break;
		}
	}


};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_ */
