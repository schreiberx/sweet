/*
 * PararealSimulationVariables.hpp
 *
 *  Created on: 18 Apr 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_

#if !SWEET_PARAREAL
#error "Parareal not activated"
#endif

#include <getopt.h>
#include <limits>



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
	double max_simulation_time = std::numeric_limits<double>::infinity();

	/**
	 * Coarse time stepping method string
	 */
	std::string coarse_timestepping_method_string = "";

	/**
	 * Coarse time stepping method order
	 */
	int coarse_timestepping_order = 1;

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
		io_long_options[io_next_free_program_option] = {"parareal-coarse-slices", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-convergence-threshold", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-verbosity", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-enabled", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-max-simulation-time", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-coarse-timestepping-method-string", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-coarse-timestepping-method-order", required_argument, 0, (int)256+io_next_free_program_option};
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
		std::cout << "	--parareal-coarse-slices=[int]              Number of coarse time slices (default=-1, auto define)" << std::endl;
		std::cout << "	--parareal-convergence-threshold=[float]    Threshold for convergence test (default=-1)" << std::endl;
		std::cout << "	--parareal-verbosity=[int]                  Verbosity level (default=0)" << std::endl;
		std::cout << "	--parareal-enabled=[0/1]                    Enable Parareal method (default=0)" << std::endl;
		std::cout << "	--parareal-max-simulation-time=[float]      Overall simulation time (default=-1)" << std::endl;
		std::cout << "	--parareal-coarse-timestepping-method-string=[string]	Identifier for coarse time stepping method (default=ln_erk)" << std::endl;
		std::cout << "	--parareal-coarse-timestepping-order=[int]	Order for coarse time stepping method (default=ln_erk)" << std::endl;
		std::cout << std::endl;
	}



	void outputConfig()
	{
		std::cout << std::endl;
		std::cout << "Parareal" << std::endl;
		std::cout << " + enabled (compile time): " << SWEET_PARAREAL << std::endl;
		std::cout << " + enabled (software flag): " << enabled << std::endl;
		std::cout << " + coarse_slices: " << coarse_slices << std::endl;
		std::cout << " + verbosity: " << verbosity << std::endl;
		std::cout << " + convergence_error_threshold: " << convergence_error_threshold << std::endl;
		std::cout << " + max_simulation_time: " << max_simulation_time << std::endl;
		std::cout << " + coarse_timestepping_method_string: " << coarse_timestepping_method_string << std::endl;
		std::cout << " + coarse_timestepping_method_order: " << coarse_timestepping_order << std::endl;
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

		case 5:
			coarse_timestepping_method_string = i_value;
			break;

		case 6:
			coarse_timestepping_method_string = atoi(i_value);
			break;

		default:
			std::cerr << "Unknown long option with id " << i_option_index << std::endl;
			break;
		}
	}


};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_ */
