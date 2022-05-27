/*
 * PararealSimulationVariables.hpp
 *
 *  Created on: 18 Apr 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_

#if !SWEET_PARAREAL
#error "Parareal not activated"
#endif


#include <unistd.h>
#include <getopt.h>



/**
 * Simulation variables which are specific to Parareal
 */
class Parareal_SimulationVariables
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
#if SWEET_PARAREAL==2
	int coarse_slices_per_proc = -1;
#endif

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
	std::string coarse_timestepping_method = "";

	/**
	 * Coarse time stepping method order
	 */
	int coarse_timestepping_order = 1;

	/**
	 * Coarse time stepping method order
	 */
	int coarse_timestepping_order2 = 1;

	/**
	 * Time step size for coarse propagator
	 * If == -1 then coarse timestep = time slice length
	 */
	double coarse_timestep_size = -1;

	/**
	 * Read reference csv files
	 * in order to compute and store errors
	 * instead of storing parareal iterations
	 */
	bool load_ref_csv_files = false;
	std::string path_ref_csv_files = "";

	/**
	 * Read csv files of fine simulation
	 * in order to compute and store errors
	 * instead of storing parareal iterations
	 */
	bool load_fine_csv_files = false;
	std::string path_fine_csv_files = "";

	/**
	 * Store parareal iterations
	 * May require too much disk storage!
	 */
	bool store_iterations = true;

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
		if (io_next_free_program_option+8 > i_max_options)
		{
			std::cerr << "Max number of program options exceeded" << std::endl;
			exit(-1);
		}

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

		io_long_options[io_next_free_program_option] = {"parareal-coarse-timestepping-method", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-coarse-timestepping-order", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-coarse-timestepping-order2", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-coarse-timestep-size", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-load-ref-csv-files", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-path-ref-csv-files", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-load-fine-csv-files", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-path-fine-csv-files", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"parareal-store-iterations", required_argument, 0, (int)256+io_next_free_program_option};
		io_next_free_program_option++;
	}


	/**
	 * Callback method to setup the values for the option with given index.
	 */
	void printOptions()
	{
		std::cout << std::endl;
		std::cout << "Parareal options:" << std::endl;
		std::cout << "	--parareal-coarse-slices=[int]              Number of coarse time slices (default=-1, auto define)" << std::endl;
		std::cout << "	--parareal-convergence-threshold=[float]    Threshold for convergence test (default=-1)" << std::endl;
		std::cout << "	--parareal-verbosity=[int]                  Verbosity level (default=0)" << std::endl;
		std::cout << "	--parareal-enabled=[0/1]                    Enable Parareal method (default=0)" << std::endl;
		std::cout << "	--parareal-max-simulation-time=[float]      Overall simulation time (default=-1)" << std::endl;
		std::cout << "	--parareal-coarse-timestepping-method=[string]	Identifier for coarse time stepping method (default=ln_erk)" << std::endl;
		std::cout << "	--parareal-coarse-timestepping-order=[int]	Order for coarse time stepping method (default=1)" << std::endl;
		std::cout << "	--parareal-coarse-timestepping-order2=[int]	Order for coarse time stepping method (default=1)" << std::endl;
		std::cout << "	--parareal-coarse-timestep_size=[float]	Time step for coarse propagation (default=-1)" << std::endl;
		std::cout << "	--parareal-load-ref-csv-files=[0/1]	Load physical reference files (default=0)" << std::endl;
		std::cout << "	--parareal-path-ref-csv-files=[0/1]	Path containing ref csv files (default="")" << std::endl;
		std::cout << "	--parareal-load-fine-csv-files=[0/1]	Load physical files of fine simulation (default=0)" << std::endl;
		std::cout << "	--parareal-path-fine-csv-files=[0/1]	Path containing fine csv files (default="")" << std::endl;
		std::cout << "	--parareal-store-iterations=[0/1]	Store physical files at each iteration (default=1)" << std::endl;
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
		std::cout << " + coarse_timestepping_method: " << coarse_timestepping_method << std::endl;
		std::cout << " + coarse_timestepping_method_order: " << coarse_timestepping_order << std::endl;
		std::cout << " + coarse_timestepping_method_order2: " << coarse_timestepping_order2 << std::endl;
		std::cout << " + coarse_timestep_size: " << coarse_timestep_size << std::endl;
		std::cout << " + load_ref_csv_files: " << load_ref_csv_files << std::endl;
		std::cout << " + path_ref_csv_files: " << path_ref_csv_files << std::endl;
		std::cout << " + load_fine_csv_files: " << load_fine_csv_files << std::endl;
		std::cout << " + path_fine_csv_files: " << path_fine_csv_files << std::endl;
		std::cout << " + store_iterations: " << store_iterations << std::endl;
		std::cout << std::endl;
	}



	/**
	 * Callback method to setup the values for the option with given index.
	 *
	 * \return Number of processed options if nothing was found or 0 in case of a found option
	 */
	int setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{
		switch(i_option_index)
		{
		case 0:
			coarse_slices = atoi(i_value);
			return -1;

		case 1:
			convergence_error_threshold = atof(i_value);
			return -1;

		case 2:
			verbosity = atoi(i_value);
			return -1;

		case 3:
			enabled = atoi(i_value);
			return -1;

		case 4:
			max_simulation_time = atof(i_value);
			return -1;

		case 5:
			coarse_timestepping_method = i_value;
			return -1;

		case 6:
			coarse_timestepping_order = atoi(i_value);
			return -1;

		case 7:
			coarse_timestepping_order2 = atoi(i_value);
			return -1;

		case 8:
			coarse_timestep_size = atof(i_value);
			return -1;

		case 9:
			load_ref_csv_files = atoi(i_value);
			return -1;

		case 10:
			path_ref_csv_files = i_value;
			return -1;

		case 11:
			load_fine_csv_files = atoi(i_value);
			return -1;

		case 12:
			path_fine_csv_files = i_value;
			return -1;

		case 13:
			store_iterations = atoi(i_value);
			return -1;

		}

		return 14;
	}


};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_ */
