#ifndef SRC_XBRAID_SIMULATION_VARIABLES_HPP_
#define SRC_XBRAID_SIMULATION_VARIABLES_HPP_

#include <unistd.h>
#include <getopt.h>
#include <string>

struct XBraid_SimulationVariables
{

	/**
	 * Max levels for XBraid solver
	 */
	int xbraid_max_levels = 15;

	/**
	 * (Boolean) Whether to skip all work on first down cycle
	 */
	int xbraid_skip = 1;

	/**
	 * Minimum possible coarse grid size
	 */
	int xbraid_min_coarse = 3;

	/**
	 * Number of CF relaxation sweeps on all levels
	 */
	int xbraid_nrelax = 1;

	/**
	 * Number of CF relaxation only for level 0 -- overrides nrelax
	 */
	int xbraid_nrelax0 = -1;

	/**
	 * Halting tolerance
	 */
	double xbraid_tol = 1e-9;

	/**
	 * Halting norm to use (see docstring below)
	 */
	int xbraid_tnorm = 2;

	/**
	 * Coarsening factor
	 */
	int xbraid_cfactor = 2;

	/**
	 * Coarsening factor for only level 0 -- overrides cfactor
	 */
	int xbraid_cfactor0 = -1;

	/**
	 * Maximum number of iterations
	 */
	int xbraid_max_iter = 100;

	/**
	 * (Boolean) Do FMG cycle (1) or V cycle (0)
	 */
	int xbraid_fmg = 0;

	/**
	 * (Boolean) Use user-defined residual
	 */
	int xbraid_res = 0;

	/**
	 * Full storage on leveles >= 'storage'
	 */ 
	int xbraid_storage = -1;

	/**
	 * Level of XBraid printing to the screen
	 */
	int xbraid_print_level = 2;

	/**
	 * Frequency of calls to Access routine: 1 is for only after simulation and only on the fines time-grid;
	 * 					 2 is for every iteration and every level 
	 */
	int xbraid_access_level = 1;

	/**
	 * Run no simulation, only run wrapper tests
	 */
	int xbraid_run_wrapper_tests = 0;

	/**
	 * Do not compute full residual from user routine each iteration
	 */
	int xbraid_fullrnorm = 0;

	/**
	 * Use the solution from sequential time stepping as the intial guess
	 */
	int xbraid_use_seq_soln = 0;

	/**
	 * Use a random initial guess (1) or a zero initial guess (0)
	 */
	int xbraid_use_rand = 1;

	/**
	 * Number of processors in time
	 */
	int xbraid_pt = 1;


	void outputConfig()
	{
		std::cout << std::endl;
		std::cout << "XBraid:" << std::endl;
		std::cout << " + xbraid_max_levels: "           << xbraid_max_levels         << std::endl;
		std::cout << " + xbraid_skip: "                 << xbraid_skip               << std::endl;
		std::cout << " + xbraid_min_coarse: "           << xbraid_min_coarse         << std::endl;
		std::cout << " + xbraid_nrelax: "               << xbraid_nrelax             << std::endl;
		std::cout << " + xbraid_nrelax0: "              << xbraid_nrelax0            << std::endl;
		std::cout << " + xbraid_tol: "                  << xbraid_tol                << std::endl;
		std::cout << " + xbraid_tnorm: "                << xbraid_tnorm              << std::endl;
		std::cout << " + xbraid_cfactor: "              << xbraid_cfactor            << std::endl;
		std::cout << " + xbraid_cfactor0: "             << xbraid_cfactor0           << std::endl;
		std::cout << " + xbraid_max_iter: "             << xbraid_max_iter           << std::endl;
		std::cout << " + xbraid_fmg: "                  << xbraid_fmg                << std::endl;
		std::cout << " + xbraid_res: "                  << xbraid_res                << std::endl;
		std::cout << " + xbraid_storage: "              << xbraid_storage_           << std::endl;
		std::cout << " + xbraid_print_level: "          << xbraid_print_level        << std::endl;
		std::cout << " + xbraid_access_level: "         << xbraid_access_level       << std::endl;
		std::cout << " + xbraid_run_wrapper_tests: "    << xbraid_run_wrapper_tests  << std::endl;
		std::cout << " + xbraid_fullrnorm: "            << xbraid_fullrnorm          << std::endl;
		std::cout << " + xbraid_use_seq_soln: "         << xbraid_use_seq_soln       << std::endl;
		std::cout << " + xbraid_use_rand: "             << xbraid_use_rand           << std::endl;
		std::cout << " + xbraid_pt: "                   << xbraid_pt                 << std::endl;
		std::cout                                                                    << std::endl;
	}

	void printOptions()
	{
		std::cout << ""                                                                                       << std::endl;
		std::cout << "XBraid:"                                                                                << std::endl;
		std::cout << "	--xbraid-max-levels [int]          XBraid parameter max_levels, default: 15"          << std::endl;
		std::cout << "	--xbraid-skip [int]                XBraid parameter skip, default: 1"                 << std::endl;
		std::cout << "	--xbraid-min-coarse [int]          XBraid parameter min_coarse, default: 3"           << std::endl;
		std::cout << "	--xbraid-nrelax [int]              XBraid parameter nrelax, default: 1"               << std::endl;
		std::cout << "	--xbraid-nrelax0 [int]             XBraid parameter nrelax0, default: -1"             << std::endl;
		std::cout << "	--xbraid-tol [float]               XBraid parameter tol, default: 1e-9"               << std::endl;
		std::cout << "	--xbraid-tnorm [int]               XBraid parameter tnorm, default: 2"                << std::endl;
		std::cout << "	--xbraid-cfactor [int]             XBraid parameter cfactor, default: 2"              << std::endl;
		std::cout << "	--xbraid-cfactor0 [int]            XBraid parameter cfactor0, default: -1"            << std::endl;
		std::cout << "	--xbraid-max-iter [int]            XBraid parameter max_iter, default: 100"           << std::endl;
		std::cout << "	--xbraid-fmg [int]                 XBraid parameter fmg, default: 0"                  << std::endl;
		std::cout << "	--xbraid-res [int]                 XBraid parameter res, default: 0"                  << std::endl;
		std::cout << "	--xbraid-storage [int]             XBraid parameter storage, default: -1"             << std::endl;
		std::cout << "	--xbraid-print-level [int]         XBraid parameter print_level, default: 2"          << std::endl;
		std::cout << "	--xbraid-access-level [int]        XBraid parameter access_level, default: 1"         << std::endl;
		std::cout << "	--xbraid-run-wrapper-tests [int]   XBraid parameter run_wrapper_tests, default: 0"    << std::endl;
		std::cout << "	--xbraid-fullrnorm [int]           XBraid parameter fullrnorm, default: 0"            << std::endl;
		std::cout << "	--xbraid-use-seq-soln [int]        XBraid parameter use_seq_soln, default: 0"         << std::endl;
		std::cout << "	--xbraid-use-rand [int]            XBraid parameter use_rand, default: 1"             << std::endl;
		std::cout << "	--xbraid-pt [int]                  XBraid parameter pt, default: 1"                   << std::endl;
		std::cout << ""                                                                                       << std::endl;
	}

	void setup_longOptionList(
					struct option io_long_options[],		///< string and meta information for long options
					int &io_next_free_program_option,	///< number of free options, has to be increased for each new option
					int i_max_options					///< maximum number of options
	)
	{
		io_long_options[io_next_free_program_option] = {"xbraid-max-levels", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-skip", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-min-coarse", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-nrelax", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-nrelax0", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-tol", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-tnorm", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-cfactor", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-cfactor0", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-max-iter", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-fmg", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-res", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-storage", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-print-level", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-access-level", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-run-wrapper-tests", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-fullrnorm", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-use-seq-soln", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-use-rand", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"xbraid-pt", required_argument, 0, 256+io_next_free_program_option};
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
			case 0:  xbraid_max_levels                = atoi(optarg);	return -1;
			case 1:  xbraid_skip                      = atoi(optarg);	return -1;
			case 2:  xbraid_min_coarse                = atoi(optarg);	return -1;
			case 3:  xbraid_nrelax                    = atoi(optarg);	return -1;
			case 4:  xbraid_nrelax0                   = atoi(optarg);	return -1;
			case 5:  xbraid_tol                       = atof(optarg);	return -1;
			case 6:  xbraid_tnorm                     = atoi(optarg);	return -1;
			case 7:  xbraid_cfactor                   = atoi(optarg);	return -1;
			case 8:  xbraid_cfactor0                  = atoi(optarg);	return -1;
			case 9:  xbraid_max_iter                  = atoi(optarg);	return -1;
			case 10: xbraid_fmg                       = atoi(optarg);	return -1;
			case 11: xbraid_res                       = atoi(optarg);	return -1;
			case 12: xbraid_storage                   = atoi(optarg);	return -1;
			case 13: xbraid_print_level               = atoi(optarg);	return -1;
			case 14: xbraid_access_level              = atoi(optarg);	return -1;
			case 15: xbraid_rn_wrapper_tests          = atoi(optarg);	return -1;
			case 16: xbraid_fullrnorm                 = atoi(optarg);	return -1;
			case 17: xbraid_use_seq_soln              = atoi(optarg);	return -1;
			case 18: xbraid_use_rand                  = atoi(optarg);	return -1;
			case 19: xbraid_pt                        = atoi(optarg);	return -1;
		}
		return 20;
	}

};

#endif /* SRC_XBRAID_SIMULATION_VARIABLES_HPP_ */


