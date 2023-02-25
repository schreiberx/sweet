#ifndef SRC_XBRAID_SIMULATION_VARIABLES_HPP_
#define SRC_XBRAID_SIMULATION_VARIABLES_HPP_

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <sweet/core/ProgramArguments.hpp>
#include "../sweet/shacks/ShackInterface.hpp"

struct XBraid_SimulationVariables	:
	public ShackInterface
{

	/**
	 * Enable XBraid
	 */
	bool xbraid_enabled = false;

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
	 * Number of V-cycles at the end of a F-cycle
	 */
	int xbraid_fmg_vcyc = 0;


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

	/**
	 * Timestepping methods for all levels, separated by comma
	 */
	std::string xbraid_timestepping_method;

	/**
	 * Timestepping orders for all levels, separated by comma
	 */
	std::string xbraid_timestepping_order;

	/**
	 * Timestepping orders for all levels, separated by comma
	 */
	std::string xbraid_timestepping_order2;

	/**
	 * Viscosity order for all levels, separated by comma
	 */
	std::string xbraid_viscosity_order = "2";

	/**
	 * Diffusion coefficient for all levels, separated by comma
	 */
	std::string xbraid_viscosity_coefficient = "0";

	/**
	 * Verbosity of XBraid controller
	 */
	int xbraid_verbosity = 0;

	/**
	 * Read reference csv files
	 * in order to compute and store errors
	 * instead of storing XBraid iterations
	 */
	bool xbraid_load_ref_csv_files = false;
	std::string xbraid_path_ref_csv_files = "";

	/**
	 * Read csv files of fine simulation
	 * in order to compute and store errors
	 * instead of storing XBraid iterations
	 */
	bool xbraid_load_fine_csv_files = false;
	std::string xbraid_path_fine_csv_files = "";

	/**
	 * Store XBraid iterations
	 * May require too much disk storage!
	 */
	bool xbraid_store_iterations = true;

	/**
	 * Spatial coarsen between levels
	 * Proportionally to the timestep size
	 */
	int xbraid_spatial_coarsening = 0;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << ""                                                                                                         << std::endl;
		std::cout << "XBraid:"                                                                                                  << std::endl;
		std::cout << "	--xbraid-enabled [0/1]                       XBraid parameter enabled, default: 0"                      << std::endl;
		std::cout << "	--xbraid-max-levels [int]                    XBraid parameter max_levels, default: 15"                  << std::endl;
		std::cout << "	--xbraid-skip [0/1]                          XBraid parameter skip, default: 1"                         << std::endl;
		std::cout << "	--xbraid-min-coarse [int]                    XBraid parameter min_coarse, default: 3"                   << std::endl;
		std::cout << "	--xbraid-nrelax [int]                        XBraid parameter nrelax, default: 1"                       << std::endl;
		std::cout << "	--xbraid-nrelax0 [int]                       XBraid parameter nrelax0, default: -1"                     << std::endl;
		std::cout << "	--xbraid-tol [float]                         XBraid parameter tol, default: 1e-9"                       << std::endl;
		std::cout << "	--xbraid-tnorm [int]                         XBraid parameter tnorm, default: 2"                        << std::endl;
		std::cout << "	--xbraid-cfactor [int]                       XBraid parameter cfactor, default: 2"                      << std::endl;
		std::cout << "	--xbraid-cfactor0 [int]                      XBraid parameter cfactor0, default: -1"                    << std::endl;
		std::cout << "	--xbraid-max-iter [int]                      XBraid parameter max_iter, default: 100"                   << std::endl;
		std::cout << "	--xbraid-fmg [int]                           XBraid parameter fmg, default: 0"                          << std::endl;
		std::cout << "	--xbraid-fmg-vcyc [int]                      XBraid parameter fmg_vcyc, default: 0"                     << std::endl;
		std::cout << "	--xbraid-res [int]                           XBraid parameter res, default: 0"                          << std::endl;
		std::cout << "	--xbraid-storage [int]                       XBraid parameter storage, default: -1"                     << std::endl;
		std::cout << "	--xbraid-print-level [int]                   XBraid parameter print_level, default: 2"                  << std::endl;
		std::cout << "	--xbraid-access-level [int]                  XBraid parameter access_level, default: 1"                 << std::endl;
		std::cout << "	--xbraid-run-wrapper-tests [0/1]             XBraid parameter run_wrapper_tests, default: 0"            << std::endl;
		std::cout << "	--xbraid-fullrnorm [int]                     XBraid parameter fullrnorm, default: 0"                    << std::endl;
		std::cout << "	--xbraid-use-seq-soln [0/1]                  XBraid parameter use_seq_soln, default: 0"                 << std::endl;
		std::cout << "	--xbraid-use-rand [0/1]                      XBraid parameter use_rand, default: 1"                     << std::endl;
		std::cout << "	--xbraid-pt [int]                            XBraid parameter pt, default: 1"                           << std::endl;
		std::cout << "	--xbraid-timestepping-method [string]        XBraid parameter timestepping-method, default: ''"         << std::endl;
		std::cout << "	--xbraid-timestepping-order [string]         XBraid parameter timestepping-order, default: ''"          << std::endl;
		std::cout << "	--xbraid-timestepping-order2 [string]        XBraid parameter timestepping-order2, default: ''"         << std::endl;
		std::cout << "	--xbraid-viscosity-order [string]            XBraid parameter viscosity-order, default: '2'"            << std::endl;
		std::cout << "	--xbraid-viscosity-coefficient [string]      XBraid parameter viscosity-coefficient, default: '0'"      << std::endl;
		std::cout << "	--xbraid-verbosity [int]                     XBraid parameter verbosity, default: 0"                    << std::endl;
		std::cout << "	--xbraid-load-ref-csv-files [0/1]            XBraid parameter load_ref_csv_files, default: 0"           << std::endl;
		std::cout << "	--xbraid-path-ref-csv-files [string]         XBraid parameter path_ref_csv_files, default: ''"          << std::endl;
		std::cout << "	--xbraid-load-fine-csv-files [0/1]           XBraid parameter load_fine_csv_files, default: 0"          << std::endl;
		std::cout << "	--xbraid-path-fine-csv-files [string]        XBraid parameter path_fine_csv_files, default: ''"         << std::endl;
		std::cout << "	--xbraid-store-iterations [0/1]              XBraid parameter store_iterations, default: 0"             << std::endl;
		std::cout << "	--xbraid-spatial-coarsening [0/1]            XBraid parameter spatial_coarsening, default: 0"           << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--xbraid-enabled", xbraid_enabled);
		i_pa.getArgumentValueByKey("--xbraid-max-levels", xbraid_max_levels);
		i_pa.getArgumentValueByKey("--xbraid-skip", xbraid_skip);
		i_pa.getArgumentValueByKey("--xbraid-min-coarse", xbraid_min_coarse);
		i_pa.getArgumentValueByKey("--xbraid-nrelax", xbraid_nrelax);
		i_pa.getArgumentValueByKey("--xbraid-nrelax0", xbraid_nrelax0);
		i_pa.getArgumentValueByKey("--xbraid-tol", xbraid_tol);
		i_pa.getArgumentValueByKey("--xbraid-tnorm", xbraid_tnorm);
		i_pa.getArgumentValueByKey("--xbraid-cfactor", xbraid_cfactor);
		i_pa.getArgumentValueByKey("--xbraid-cfactor0", xbraid_cfactor0);
		i_pa.getArgumentValueByKey("--xbraid-max-iter", xbraid_max_iter);
		i_pa.getArgumentValueByKey("--xbraid-fmg", xbraid_fmg);
		i_pa.getArgumentValueByKey("--xbraid-fmg-vcyc", xbraid_fmg_vcyc);
		i_pa.getArgumentValueByKey("--xbraid-res", xbraid_res);
		i_pa.getArgumentValueByKey("--xbraid-storage", xbraid_storage);
		i_pa.getArgumentValueByKey("--xbraid-print-level", xbraid_print_level);
		i_pa.getArgumentValueByKey("--xbraid-access-level", xbraid_access_level);
		i_pa.getArgumentValueByKey("--xbraid-run-wrapper-tests", xbraid_run_wrapper_tests);
		i_pa.getArgumentValueByKey("--xbraid-fullrnorm", xbraid_fullrnorm);
		i_pa.getArgumentValueByKey("--xbraid-use-seq-soln", xbraid_use_seq_soln);
		i_pa.getArgumentValueByKey("--xbraid-use-rand", xbraid_use_rand);
		i_pa.getArgumentValueByKey("--xbraid-pt", xbraid_pt);
		i_pa.getArgumentValueByKey("--xbraid-timestepping-method", xbraid_timestepping_method);
		i_pa.getArgumentValueByKey("--xbraid-timestepping-order", xbraid_timestepping_order);
		i_pa.getArgumentValueByKey("--xbraid-timestepping-order2", xbraid_timestepping_order2);
		i_pa.getArgumentValueByKey("--xbraid-viscosity-order", xbraid_viscosity_order);
		i_pa.getArgumentValueByKey("--xbraid-viscosity-coefficient", xbraid_viscosity_coefficient);
		i_pa.getArgumentValueByKey("--xbraid-verbosity", xbraid_verbosity);
		i_pa.getArgumentValueByKey("--xbraid-load-ref-csv-files", xbraid_load_ref_csv_files);
		i_pa.getArgumentValueByKey("--xbraid-path-ref-csv-files", xbraid_path_ref_csv_files);
		i_pa.getArgumentValueByKey("--xbraid-load-fine-csv-files", xbraid_load_fine_csv_files);
		i_pa.getArgumentValueByKey("--xbraid-path-fine-csv-files", xbraid_path_fine_csv_files);
		i_pa.getArgumentValueByKey("--xbraid-store-iterations", xbraid_store_iterations);
		i_pa.getArgumentValueByKey("--xbraid-spatial-coarsening", xbraid_spatial_coarsening);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "XBraid:" << std::endl;
		std::cout << " + xbraid_enabled: "                      << xbraid_enabled                      << std::endl;
		std::cout << " + xbraid_max_levels: "                   << xbraid_max_levels                   << std::endl;
		std::cout << " + xbraid_skip: "                         << xbraid_skip                         << std::endl;
		std::cout << " + xbraid_min_coarse: "                   << xbraid_min_coarse                   << std::endl;
		std::cout << " + xbraid_nrelax: "                       << xbraid_nrelax                       << std::endl;
		std::cout << " + xbraid_nrelax0: "                      << xbraid_nrelax0                      << std::endl;
		std::cout << " + xbraid_tol: "                          << xbraid_tol                          << std::endl;
		std::cout << " + xbraid_tnorm: "                        << xbraid_tnorm                        << std::endl;
		std::cout << " + xbraid_cfactor: "                      << xbraid_cfactor                      << std::endl;
		std::cout << " + xbraid_cfactor0: "                     << xbraid_cfactor0                     << std::endl;
		std::cout << " + xbraid_max_iter: "                     << xbraid_max_iter                     << std::endl;
		std::cout << " + xbraid_fmg: "                          << xbraid_fmg                          << std::endl;
		std::cout << " + xbraid_fmg_vcyc: "                     << xbraid_fmg_vcyc                     << std::endl;
		std::cout << " + xbraid_res: "                          << xbraid_res                          << std::endl;
		std::cout << " + xbraid_storage: "                      << xbraid_storage                      << std::endl;
		std::cout << " + xbraid_print_level: "                  << xbraid_print_level                  << std::endl;
		std::cout << " + xbraid_access_level: "                 << xbraid_access_level                 << std::endl;
		std::cout << " + xbraid_run_wrapper_tests: "            << xbraid_run_wrapper_tests            << std::endl;
		std::cout << " + xbraid_fullrnorm: "                    << xbraid_fullrnorm                    << std::endl;
		std::cout << " + xbraid_use_seq_soln: "                 << xbraid_use_seq_soln                 << std::endl;
		std::cout << " + xbraid_use_rand: "                     << xbraid_use_rand                     << std::endl;
		std::cout << " + xbraid_pt: "                           << xbraid_pt                           << std::endl;
		std::cout << " + xbraid_timestepping_method: "          << xbraid_timestepping_method          << std::endl;
		std::cout << " + xbraid_timestepping_order: "           << xbraid_timestepping_order           << std::endl;
		std::cout << " + xbraid_timestepping_order2: "          << xbraid_timestepping_order2          << std::endl;
		std::cout << " + xbraid_viscosity_order: "              << xbraid_viscosity_order              << std::endl;
		std::cout << " + xbraid_viscosity_coefficient: "        << xbraid_viscosity_coefficient        << std::endl;
		std::cout << " + xbraid_verbosity: "                    << xbraid_verbosity                    << std::endl;
		std::cout << " + xbraid_load_ref_csv_files: "           << xbraid_load_ref_csv_files           << std::endl;
		std::cout << " + xbraid_path_ref_csv_files: "           << xbraid_path_ref_csv_files           << std::endl;
		std::cout << " + xbraid_load_fine_csv_files: "          << xbraid_load_fine_csv_files          << std::endl;
		std::cout << " + xbraid_path_fine_csv_files: "          << xbraid_path_fine_csv_files          << std::endl;
		std::cout << " + xbraid_store_iterations: "             << xbraid_store_iterations             << std::endl;
		std::cout << " + xbraid_spatial_coarsening: "           << xbraid_spatial_coarsening           << std::endl;
	}
};

#endif /* SRC_XBRAID_SIMULATION_VARIABLES_HPP_ */


