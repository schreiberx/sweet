/*
 * PararealSimulationVariables.hpp
 *
 *  Created on: 18 Apr 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONVARIABLES_HPP_

#if !SWEET_PARAREAL
#error "Parareal not activated"
#endif


#include <unistd.h>
#include <getopt.h>

#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>


/**
 * Simulation variables which are specific to Parareal
 */
class Parareal_SimulationVariables	:
	public sweet::ClassDictionaryInterface
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
	 * Spatial coarsening between the fine and coarse levels
	 */
	bool spatial_coarsening = false;

	/**
	 * Maximum number of parareal iterations
	 */
	int max_iter = -1;
public:

	void printProgramArguments(const std::string& i_prefix = "")
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
		std::cout << "	--parareal-spatial-coarsening=[0/1]	Spatial coarsening between the fine and coarse levels (default=0)" << std::endl;
		std::cout << "	--parareal-max-iter=[int]	Maximum number of parareal iterations (default=-1)" << std::endl;
		std::cout << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--parareal-coarse-slices", coarse_slices);
		i_pa.getArgumentValueByKey("--parareal-convergence-threshold", convergence_error_threshold);
		i_pa.getArgumentValueByKey("--parareal-verbosity", verbosity);
		i_pa.getArgumentValueByKey("--parareal-enabled", enabled);
		i_pa.getArgumentValueByKey("--parareal-max-simulation-time", max_simulation_time);
		i_pa.getArgumentValueByKey("--parareal-coarse-timestepping-method", coarse_timestepping_method);
		i_pa.getArgumentValueByKey("--parareal-coarse-timestepping-order", coarse_timestepping_order);
		i_pa.getArgumentValueByKey("--parareal-coarse-timestepping-order2", coarse_timestepping_order2);
		i_pa.getArgumentValueByKey("--parareal-coarse-timestep-size", coarse_timestep_size);
		i_pa.getArgumentValueByKey("--parareal-load-ref-csv-files", load_ref_csv_files);
		i_pa.getArgumentValueByKey("--parareal-path-ref-csv-files", path_ref_csv_files);
		i_pa.getArgumentValueByKey("--parareal-load-fine-csv-files", load_fine_csv_files);
		i_pa.getArgumentValueByKey("--parareal-path-fine-csv-files", path_fine_csv_files);
		i_pa.getArgumentValueByKey("--parareal-store-iterations", store_iterations);
		i_pa.getArgumentValueByKey("--parareal-spatial-coarsening", spatial_coarsening);
		i_pa.getArgumentValueByKey("--parareal-max-iter", max_iter);

		if (max_simulation_time <= 0)
			return error.set("You need to use --parareal-max-simulation-time with parareal");

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
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
		std::cout << " + spatial coarsening: " << spatial_coarsening << std::endl;
		std::cout << " + max_iter: " << max_iter << std::endl;
		std::cout << std::endl;
	}
};


#endif
