/*
 * SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SRC_SIMULATION_VARIABLES_HPP_
#define SRC_SIMULATION_VARIABLES_HPP_

#include <unistd.h>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <sweet/StringSplit.hpp>
#include <sweet/SWEETError.hpp>
#include <sweet/TransformationPlans.hpp>
#include <sweet/defaultPrecompilerValues.hpp>
#include <sweet/dict/Dict.hpp>
#include <sweet/ProgramArguments.hpp>
#include "shacks/ShackInterface.hpp"

#if SWEET_THREADING
#include <omp.h>
#endif

#ifndef SWEET_USE_SPHERE_SPECTRAL_SPACE
#	define SWEET_USE_SPHERE_SPECTRAL_SPACE 1
#endif

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
#	include <sweet/sphere/SphereData_Spectral.hpp>
#	include <sweet/sphere/SphereData_Physical.hpp>
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal_SimulationVariables.hpp>
#endif

#if SWEET_LIBPFASST
#       include <libpfasst/LibPFASST_SimulationVariables.hpp>
#endif

#if SWEET_XBRAID
#       include <xbraid/XBraid_SimulationVariables.hpp>
#endif


/*
 * REXI program parameters
 */
#include <rexi/REXI_SimulationVariables.hpp>

/*
 * Polvani program parameters
 */
#include <../programs/swe_plane_benchmarks/SWE_bench_Polvani_SimulationVariables.hpp>

#include "shacksShared/ShackBenchmark.hpp"
#include "shacksShared/ShackDiagnostics.hpp"
#include "shacksShared/ShackDiscretization.hpp"
#include "shacksShared/ShackIOData.hpp"
#include "shacksShared/ShackMisc.hpp"
#include "shacksShared/ShackShackParallelization.hpp"
#include "shacksShared/ShackSDC.hpp"
#include "shacksShared/ShackSimulationCoefficients.hpp"
#include "shacksShared/ShackTimestepControl.hpp"


#include <sweet/ProgramArguments.hpp>

/**
 * program parameters without specific association.
 * These variables can be used differently for each program
 */
class UserDefined
{
public:
	std::string var[20];
};




/**
 * This class exists for convenience reasons.
 *
 * It offers a common structure for the used variables.
 */
class SimulationVariables
{

public:

	Diagnostics diag;
	ShackIOData iodata;
	Benchmark benchmark;
	SimulationCoefficients sim;
	Discretization disc;
	UserDefined user_defined;
	Misc misc;
	ShackParallelization parallelization;
	TimestepControl timecontrol;
	SDC sdc;

	EXP_SimulationVariables rexi;
	SWEPolvani_SimulationVariables swe_polvani;


public:
#if SWEET_PARAREAL
	Parareal_SimulationVariables parareal;
#endif

#if SWEET_LIBPFASST
	LibPFASST_SimulationVariables libpfasst;
#endif

#if SWEET_XBRAID
	XBraid_SimulationVariables xbraid;
#endif


	void outputConfig()
	{
		sim.printClass();
		disc.printClass();
		benchmark.printClass();
		iodata.printClass();
		timecontrol.printClass();

		rexi.printClass();
		sdc.printClass();
		swe_polvani.printClass();
		misc.printClass();
		parallelization.printClass();
		diag.printClass();

#if SWEET_PARAREAL
		parareal.printClass();
#endif

#if SWEET_LIBPFASST
		libpfasst.printClass();
#endif

#if SWEET_XBRAID
		xbraid.printClass();
#endif

	}


	/**
	 * update variables which are based on others
	 */
	void reset()
	{
		if (timecontrol.max_simulation_time < 0)
			SWEETError("timecontrol.max_simulation_time < 0");

		if (timecontrol.max_timesteps_nr < 0)
			SWEETError("timecontrol.max_timesteps_nr < 0");

		timecontrol.current_timestep_nr = 0;
		timecontrol.current_simulation_time = 0;
		timecontrol.current_timestep_size = timecontrol.setup_timestep_size;

		if ((disc.space_res_physical[0] != -1) && (disc.space_res_physical[1] != -1))
			if ((disc.space_res_physical[0] & 1) || (disc.space_res_physical[1] & 1))
				std::cout << "WARNING: Typically there are only even resolutions supported!" << std::endl;

		if (benchmark.random_seed >= 0)
			srandom(benchmark.random_seed);
	}



	void print_params(const char *i_user_defined_program_parameters[])
	{
		sim.printProgramArguments();
		benchmark.printProgramArguments();
		parallelization.printProgramArguments();
		disc.printProgramArguments();
		iodata.printProgramArguments();
		misc.printProgramArguments();
		diag.printProgramArguments();
		timecontrol.printProgramArguments();


		rexi.printProgramArguments();
		sdc.printProgramArguments();
		swe_polvani.printProgramArguments();


#if SWEET_PARAREAL
		parareal.printProgramArguments();
#endif

#if SWEET_LIBPFASST
		libpfasst.printProgramArguments();
#endif

#if SWEET_XBRAID
		xbraid.printProgramArguments();
#endif

		if (i_user_defined_program_parameters != nullptr)
		{
			std::cout << "" << std::endl;
			std::cout << "User defined program parameters:" << std::endl;

			for (int i = 0; i_user_defined_program_parameters[i] != nullptr; i++)
			{
				std::cout << "	--" << i_user_defined_program_parameters[i] << "	\t(see program for description)" << std::endl;
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;
	}


	/**
	 * setup the variables based on program parameters
	 *
	 *
	 * Example for user_defined_program_parameters:
	 * const char *user_defined_program_parameters[] = {{"sweet-file-dict"}, nullptr};
	 */
	bool setupFromMainParameters(
			int i_argc,					///< argc from main()
			char *const i_argv[],		///< argv from main()
			const char *user_defined_program_parameters[] = nullptr,			///< list of strings of simulation-specific program parameters (without --), has to be terminated by nullptr
			bool i_run_prog_parameter_validation = true
	)
	{
		sweet::ProgramArguments programArguments;

		/*
		 * Parse program arguments
		 */
		
		bool argsOK = programArguments.setup(i_argc, i_argv);

		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
			return false;
		
		if (!argsOK)
			SWEETError("ERROR: "+programArguments.error.get());

		sim.processProgramArguments(programArguments);
		benchmark.processProgramArguments(programArguments);
		parallelization.processProgramArguments(programArguments);
		disc.processProgramArguments(programArguments);
		iodata.processProgramArguments(programArguments);
		misc.processProgramArguments(programArguments);
		diag.processProgramArguments(programArguments);
		timecontrol.processProgramArguments(programArguments);

		rexi.processProgramArguments(programArguments);
		sdc.processProgramArguments(programArguments);
		swe_polvani.processProgramArguments(programArguments);


#if SWEET_PARAREAL
		parareal.processProgramArguments(programArguments);
#endif

#if SWEET_LIBPFASST
		libpfasst.processProgramArguments(programArguments);
#endif

#if SWEET_XBRAID
		xbraid.processProgramArguments(programArguments);
#endif
		return true;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
