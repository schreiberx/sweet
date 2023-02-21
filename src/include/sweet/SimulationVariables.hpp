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
#include "shacksShared/ShackShackIOData.hpp"
#include "shacksShared/ShackMisc.hpp"
#include "shacksShared/ShackShackParallelization.hpp"
#include "shacksShared/ShackSDC.hpp"
#include "shacksShared/ShackSimulationCoefficients.hpp"
#include "shacksShared/ShackTimestepControl.hpp"


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
#if SWEET_PARAREAL
	Parareal_SimulationVariables parareal;
#endif

#if SWEET_LIBPFASST
	LibPFASST_SimulationVariables libpfasst;
#endif

#if SWEET_XBRAID
	XBraid_SimulationVariables xbraid;
#endif


public:
	EXP_SimulationVariables rexi;
	SWEPolvani_SimulationVariables swe_polvani;


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



	void outputConfig()
	{
		sim.outputConfig();
		disc.outputConfig();
		benchmark.outputConfig();
		iodata.outputConfig();
		timecontrol.outputConfig();

		rexi.outputConfig();
		sdc.outputConfig();
		swe_polvani.outputConfig();
		misc.outputConfig();
		parallelization.outputConfig();
		diag.outputConfig();

#if SWEET_PARAREAL
		parareal.outputConfig();
#endif

#if SWEET_LIBPFASST
		libpfasst.outputConfig();
#endif

#if SWEET_XBRAID
		xbraid.outputConfig();
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

		sim.outputProgParams();
		benchmark.outputProgParams();
		parallelization.outputProgParams();
		disc.outputProgParams();
		iodata.outputProgParams();


		misc.outputProgParams();
		rexi.outputProgParams();
		sdc.outputProgParams();
		swe_polvani.outputProgParams();


#if SWEET_PARAREAL
		parareal.printOptions();
#endif

#if SWEET_LIBPFASST
		libpfasst.printOptions();
#endif

#if SWEET_XBRAID
		xbraid.printOptions();
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
		const int max_options = 200;
		struct option long_options[max_options+1];

		for (std::size_t i = 0; i < max_options+1; i++)
		{
			long_options[i].flag = 0;
			long_options[i].has_arg = 0;
			long_options[i].name = 0;
			long_options[i].val = 0;
		}

		int next_free_program_option = 0;

        int benchmark_start_option_index = next_free_program_option;
		benchmark.setup_longOptionsList(long_options, next_free_program_option);

        int sim_start_option_index = next_free_program_option;
		sim.setup_longOptionsList(long_options, next_free_program_option);

        int iodata_start_option_index = next_free_program_option;
		iodata.setup_longOptionsList(long_options, next_free_program_option);

        int misc_start_option_index = next_free_program_option;
		misc.setup_longOptionsList(long_options, next_free_program_option);

        int parallelization_start_option_index = next_free_program_option;
        parallelization.setup_longOptionsList(long_options, next_free_program_option);

        int disc_start_option_index = next_free_program_option;
		disc.setup_longOptionsList(long_options, next_free_program_option);

        int timecontrol_start_option_index = next_free_program_option;
		timecontrol.setup_longOptionsList(long_options, next_free_program_option);

#if SWEET_PARAREAL
        int parareal_start_option_index = next_free_program_option;
        parareal.setup_longOptionList(
        		long_options,
				next_free_program_option,	///< also updated (IO)
				max_options
			);
#endif

#if SWEET_LIBPFASST
        int libpfasst_start_option_index = next_free_program_option;
        libpfasst.setup_longOptionList(
				       long_options,
				       next_free_program_option,	///< also updated (IO)
				       max_options
				       );
#endif

#if SWEET_XBRAID
        int xbraid_start_option_index = next_free_program_option;
        xbraid.setup_longOptionList(
				       long_options,
				       next_free_program_option,	///< also updated (IO)
				       max_options
				       );
#endif


        int rexi_start_option_index = next_free_program_option;
        rexi.setup_longOptionList(
        		long_options,
				next_free_program_option,	///< also updated (IO)
				max_options
			);

		int sdc_start_option_index = next_free_program_option;
		sdc.setup_longOptionsList(long_options, next_free_program_option);

        int swe_polvani_start_option_index = next_free_program_option;
        swe_polvani.setup_longOptionList(
        		long_options,
				next_free_program_option,	///< also updated (IO)
				max_options
			);

        // Test dummy object
        long_options[next_free_program_option] = {"dummy", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        if (user_defined_program_parameters != nullptr)
        {
			int opt_nr;
			for (opt_nr = next_free_program_option; opt_nr < max_options; opt_nr++)
			{
				if (user_defined_program_parameters[opt_nr-next_free_program_option] == nullptr)
					break;

				long_options[opt_nr].name = user_defined_program_parameters[opt_nr-next_free_program_option];
				long_options[opt_nr].has_arg = required_argument;
				long_options[opt_nr].flag = 0;
				long_options[opt_nr].val = 256+opt_nr;
			}

			if (opt_nr == max_options)
			{
				SWEETError("Max number of arguments reached. Reduce number of program arguments");
			}
        }

		// index into long_options for argument to be determined
		int option_index = 0;

		int opt;
		while (1)
		{
			opt = getopt_long(
				i_argc, i_argv,
				"N:M:n:m:u:U:s:X:Y:f:F:b:x:y:t:i:T:v:V:O:o:H:r:a:R:W:F:S:g:G:d:zh",
				long_options, &option_index
			);

			if (opt == -1)
				break;

			if (opt == '?')
			{
				//print_params();
				std::cerr << std::endl;
				std::cerr << "Error while processing program arguments (see above)" << std::endl;
				std::cerr << std::endl;
				std::cerr << "Please use -h option as first argument to see available parameters" << std::endl;
				std::cerr << std::endl;
				return false;
			}

			/*
			 * LONG OPTIONS?
			 */
			if (opt >= 256)
			{
				int i = opt-256;

				if (i < next_free_program_option)
				{
					int c = 0;

					{
						int retval = benchmark.setup_longOptionValue(i-benchmark_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					{
						int retval = sim.setup_longOptionValue(i-sim_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					{
						int retval = iodata.setup_longOptionValue(i-iodata_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}


					{
						int retval = misc.setup_longOptionValue(i-misc_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					{
						int retval = parallelization.setup_longOptionValue(i-parallelization_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					{
						int retval = disc.setup_longOptionValue(i-disc_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					{
						int retval = timecontrol.setup_longOptionValue(i-timecontrol_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

#if SWEET_PARAREAL
					{
						int retval = parareal.setup_longOptionValue(i-parareal_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}
#endif

#if SWEET_LIBPFASST
					{
						int retval = libpfasst.setup_longOptionValue(i-libpfasst_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}
#endif

#if SWEET_XBRAID
					{
						int retval = xbraid.setup_longOptionValue(i-xbraid_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}
#endif


					{
						int retval = rexi.setup_longOptionValue(i-rexi_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					{
						int retval = sdc.setup_longOptionValue(i-sdc_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					{
						int retval = swe_polvani.setup_longOptionValue(i-swe_polvani_start_option_index, optarg);
						if (retval == -1)
							continue;
						c += retval;
					}

					c++;

					/*
					 * This can be tested with the --dummy parameter
					 */
					if (c != next_free_program_option-1)
					{
						outputConfig();
						std::cout << "TEST TEST" << std::endl;
						std::cout << (int)c << std::endl;
						std::cout << (int)next_free_program_option-1 << std::endl;
						SWEETError("Inconsistent processing of arguments");
					}
				}
				else
				{
					int user_defined_id = i-next_free_program_option;

					if (user_defined_id >= max_options)
					{
						std::cout << std::endl;
						std::cout << "SERIOUS ERROR" << std::endl;
						std::cout << " + long option: " << i_argv[option_index] << std::endl;
						std::cout << " + user_defined_id " << user_defined_id << std::endl;
						std::cout << std::endl;
						exit(1);
					}
					user_defined.var[i-next_free_program_option] = optarg;
				}
				continue;
			}

			// short options from hereon
			if (optarg != nullptr)
			{
				if (optarg[0] == '=')
				{
					std::cerr << "Short option parameters may not be specified with an equal '=' sign!" << std::endl;
					SWEETError("Exit");
				}
			}

			switch (opt)
			{
			/*
			 * SHORT OPTIONS
			 */
			case 'd':
				iodata.output_floating_point_precision = atoi(optarg);
				break;

			case 'N':
				{
					int c = StringSplit::split2int(optarg, &disc.space_res_physical[0], &disc.space_res_physical[1]);
					if (c == 1)
						disc.space_res_physical[1] = disc.space_res_physical[0];
				}
				break;

			case 'M':
				{
					int c = StringSplit::split2int(optarg, &disc.space_res_spectral[0], &disc.space_res_spectral[1]);
					if (c == 1)
						disc.space_res_spectral[1] = disc.space_res_spectral[0];
				}
				break;

			case 'n':
				disc.space_res_physical[0] = atoi(optarg);
				break;

			case 'm':
				disc.space_res_physical[1] = atoi(optarg);
				break;

			case 'r':
				benchmark.object_scale = atof(optarg);
				break;

			case 't':
				timecontrol.max_simulation_time = atof(optarg);
				break;

			case 'T':
				timecontrol.max_timesteps_nr = atoi(optarg);
				break;

			case 'u':
				sim.viscosity = atof(optarg);
				break;

			case 'U':
				sim.viscosity_order = atoi(optarg);
				break;

			case 'S':
				disc.space_use_spectral_basis_diffs = atoi(optarg);
				break;

			case 'X':
				sim.plane_domain_size[0] = atof(optarg);
				break;

			case 'Y':
				sim.plane_domain_size[1] = atof(optarg);
				break;

			case 'x':
				benchmark.object_coord_x = atof(optarg);
				break;

			case 'y':
				benchmark.object_coord_y = atof(optarg);
				break;

			case 'f':
				sim.plane_rotating_f0 = atof(optarg);
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
				sim.sphere_rotating_coriolis_omega = atof(optarg);
				sim.sphere_fsphere_f0 = atof(optarg);
#endif
				break;

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
			case 'F':
				sim.sphere_use_fsphere = atoi(optarg);
				break;

			case 'a':
				sim.sphere_radius = atof(optarg);
				break;
#endif

			case 'G':
				misc.gui_enabled = atoi(optarg);
				break;

			case 'g':
				sim.gravitation = atof(optarg);
				break;

			case 'v':
				misc.verbosity = atoi(optarg);
				break;

			case 'O':
				iodata.output_file_name = optarg;
				if (iodata.output_file_name == "-")
					iodata.output_file_name = "";
				break;

			case 'o':
				iodata.output_each_sim_seconds = atof(optarg);
				break;

			case 'H':
				sim.h0 = atof(optarg);
				break;

			case 'R':
				disc.timestepping_order = atoi(optarg);
				break;

			case 'i':
				iodata.setup_initial_condition_filenames(optarg);
				break;

			case 'h':
				print_params(user_defined_program_parameters);
				return false;

			default:
				print_params(user_defined_program_parameters);

				std::cerr << "The option '-";
				std::cerr << (char)opt;
				std::cerr << "' was specified to be available, but it's parameter detection is not implemented." << std::endl;
				std::cerr << "Please contact the SWEET developer" << std::endl;

				SWEETError("Exit");
				return false;
			}
		}



		if (i_run_prog_parameter_validation)
		{
			if (	(disc.space_res_physical[0] == 0 || disc.space_res_physical[1] == 0)	&&
					(disc.space_res_spectral[0] == 0 || disc.space_res_spectral[1] == 0)
			)
			{
				SWEETError("Select physical resolution or spectral modes (use -N (or -n, -m) for physical and -M for spectral) ");
			}

			if (iodata.output_file_mode == "default")
			{
#if 1
				iodata.output_file_mode = "bin";

				if (iodata.output_file_name == "X")
					iodata.output_file_name = "output_%s_t%020.8f.sweet";
#else
				iodata.output_file_mode = "csv";

				if (iodata.output_file_name == "X")
					iodata.output_file_name = "output_%s_t%020.8f.csv";
#endif
			}
			else
			{
				if (iodata.output_file_name == "X")
				{
					if (iodata.output_file_mode == "csv")
						iodata.output_file_name = "output_%s_t%020.8f.csv";
					else if (iodata.output_file_mode == "bin")
						iodata.output_file_name = "output_%s_t%020.8f.sweet";
					else if (iodata.output_file_mode == "csv_spec_evol")
						iodata.output_file_name = "output_%s_t%020.8f.txt";
					else
						SWEETError("Unknown filemode '"+iodata.output_file_mode+"'");
				}
			}
		}

		reset();

#if SWEET_PARAREAL
		// if max simulation time was not set for parareal, copy max simulation time from default parameters to parareal parameters.
		if (parareal.max_simulation_time <= 0)
			parareal.max_simulation_time = timecontrol.max_simulation_time;
#endif

		/*
		 * WARNING: the precision of std::cout and std::cerr is set here.
		 * This is not related to the simulation variables but makes it very convenient
		 * to specify it in all other programs.
		 */

		if (iodata.output_file_name == "-")
			iodata.output_file_name = "";

		if (iodata.output_floating_point_precision >= 0)
		{
			std::cout << std::setprecision(iodata.output_floating_point_precision);
			std::cerr << std::setprecision(iodata.output_floating_point_precision);
		}

		if (disc.timestepping_order2 <= 0)
		{
			disc.timestepping_order2 = disc.timestepping_order;
		}

		return true;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
