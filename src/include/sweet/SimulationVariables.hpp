/*
 * SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_SIMULATION_VARIABLES_HPP_
#define SRC_SIMULATION_VARIABLES_HPP_

#include <unistd.h>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>



/**
 * This class exists for convenience reasons.
 *
 * It offers a common structure for the used variables.
 */
class SimulationVariables
{
public:
	/**
	 * Diagnostic variables
	 */
	struct Diagnostics
	{
		/// total mass
		double total_mass = 0;
		/// total energy
		double total_energy = 0;
		/// total potential enstropy
		double total_potential_enstrophy = 0;

	} diag;



public:
	/**
	 * Values and parameters to setup simulations
	 */
	struct Setup
	{
		/// average height for initialization
		double h0 = 1000.0;

		/// setup scenario
		int scenario = 1;

		/// radius
		double radius_scale = 1;

		/// setup coordinate of e.g. radial breaking dam, x-placement \in [0;1]
		double coord_x = 0.5;
		/// setup coordinate of e.g. radial breaking dam, y-placement \in [0;1]
		double coord_y = 0.5;

		/// Frequency multiplier for wave-like scenario
		double initial_freq_x_mul = 2.0;
		double initial_freq_y_mul = 1.0;

		/// filenames of input data for setup (this has to be setup by each application individually)
		std::vector<std::string> input_data_filenames;

		/// use "BINARY;filename1;filename2" to specify that the binary files should be read in binary format
		bool input_data_binary = false;

		void setup_initial_condition_filenames(
				const std::string i_string
		)
		{
			std::size_t last_pos = 0;
			for (std::size_t pos = 0; i_string[pos] != '\0'; pos++)
			{
				if (i_string[pos] != ';')
					continue;

				input_data_filenames.push_back(i_string.substr(last_pos, pos-last_pos));
				last_pos = pos+1;
			}

			input_data_filenames.push_back(i_string.substr(last_pos));

			if (input_data_filenames.size() > 0)
			{
				if (input_data_filenames[0] == "BINARY")
				{
					input_data_filenames.erase(input_data_filenames.begin());
					input_data_binary = true;
				}
			}
		}
	} setup;


	/**
	 * simulation coefficients
	 */
	struct Coefficients
	{
		/// gravitational constant
		double g = 9.81;

		/// For more information on viscosity,
		/// see 13.3.1 "Generic Form of the Explicit Diffusion Mechanism"
		/// in "Numerical Techniques for Global Atmospheric Models"

		/// viscosity-term on velocities with 2nd order diff operator
		double viscosity = 0.0;

		/// hyper viscosity-term on velocities with 4th order diff operator
		int viscosity_order = 2;

#if 0
		/// viscosity-term on velocities with 2nd order diff operator on potential
		double potential_viscosity = 0.0;

		/// hyper viscosity-term on velocities with 4th order diff operator on potential
		double potential_viscosity_order = 0.0;
#endif
		/// CFL condition
		double CFL = 0.05;

		/// Coriolis frequency f0
		double f0 = 0.0;

		/// Beta coefficient for f(y_N) = f0 + y_N*beta
		/// here, y_N is the normalized y coordinate \in [0;1]
		double beta = 0.0;

		/// zero out the v-component at the top and bottom layer
		bool top_bottom_zero_v_velocity = false;

		/// domain size
//		double domain_size[2] = {1000.0*1000.0, 1000.0*1000.0};
		double domain_size[2] = {1.0, 1.0};
	} sim;


	/**
	 * This class stored the discretization-related parameters
	 *
	 * resolution / timestepping
	 */
	struct Discretization
	{
		/// resolution
		std::size_t res[2] = {128, 128};

		/// size of cell (hx, hy)
		/// this is computed based on disc.res and sim.domain_size
		double cell_size[2] = {0,0};

		/// use leapfrog like update? (predictor / corrector intermixing h and v,u updates)
		bool timestepping_leapfrog_like_update = false;

		/// use up/downwinding for the advection of h
		bool timestepping_up_and_downwinding = false;

		/// order of Runge-Kutta scheme for time stepping
		double timestepping_runge_kutta_order = 1;

		// use spectral differential operators
		bool use_spectral_basis_diffs = true;
	} disc;


	/**
	 * program parameters without specific association.
	 * These variables can be used different for each program
	 */
	struct Bogus
	{
		double var[20] =
		{
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity()
		};
	} bogus;


	/**
	 * Miscellaneous variables
	 */
	struct Misc
	{
		/// set verbosity of simulation
		int verbosity = 0;

		/// precision for floating point output to std::cout and std::endl
		int output_floating_point_precision = 12;

		/// activate GUI mode?
		bool gui_enabled = (SWEET_GUI == 0 ? false : true);

		/// output verbose information every given period of simulation time.
		double be_verbose_after_this_simulation_time_period = 0;

		/// prefix of filename for output of data
		std::string output_file_name_prefix;

		/// prefix of filename for output of data
		double output_each_sim_seconds = -1;

		/// Last simulation seconds of output
		double output_next_sim_seconds = 0;

		/// id for visualization
		int vis_id = 0;
	} misc;


	/**
	 * timestepping
	 */
	struct TimestepControl
	{
		/// Continue running simulation timestepping.
		/// This is beneficial to pause simulations if driven interactively.
		bool run_simulation_timesteps = true;

		/// number of simulated time steps
		int current_timestep_nr = 0;

		/// current time step size
		double current_timestep_size = -1;

		/// time in simulation
		double current_simulation_time = 0;

		/// maximum number of time steps to simulate
		int max_timesteps_nr = -1;

		/// maximum simulation time to execute the simulation for
		double max_simulation_time = -1;
	} timecontrol;



	/**
	 * update variables which are based on others
	 */
	void reset()
	{
		disc.cell_size[0] = sim.domain_size[0]/(double)disc.res[0];
		disc.cell_size[1] = sim.domain_size[1]/(double)disc.res[1];

		timecontrol.current_timestep_size = -1;
		timecontrol.current_timestep_nr = 0;
		timecontrol.current_simulation_time = 0;

		if ((disc.res[0] & 1) || (disc.res[1] & 1))
			std::cout << "WARNING: Typically there are only even resolutions supported!" << std::endl;
	}



	/**
	 * setup the variables based on program parameters
	 */
	bool setupFromMainParameters(
			int i_argc,				///< argc from main()
			char *i_argv[],			///< argv from main()
			const char *bogus_var_names[] = nullptr			///< list of strings of simulation-specific variables, has to be terminated by nullptr
	)
	{
		int next_free_program_option = 2;
		const int max_options = 30;
        static struct option long_options[max_options+1] = {
    			{"test-initial-freq-x-mul", required_argument, 0, 256+'a'+0}, // 0
    			{"test-initial-freq-y-mul", required_argument, 0, 256+'a'+1}, // 1
    			{"initial-coord-x", required_argument, 0, 256+'a'+2}, // 2
    			{"initial-coord-y", required_argument, 0, 256+'a'+3}, // 3
    			{0, 0, 0, 0}, // 4
    			{0, 0, 0, 0}, // 5
    			{0, 0, 0, 0}, // 6
    			{0, 0, 0, 0}, // 7
    			{0, 0, 0, 0}, // 8
    			{0, 0, 0, 0}, // 9

				{0, 0, 0, 0}, // 0
				{0, 0, 0, 0}, // 1
				{0, 0, 0, 0}, // 2
				{0, 0, 0, 0}, // 3
				{0, 0, 0, 0}, // 4
				{0, 0, 0, 0}, // 5
				{0, 0, 0, 0}, // 6
				{0, 0, 0, 0}, // 7
				{0, 0, 0, 0}, // 8
				{0, 0, 0, 0}, // 9

				{0, 0, 0, 0}, // 0
				{0, 0, 0, 0}, // 1
				{0, 0, 0, 0}, // 2
				{0, 0, 0, 0}, // 3
				{0, 0, 0, 0}, // 4
				{0, 0, 0, 0}, // 5
				{0, 0, 0, 0}, // 6
				{0, 0, 0, 0}, // 7
				{0, 0, 0, 0}, // 8
				{0, 0, 0, 0}, // 9	Option Nr. 30
				{0, 0, 0, 0} // NULL
        };


        if (bogus_var_names != nullptr)
        {
			int opt_nr;
			for (opt_nr = next_free_program_option; opt_nr < max_options; opt_nr++)
			{
				if (bogus_var_names[opt_nr-next_free_program_option] == nullptr)
					break;

				long_options[opt_nr].name = bogus_var_names[opt_nr-next_free_program_option];
				long_options[opt_nr].has_arg = required_argument;
				long_options[opt_nr].flag = 0;
				long_options[opt_nr].val = 256+'a'+opt_nr;
			}

			if (opt_nr == max_options)
			{
				std::cerr << "Max number of arguments reached. Reduce number of program arguments" << std::endl;
				exit(1);
			}
        }

		// index into long_options for determined argument
		int option_index = 0;


		int opt;
		while (1)
		{
			opt = getopt_long(	i_argc, i_argv,
							"N:n:m:C:u:U:s:X:Y:f:b:x:y:t:i:T:v:V:O:o:H:r:R:W:F:S:g:p:P:G:d:z",
							long_options, &option_index
					);

			if (opt == -1)
				break;

			/*
			 * LONG OPTIONS
			 */
			if (opt >= 256+'a' && opt <= 256+'z')
			{
				int i = opt-(256+'a');

				if (i < next_free_program_option)
				{
					switch(i)
					{
					case 0:		setup.initial_freq_x_mul = atof(optarg);	break;
					case 1:		setup.initial_freq_y_mul = atof(optarg);	break;
					case 2:		setup.coord_x = atof(optarg);	break;
					case 3:		setup.coord_y = atof(optarg);	break;
					}
				}
				else
				{
					bogus.var[i-next_free_program_option] = atof(optarg);
				}
				continue;
			}

			switch (opt)
			{
			/*
			 * SHORT OPTIONS
			 */
			case 'd':
				misc.output_floating_point_precision = atoi(optarg);
				break;

			case 'N':
				disc.res[0] = atoi(optarg);
				disc.res[1] = disc.res[0];
				break;

			case 'n':
				disc.res[0] = atoi(optarg);
				break;

			case 'm':
				disc.res[1] = atoi(optarg);
				break;

			case 'C':
				sim.CFL = atof(optarg);
				break;

			case 'r':
				setup.radius_scale = atof(optarg);
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
#if 0
			case 'p':
				sim.potential_viscosity = atof(optarg);
				break;

			case 'P':
				sim.potential_viscosity_order = atof(optarg);
				break;
#endif
			case 's':
				setup.scenario = atoi(optarg);
				break;

			case 'S':
				disc.use_spectral_basis_diffs = atoi(optarg);
				break;

			case 'X':
				sim.domain_size[0] = atof(optarg);
				break;

			case 'Y':
				sim.domain_size[1] = atof(optarg);
				break;

			case 'x':
				setup.coord_x = atof(optarg);
				break;

			case 'y':
				setup.coord_y = atof(optarg);
				break;

			case 'f':
				sim.f0 = atof(optarg);
				break;

			case 'b':
				sim.beta = atof(optarg);
				break;

			case 'z':
				sim.top_bottom_zero_v_velocity = true;
				break;

			case 'G':
				misc.gui_enabled = atoi(optarg);
				break;

			case 'g':
				sim.g = atof(optarg);
				break;

			case 'v':
				misc.verbosity = atoi(optarg);
				break;

			case 'V':
				misc.be_verbose_after_this_simulation_time_period = atof(optarg);
				break;

			case 'O':
				misc.output_file_name_prefix = optarg;
				break;

			case 'o':
				misc.output_each_sim_seconds = atof(optarg);
				break;

			case 'H':
				setup.h0 = atof(optarg);
				break;

			case 'R':
				disc.timestepping_runge_kutta_order = atoi(optarg);
				break;

			case 'W':
				disc.timestepping_up_and_downwinding = atoi(optarg);
				break;

			case 'i':
				setup.setup_initial_condition_filenames(optarg);
				break;

			case 'F':
				disc.timestepping_leapfrog_like_update = atoi(optarg);
				std::cout << "WARNING: This time stepping method produces significant errors!" << std::endl;
				std::cerr << "WARNING: This time stepping method produces significant errors!" << std::endl;
				break;

			default:
				const char *help_strings[] = {
						"Simulation runtime parameters",
						"	-X [length]	length of simulation domain in x direction",
						"	-Y [width]	width of simulation domain in y direction",
						"	-u [visc]	viscosity",
						"	-U [visc]	hyperviscosity",
						"	-p [visc]	potential viscosity",
						"	-P [visc]	potential hyperviscosity",
						"	-f [float]	f-parameter for f-plane",
						"	-g [float]	gravity",
						"",
						"Simulation setup parameters",
						"	-s [scen]	scenario id",
						"	            0: radial dam break",
						"	            1: Gaussian dam break",
						"	            2: balanced state x",
						"	            3: balanced state y",
						"	            9: h=H0, v=0, u=0",
						"	-x [float]	x coordinate for setup \\in [0;1]",
						"	-y [float]	y coordinate for setup \\in [0;1]",
						"	-H [float]	average (initial) height of water",
						"	-r [radius]	scale factor of radius for initial condition",
						"",
						"Discretization:",
						"  >Space:",
						"	-N [res]	resolution in x and y direction",
						"	-n [resx]	resolution in x direction",
						"	-m [resy]	resolution in y direction",
						"	-S [0/1]	Control Operator discretization for DataArrays",
						"               0: FD, 1: spectral derivatives, default:0",
						"  >Time:",
						"	-W [0/1]	use up- and downwinding, default:0",
						"	-F [0/1]	use leapfrog-like algorithm, default:0",
						"	-R [1-RKn]	order of Runge-Kutta method, default:1",
						"	-C [cfl]	CFL condition, use negative value for fixed time step size",
						"",
						"Control:",
						"	-t [time]	maximum simulation time",
						"	-T [stepnr]	maximum number of time steps",
						"",
						"Misc options",
						"	-v [int]	verbosity level",
						"	-V [double]	period of output",
						"	-G [0/1]	graphical user interface",
						"	-O [string]	string prefix for filename of output of simulation data",
						"	-d [int]	accuracy of floating point output",
						"	-i [file0][;file1][;file3]...	string with filenames for initial conditions",
						"	            specify BINARY; as first file name to read files as binary raw data",
				};

				std::cerr << "Usage information: " << std::endl;
				for (std::size_t i = 0; i < sizeof(help_strings)/sizeof(*help_strings); i++)
					std::cerr << help_strings[i] << std::endl;

				std::cerr << "Unknown option '" << (char)opt << "'" << std::endl;
				return false;
			}
		}

		reset();

		if (misc.verbosity > 1)
		{
			for (int i = 0; i < i_argc; i++)
				std::cout << i_argv[i] << " ";
			std::cout << std::endl;
		}

		/*
		 * WARNING: the precision of std::cout and std::cerr is set here.
		 * This is not related to the simulation variables but makes it very convenient
		 * to specify it in all other programs.
		 */

		std::cout << std::setprecision(misc.output_floating_point_precision);
		std::cerr << std::setprecision(misc.output_floating_point_precision);

		return true;
	}
};






#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
