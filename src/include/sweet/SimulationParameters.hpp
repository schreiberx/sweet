/*
 * Parameters.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_EXAMPLES_PARAMETERS_HPP_
#define SRC_EXAMPLES_PARAMETERS_HPP_

#include <unistd.h>



class SimulationParameters
{
public:
	/**
	 * SETUP
	 *
	 * values and parameters to setup simulations
	 */
	// average height for initialization
	double setup_h0 = 1000.0;

	// setup scenario
	int setup_scenario = 1;

	// radius
	double setup_radius_scale = 1;

	double setup_coord_x = 0.5;
	double setup_coord_y = 0.5;


	/**
	 * SIMULATION PARAMETERS
	 */
	// gravitation
	double sim_g = 9.81;

	// viscosity
	double sim_viscocity = 0.0;

	// viscosity
	double sim_hyper_viscocity = 0.0;

	// potential
	double sim_potential_viscocity = 0.0;

	// potential
	double sim_potential_hyper_viscocity = 0.0;

	// cfl condition
	double sim_CFL = 0.01;

	// Coriolis term
	double sim_f = 0.0;

	bool use_f_array = false;

	// domain length
	double sim_domain_size[2] = {1000.0*1000.0, 1000.0*1000.0};

	double sim_domain_size2_dbl = -1;

	/**
	 * DISCRETIZATION
	 *
	 * resolution / timestepping
	 */
	// resolution
	std::size_t res[2] = {0,0};

	double res2_dbl = -1.0;

	// use leapfrog like update? (predictor / corrector intermixing h and v,u updates)
	bool timestepping_leapfrog_like_update = false;

	// use up/downwinding for the advection of h
	bool timestepping_up_and_downwinding = false;

	// order of Runge-Kutta scheme for time stepping
	double timestepping_runge_kutta_order = 1;

	// size of time step
	double timestepping_timestep_size = 0;

	// size of cell (hx, hy)
	double sim_cell_size[2] = {0,0};


	// mass
	double diagnostics_mass = 0;
	// energy
	double diagnostics_energy = 0;
	// potential enstropy
	double diagnostics_potential_entrophy = 0;


	/**
	 * program parameters without specific association
	 */
	double bogus_var0 = std::numeric_limits<double>::infinity();
	double bogus_var1 = std::numeric_limits<double>::infinity();
	double bogus_var2 = std::numeric_limits<double>::infinity();
	double bogus_var3 = std::numeric_limits<double>::infinity();
	double bogus_var4 = std::numeric_limits<double>::infinity();

	/**
	 * set verbosity of simulation
	 */
	int verbosity = 0;

	bool gui_enabled =
#if SWEET_GUI
			true
#else
			false
#endif
			;


	/**
	 * SIM CONTROL
	 *
	 * Simulation control variables
	 */
	// run simulation timestepping
	bool run_simulation = true;

	// number of simulated time steps
	int status_timestep_nr = 0;
	int max_timesteps_nr = -1;

	// time in simulation
	double status_simulation_timestep_size = -1;
	double status_simulation_time = 0;
	double max_simulation_time = -1;

	/**
	 * visualization - which conserved quantity to visualize
	 */
	// id for visualization
	int vis_id = 0;

	// use spectral differential operators
	bool use_spectral_diffs = false;

	/**
	 * update variables which are based on others
	 */
	void reset()
	{
		sim_cell_size[0] = sim_domain_size[0]/(double)res[0];
		sim_cell_size[1] = sim_domain_size[1]/(double)res[1];

		res2_dbl = res[0]*res[1];
		sim_domain_size2_dbl = sim_domain_size[0]*sim_domain_size[1];

		status_simulation_timestep_size = -1;
		status_timestep_nr = 0;
		status_simulation_time = 0;

		if ((res[0] & 1) || (res[1] & 1))
		{
			std::cout << "Only even resolutions supported!" << std::endl;
		}
	}


	void setup(
			int i_argc,
			char *i_argv[]
	)
	{
		res[0] = 128;
		res[1] = 128;

		int opt;
		while ((opt = getopt(i_argc, i_argv, "N:n:m:C:u:U:s:X:Y:a:b:c:d:e:f:x:y:t:T:v:H:r:R:W:F:S:g:p:P:G:")) != -1)
		{
			switch (opt)
			{
			case 'N':
				res[0] = atoi(optarg);
				res[1] = res[0];
				break;

			case 'n':
				res[0] = atoi(optarg);
				break;

			case 'm':
				res[1] = atoi(optarg);
				break;

			case 'C':
				sim_CFL = atof(optarg);
				break;

			case 'r':
				setup_radius_scale = atof(optarg);
				break;

			case 't':
				max_simulation_time = atof(optarg);
				break;

			case 'T':
				max_timesteps_nr = atoi(optarg);
				break;

			case 'u':
				sim_viscocity = atof(optarg);
				break;

			case 'U':
				sim_hyper_viscocity = atof(optarg);
				break;

			case 'p':
				sim_potential_viscocity = atof(optarg);
				break;

			case 'P':
				sim_potential_hyper_viscocity = atof(optarg);
				break;

			case 's':
				setup_scenario = atoi(optarg);
				break;

			case 'S':
				use_spectral_diffs = atoi(optarg);
				break;

			case 'X':
				sim_domain_size[0] = atof(optarg);
				break;

			case 'Y':
				sim_domain_size[1] = atof(optarg);
				break;

			case 'R':
				timestepping_runge_kutta_order = atoi(optarg);
				break;

			case 'a':
				bogus_var0 = atof(optarg);
				break;

			case 'b':
				bogus_var1 = atof(optarg);
				break;

			case 'c':
				bogus_var2 = atof(optarg);
				break;

			case 'd':
				bogus_var3 = atof(optarg);
				break;

			case 'e':
				bogus_var4 = atof(optarg);
				break;

			case 'x':
				setup_coord_x = atof(optarg);
				break;

			case 'y':
				setup_coord_y = atof(optarg);
				break;

			case 'f':
				sim_f = atof(optarg);
				break;

			case 'G':
				gui_enabled = atoi(optarg);
				break;

			case 'g':
				sim_g = atof(optarg);
				break;

			case 'v':
				verbosity = atoi(optarg);
				break;

			case 'H':
				setup_h0 = atof(optarg);
				break;

			case 'W':
				timestepping_up_and_downwinding = atoi(optarg);
				break;

			case 'F':
				timestepping_leapfrog_like_update = atoi(optarg);
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
						"	-x [float]	x coordinate for setup \\in [0;1]",
						"	-y [float]	y coordinate for setup \\in [0;1]",
						"	-H [float]	average (initial) height of water",
						"	-r [radius]	scale factor of radius for initial condition",
						"",
						"Discretization:",
						"  >Space:",
						"	-N [res]	resolution in x and y direction",
						"	-n [resx]	resolution in x direction",
						"	-m [resx]	resolution in x direction",
						"	-S [0/1]	use spectral derivatives (experimental)",
						"  >Time:",
						"	-W [0/1]	use up- and downwinding",
						"	-F [0/1]	use leapfrog-like algorithm",
						"	-R [1-RKn]	order of Runge-Kutta method",
						"	-C [cfl]	CFL condition",
						"",
						"Control:",
						"	-t [time]	maximum simulation time",
						"	-T [stepnr]	maximum number of time steps",
						"",
						"Misc options",
						"	-a [float]	bogus variable a",
						"	-b [float]	bogus variable a",
						"	-v [int]	verbosity level",
				};

				std::cerr << "Usage information: " << std::endl;
				for (std::size_t i = 0; i < sizeof(help_strings)/sizeof(*help_strings); i++)
					std::cerr << help_strings[i] << std::endl;

				std::cerr << "Unknown option '" << (char)opt << "'" << std::endl;
				exit(1);
				break;
			}
		}

		reset();

		if (verbosity > 1)
		{
			for (int i = 0; i < i_argc; i++)
				std::cout << i_argv[i] << " ";
			std::cout << std::endl;
		}
	}
};






#endif /* SRC_EXAMPLES_PARAMETERS_HPP_ */
