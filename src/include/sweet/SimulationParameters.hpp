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
	int setup_scenario = 0;

	// radius
	double setup_radius = 100000;

	double setup_coord_x = 500.0*1000.0;
	double setup_coord_y = 500.0*1000.0;


	/**
	 * SIMULATION PARAMETERS
	 */
	// gravitation
	double sim_g = 9.81;

	// viscosity
	double sim_viscocity = 0.0;

	// viscosity
	double sim_hyper_viscocity = 0.0;

	// cfl condition
	double sim_CFL = 0.01;

	// Coriolis term
	double sim_f = 0.0;

	bool use_f_array = false;

	// domain length
	double sim_domain_length[2] = {1000.0*1000.0, 1000.0*1000.0};


	/**
	 * DISCRETIZATION
	 *
	 * resolution / timestepping
	 */
	// resolution
	std::size_t res[2] = {0,0};

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
	double bogus_var0 = 0;
	double bogus_var1 = 0;

	/**
	 * set verbosity of simulation
	 */
	int verbosity = 0;


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
		sim_cell_size[0] = sim_domain_length[0]/res[0];
		sim_cell_size[1] = sim_domain_length[1]/res[1];

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
		while ((opt = getopt(i_argc, i_argv, "N:n:m:C:u:U:s:X:Y:a:b:f:x:y:t:T:v:H:r:R:W:F:S:")) != -1)
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
				setup_radius = atoi(optarg);
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

			case 's':
				setup_scenario = atoi(optarg);
				break;

			case 'S':
				use_spectral_diffs = atoi(optarg);
				break;

			case 'X':
				sim_domain_length[0] = atof(optarg);
				break;

			case 'Y':
				sim_domain_length[1] = atof(optarg);
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

			case 'x':
				setup_coord_x = atof(optarg);
				break;

			case 'y':
				setup_coord_y = atof(optarg);
				break;

			case 'f':
				sim_f = atof(optarg);
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
				std::cerr << "Unknown option '" << (char)opt << "'" << std::endl;
				exit(1);
				break;
			}
		}

		reset();
	}
};






#endif /* SRC_EXAMPLES_PARAMETERS_HPP_ */
