/*
 * Parameters.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_EXAMPLES_PARAMETERS_HPP_
#define SRC_EXAMPLES_PARAMETERS_HPP_

#include <unistd.h>



class Parameters
{
public:
	//// problem size
	std::size_t N = 64;

	// resolution
	std::size_t res[2] = {N,N};

	// size of cell (hx, hy)
	double cell_size[2] = {0,0};

	// average height for initialization
	double h0 = 1000.0;

	// gravitation
	double g = 9.81;

	// cfl condition
	double sim_CFL = 0.005;

	// viscosity
	double sim_viscocity = 0.0;

	// viscosity
	double sim_hyper_viscocity = 0.0;

	// Coriolis term
	double sim_nondim_f = 0.0;

	// Coriolis term
	double sim_dim_f = 0.0;

	// domain length
	double sim_domain_length = 1000;

	// setup scenario
	int setup_scenario = 0;

	// run simulation timestepping
	bool run_simulation = true;

	double timestep_size = -1;

	// number of simulated time steps
	int timestep_nr = 0;

	// time in simulation
	double simulation_time = 0;

	// mass
	double mass = 0;
	// energy
	double energy = 0;
	// potential enstropy
	double potential_entrophy = 0;

	double bogus_var0 = 0;
	double bogus_var1 = 0;

	double init_coord_x = 0.5;
	double init_coord_y = 0.5;

	// id for visualization
	int vis_id = 0;

	int max_timesteps = -1;
	double max_simulation_time = -1;

	int verbosity = 0;

	void setup(
			int i_argc,
			char *i_argv[]
	)
	{
		int opt;
		while ((opt = getopt(i_argc, i_argv, "n:C:u:U:s:l:a:b:f:x:y:t:T:v:H:")) != -1)
		{
			switch (opt)
			{
			case 'n':
				res[0] = atoi(optarg);
				res[1] = res[0];
				break;

			case 'C':
				sim_CFL = atof(optarg);
				break;

			case 't':
				max_simulation_time = atof(optarg);
				break;

			case 'T':
				max_timesteps = atoi(optarg);
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

			case 'l':
				sim_domain_length = atof(optarg);
				break;

			case 'a':
				bogus_var0 = atof(optarg);
				break;

			case 'b':
				bogus_var1 = atof(optarg);
				break;

			case 'x':
				init_coord_x = atof(optarg);
				break;

			case 'y':
				init_coord_y = atof(optarg);
				break;

			case 'f':
				sim_nondim_f = atof(optarg);
				break;

			case 'v':
				verbosity = atoi(optarg);
				break;

			case 'H':
				h0 = atof(optarg);
				break;

			default:
				std::cerr << "Unknown option '" << (char)opt << "'" << std::endl;
				exit(1);
				break;
			}
		}

		cell_size[0] = sim_domain_length/res[0];
		cell_size[1] = sim_domain_length/res[1];


		sim_dim_f = sim_nondim_f/sim_domain_length;
	}
};






#endif /* SRC_EXAMPLES_PARAMETERS_HPP_ */
