/*
 * Parameters.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_EXAMPLES_PARAMETERS_HPP_
#define SRC_EXAMPLES_PARAMETERS_HPP_




class Parameters
{
public:
	//// problem size
	std::size_t N = 64;

	// resolution
	std::size_t res[2] = {N,N};

	// size of cell (hx, hy)
	double cell_size[2] = {0,0};

	double h0 = 1000.0;

	// gravitation
	double g = 9.81;

	// cfl condition
	double CFL = 0.005;

	// viscosity
	double viscocity = 0.0;

	// viscosity
	double hyper_viscocity = 0.0;

	// domain length
	double domain_length = 1000;

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

	void setup(
			int i_argc,
			char *i_argv[]
	)
	{
		if (i_argc > 1)
		{
			res[0] = atoi(i_argv[1]);
			res[1] = res[0];
		}

		if (i_argc > 2)
			CFL = atof(i_argv[2]);

		if (i_argc > 3)
			viscocity = atof(i_argv[3]);

		if (i_argc > 4)
			hyper_viscocity = atof(i_argv[4]);

		if (i_argc > 5)
			setup_scenario = atoi(i_argv[5]);

		if (i_argc > 6)
			domain_length = atof(i_argv[6]);

		if (i_argc > 7)
			bogus_var0 = atof(i_argv[7]);

		if (i_argc > 8)
			bogus_var1 = atof(i_argv[8]);

		cell_size[0] = domain_length/res[0];
		cell_size[1] = domain_length/res[1];
	}
};






#endif /* SRC_EXAMPLES_PARAMETERS_HPP_ */
