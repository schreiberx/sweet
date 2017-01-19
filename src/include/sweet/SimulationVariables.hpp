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
#include <limits>
#include <sweet/sweetmath.hpp>
#include <sweet/FatalError.hpp>


#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#ifndef SWEET_GUI
#	define SWEET_GUI 1
#endif

#ifndef SWEET_PFASST_CPP
#	define SWEET_PFASST_CPP 1
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal_SimulationVariables.hpp>
#endif


/**
 * This class exists for convenience reasons.
 *
 * It offers a common structure for the used variables.
 */
class SimulationVariables
{
public:
#if SWEET_PARAREAL
	PararealSimulationVariables parareal;
#endif



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

		void output()
		{
			std::cout << std::endl;
			std::cout << "DIAGNOSTICS:" << std::endl;
			std::cout << " + total_mass: " << total_mass << std::endl;
			std::cout << " + total_energy: " << total_energy << std::endl;
			std::cout << " + total_potential_enstrophy: " << total_potential_enstrophy << std::endl;
			std::cout << std::endl;
		}
	} diag;



public:
	/**
	 * Values and parameters to setup simulations
	 */
	struct Setup
	{
		/// setup scenario
		int benchmark_scenario_id = -1;

		/// radius
		double radius_scale = 1;

		/// setup coordinate of e.g. radial breaking dam, x-placement \in [0;1]
		double setup_coord_x = 0.5;

		/// setup coordinate of e.g. radial breaking dam, y-placement \in [0;1]
		double setup_coord_y = 0.5;

		/// rotation angle for advection equation
		double advection_rotation_angle = 0;

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

		void output()
		{
			std::cout << std::endl;
			std::cout << "SETUP:" << std::endl;
			std::cout << " + benchmark_scenario_id: " << benchmark_scenario_id << std::endl;
			std::cout << " + radius_scale: " << radius_scale << std::endl;
			std::cout << " + setup_coord_x: " << setup_coord_x << std::endl;
			std::cout << " + setup_coord_y: " << setup_coord_y << std::endl;
			std::cout << " + advection_rotation_angle: " << advection_rotation_angle << std::endl;
			std::cout << " + input_data_filenames:" << std::endl;
			for (int i = 0; i < input_data_filenames.size(); i++)
				std::cout << "    - filename " << i << " " << input_data_filenames[i] << std::endl;
			std::cout << " + input_data_binary: " << input_data_binary << std::endl;
			std::cout << std::endl;
		}
	} setup;



	/**
	 * simulation coefficients
	 */
	struct SimulationCoefficients
	{
		/// average height for initialization
		double h0 = 1000.0;


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
		/// 7.292 × 10 −5
		double f0 = 0.00007292;

		/// Beta coefficient for f(y_N) = f0 + y_N*beta
		/// here, y_N is the normalized y coordinate \in [0;1]
		double beta = 0.0;

		// constants from Galwesky et al. paper

		/**
		 * Earth radius for simulations on the sphere
		 */
		double earth_radius = 6.37122e6;

		/**
		 * Coriolis effect
		 */
		double &coriolis_omega = f0;

		/**
		 * Gravitational constant
		 */
		double gravitation = 9.80616;

		/// zero out the v-component at the top and bottom layer
		bool top_bottom_zero_v_velocity = false;

		/// domain size
		double domain_size[2] = {1.0, 1.0};


		void output()
		{
			std::cout << std::endl;
			std::cout << "SIMULATION COEFFICIENTS:" << std::endl;
			std::cout << " + h0: " << h0 << std::endl;
			std::cout << " + viscosity: " << viscosity << std::endl;
			std::cout << " + viscosity_order: " << viscosity_order << std::endl;
			std::cout << " + CFL: " << CFL << std::endl;
			std::cout << " + f0: " << f0 << std::endl;
			std::cout << " + beta: " << beta << std::endl;
			std::cout << " + earth_radius: " << earth_radius << std::endl;
			std::cout << " + coriolis_omega: " << coriolis_omega << std::endl;
			std::cout << " + gravitation: " << gravitation << std::endl;
			std::cout << " + top_bottom_zero_v_velocity: " << top_bottom_zero_v_velocity << std::endl;
			std::cout << " + domain_size: " << domain_size[0] << " x " << domain_size[1] << std::endl;
			std::cout << std::endl;
		}
	} sim;


	/**
	 * This class stored the discretization-related parameters
	 *
	 * resolution / timestepping
	 */
	struct Discretization
	{
		/// resolution in physical space (grid cells)
		int res_physical[2] = {0, 0};

		/// resolution in spectral space (number of modes)
		int res_spectral[2] = {0, 0};

		/// size of cell (hx, hy)
		/// this is computed based on disc.res and sim.domain_size
		double cell_size[2] = {0,0};

		/// use up/downwinding for the advection of h
		bool timestepping_up_and_downwinding = false;

		/// order of Runge-Kutta scheme for time stepping
		double timestepping_runge_kutta_order = 4;

		/// use spectral differential operators
		bool use_spectral_basis_diffs =
#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_SPHERE_SPECTRAL_SPACE
				true;
#else
				false;
#endif

		void output()
		{
			std::cout << std::endl;
			std::cout << "DISCRETIZATION:" << std::endl;
			std::cout << " + res_physical: " << res_physical[0] << " x " << res_physical[1] << std::endl;
			std::cout << " + res_spectral: " << res_spectral[0] << " x " << res_spectral[1] << std::endl;
			std::cout << " + cell_size: " << res_physical[0] << " x " << cell_size[1] << std::endl;
			std::cout << " + timestepping_up_and_downwinding: " << timestepping_up_and_downwinding << std::endl;
			std::cout << " + timestepping_runge_kutta_order: " << timestepping_runge_kutta_order << std::endl;
			std::cout << " + use_spectral_basis_diffs: " << use_spectral_basis_diffs << std::endl;
			std::cout << std::endl;
		}
	} disc;



	/**
	 * REXI
	 */
	struct REXI
	{
		/**
		 * Activate REXI instead of standard time stepping
		 */
		int use_rexi = 0;

		/**
		 * REXI parameter h
		 */
		double rexi_h = 0.2;

		/**
		 * REXI parameter M
		 */
		int rexi_M = 128;

		/**
		 * REXI parameter L
		 */
		int rexi_L = 0;

		/**
		 * Use only half of the poles for REXI
		 */
		bool rexi_use_half_poles = true;

		/**
		 * Extend modes for certain operations
		 */
		int rexi_use_extended_modes = 0;

		/**
		 * Normalize REXI for geostrophic balance
		 */
		bool rexi_normalization = true;

		void output()
		{
			std::cout << std::endl;
			std::cout << "REXI:" << std::endl;
			std::cout << " + use_rexi: " << use_rexi << std::endl;
			std::cout << " + rexi_h: " << rexi_h << std::endl;
			std::cout << " + rexi_M: " << rexi_M << std::endl;
			std::cout << " + rexi_L: " << rexi_L << std::endl;
			std::cout << " + rexi_use_half_poles: " << rexi_use_half_poles << std::endl;
			std::cout << " + rexi_use_extended_modes: " << rexi_use_extended_modes << std::endl;
			std::cout << " + rexi_normalization: " << rexi_normalization << std::endl;
			std::cout << std::endl;
		}
	} rexi;


#if SWEET_PFASST_CPP
	struct Pfasst
	{
		int nlevels;
		int nnodes;
		int nspace;
		int nsteps;
		int niters;
		double dt;
#if 0
		// TODO: encapsulate this in sweet SimVars
		auto const nlevels   = pfasst::config::get_value<int>("nlevels", 1);
		auto const nnodes    = pfasst::config::get_value<int>("nnodes", 3);
		auto const nspace    = pfasst::config::get_value<int>("nspace", 8193);
		auto const nsteps    = pfasst::config::get_value<int>("nsteps", 16);
		auto const niters    = pfasst::config::get_value<int>("niters", 4);
		auto const dt        = pfasst::config::get_value<double>("dt", 0.1);
#endif
		void output()
		{
			std::cout << std::endl;
			std::cout << "PFASST:" << std::endl;
			std::cout << " [TODO]" << std::endl;
			std::cout << std::endl;
		}
	} pfasst;
#endif


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
		void output()
		{
			std::cout << std::endl;
			std::cout << "MISC:" << std::endl;
			std::cout << " + verbosity: " << verbosity << std::endl;
			std::cout << " + output_floating_point_precision: " << output_floating_point_precision << std::endl;
			std::cout << " + gui_enabled: " << gui_enabled << std::endl;
			std::cout << " + be_verbose_after_this_simulation_time_period: " << be_verbose_after_this_simulation_time_period << std::endl;
			std::cout << " + output_file_name_prefix: " << output_file_name_prefix << std::endl;
			std::cout << " + output_each_sim_seconds: " << output_each_sim_seconds << std::endl;
			std::cout << " + output_next_sim_seconds: " << output_next_sim_seconds << std::endl;
			std::cout << " + vis_id: " << vis_id << std::endl;
			std::cout << " + use_nonlinear_equations: " << use_nonlinear_equations << std::endl;
			std::cout << " + sphere_use_robert_functions: " << sphere_use_robert_functions << std::endl;
			std::cout << " + output_time_scale: " << output_time_scale << std::endl;
			std::cout << std::endl;
		}

		/// set verbosity of simulation
		int verbosity = 0;

		/// precision for floating point output to std::cout and std::endl
		int output_floating_point_precision = 18;

		/// activate GUI mode?
		bool gui_enabled = (SWEET_GUI == 0 ? false : true);

		/// output verbose information every given period of simulation time.
		double be_verbose_after_this_simulation_time_period = 0;

		/// prefix of filename for output of data
		std::string output_file_name_prefix = "prog_%s_t%020.8f.csv";

		/// prefix of filename for output of data
		double output_each_sim_seconds = -1;

		/// Simulation seconds for next output
		double output_next_sim_seconds = 0;

		/// id for visualization
		int vis_id = 0;

		/// Use non-linear equations for simulations
		int use_nonlinear_equations = 1;


		/// Use robert function formulation on the sphere
		bool sphere_use_robert_functions = true;

		/// time scaling for output
		/// e.g. use scaling by 1.0/(60*60) to output days instead of seconds
		double output_time_scale = 1.0;

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


		void output()
		{
			std::cout << std::endl;
			std::cout << "TIMECONTROL:" << std::endl;
			std::cout << " + run_simulation_timesteps: " << run_simulation_timesteps << std::endl;
			std::cout << " + current_timestep_nr: " << current_timestep_nr << std::endl;
			std::cout << " + current_timestep_size: " << current_timestep_size << std::endl;
			std::cout << " + current_simulation_time: " << current_simulation_time << std::endl;
			std::cout << " + max_timesteps_nr: " << max_timesteps_nr << std::endl;
			std::cout << " + max_simulation_time: " << max_simulation_time << std::endl;
			std::cout << std::endl;
		}
	} timecontrol;


	void output()
	{
		sim.output();
		disc.output();
		setup.output();
		timecontrol.output();

		rexi.output();
		misc.output();
		diag.output();

#if SWEET_PARAREAL
		parareal.output();
#endif

		pfasst.output();
	}


	/**
	 * update variables which are based on others
	 */
	void reset()
	{
		disc.cell_size[0] = sim.domain_size[0]/(double)disc.res_physical[0];
		disc.cell_size[1] = sim.domain_size[1]/(double)disc.res_physical[1];

		timecontrol.current_timestep_size = -1;
		timecontrol.current_timestep_nr = 0;
		timecontrol.current_simulation_time = 0;

		if ((disc.res_physical[0] & 1) || (disc.res_physical[1] & 1))
			std::cout << "WARNING: Typically there are only even resolutions supported!" << std::endl;
	}



	/**
	 * setup the variables based on program parameters
	 */
	bool setupFromMainParameters(
			int i_argc,					///< argc from main()
			char *const i_argv[],		///< argv from main()
			const char *bogus_var_names[] = nullptr			///< list of strings of simulation-specific variables, has to be terminated by nullptr
	)
	{
		const int max_options = 50;
        static struct option long_options[max_options+1] = {
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

				{0, 0, 0, 0}, // 0
				{0, 0, 0, 0}, // 1
				{0, 0, 0, 0}, // 2
				{0, 0, 0, 0}, // 3
				{0, 0, 0, 0}, // 4
				{0, 0, 0, 0}, // 5
				{0, 0, 0, 0}, // 6
				{0, 0, 0, 0}, // 7
				{0, 0, 0, 0}, // 8
				{0, 0, 0, 0}, // 9	Option Nr. 40

				{0, 0, 0, 0}, // 0
				{0, 0, 0, 0}, // 1
				{0, 0, 0, 0}, // 2
				{0, 0, 0, 0}, // 3
				{0, 0, 0, 0}, // 4
				{0, 0, 0, 0}, // 5
				{0, 0, 0, 0}, // 6
				{0, 0, 0, 0}, // 7
				{0, 0, 0, 0}, // 8
				{0, 0, 0, 0}, // 9	Option Nr. 50

				{0, 0, 0, 0} // NULL
        };

		int next_free_program_option = 0;

        long_options[next_free_program_option] = {"initial-coord-x", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"initial-coord-y", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"rexi-h", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"rexi-m", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"rexi-l", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"rexi-half", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"rexi", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"rexi-normalization", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"rexi-ext-modes", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"nonlinear", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"use-robert-functions", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"sphere-advection-angle", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

#if SWEET_PFASST_CPP
        long_options[next_free_program_option] = {"pfasst-nlevels", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"pfasst-nnodes", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"pfasst-nspace", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"pfasst-nsteps", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"pfasst-niters", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"pfasst-dt", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;
#endif


#if SWEET_PARAREAL
        int parareal_start_option_index = next_free_program_option;
        parareal.setup_longOptionList(long_options, next_free_program_option, max_options);
#endif

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
				long_options[opt_nr].val = 256+opt_nr;
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
			opt = getopt_long(
							i_argc, i_argv,
							"N:M:n:m:C:u:U:s:X:Y:f:b:x:y:t:i:T:v:V:O:o:H:r:a:R:W:F:S:g:G:d:z",
							long_options, &option_index
					);

			if (opt == -1)
				break;


			/*
			 * LONG OPTIONS
			 */
			if (opt >= 256)
			{
				int i = opt-256;

				if (i < next_free_program_option)
				{
					switch(i)
					{
						case 0:		setup.setup_coord_x = atof(optarg);	break;
						case 1:		setup.setup_coord_y = atof(optarg);	break;

						case 2:		rexi.rexi_h = atof(optarg);	break;
						case 3:		rexi.rexi_M = atoi(optarg);	break;
						case 4:		rexi.rexi_L = atoi(optarg);	break;
						case 5:		rexi.rexi_use_half_poles = atoi(optarg);	break;
						case 6:		rexi.use_rexi = atoi(optarg);	break;
						case 7:		rexi.rexi_normalization = atoi(optarg);	break;
						case 8:		rexi.rexi_use_extended_modes = atoi(optarg);	break;

						case 9:		misc.use_nonlinear_equations = atoi(optarg);	break;
						case 10:	misc.sphere_use_robert_functions = atoi(optarg);	break;

						case 11:	setup.advection_rotation_angle = atof(optarg);	break;

						case 12:	pfasst.nlevels = atoi(optarg);	break;
						case 13:	pfasst.nnodes = atoi(optarg);	break;
						case 14:	pfasst.nspace = atoi(optarg);	break;
						case 15:	pfasst.nsteps = atoi(optarg);	break;
						case 16:	pfasst.niters = atoi(optarg);	break;
						case 17:	pfasst.dt = atof(optarg);	break;

						default:
#if SWEET_PARAREAL
						parareal.setup_longOptionValue(i-parareal_start_option_index, optarg);
#endif
						break;
					}
				}
				else
				{
					int bogus_id = i-next_free_program_option;

					if (bogus_id >= max_options)
					{
						std::cout << std::endl;
						std::cout << "SERIOUS ERROR" << std::endl;
						std::cout << " + long option: " << i_argv[option_index] << std::endl;
						std::cout << " + bogus id " << bogus_id << std::endl;
						std::cout << std::endl;
						exit(1);
					}
					bogus.var[i-next_free_program_option] = atof(optarg);
				}
				continue;
			}


			if (optarg != nullptr)
			{
				if (optarg[0] == '=')
				{
					std::cerr << "Short option parameters may not be specified with an equal '=' sign!" << std::endl;
					exit(-1);
				}
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
				disc.res_physical[0] = atoi(optarg);
				disc.res_physical[1] = disc.res_physical[0];
				break;

			case 'M':
				disc.res_spectral[0] = atoi(optarg);
				disc.res_spectral[1] = disc.res_spectral[0];
				break;

			case 'n':
				disc.res_physical[0] = atoi(optarg);
				break;

			case 'm':
				disc.res_physical[1] = atoi(optarg);
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

			case 's':
				setup.benchmark_scenario_id = atoi(optarg);
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
				setup.setup_coord_x = atof(optarg);
				break;

			case 'y':
				setup.setup_coord_y = atof(optarg);
				break;

			case 'f':
				sim.f0 = atof(optarg);
				break;

			case 'a':
				sim.earth_radius = atof(optarg);
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
				sim.gravitation = atof(optarg);
				break;

			case 'v':
				misc.verbosity = atoi(optarg);
				break;

			case 'V':
				misc.be_verbose_after_this_simulation_time_period = atof(optarg);
				break;

			case 'O':
				misc.output_file_name_prefix = optarg;
				if (misc.output_file_name_prefix == "-")
					misc.output_file_name_prefix = "";
				break;

			case 'o':
				misc.output_each_sim_seconds = atof(optarg);
				break;

			case 'H':
				sim.h0 = atof(optarg);
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


			default:
				const char *help_strings[] = {
						"Simulation runtime parameters",
						"	-X [length]	length of simulation domain in x direction, default=1",
						"	-Y [width]	width of simulation domain in y direction, default=1",
						"	-u [visc]	viscosity, , default=0",
						"	-U [visc]	viscosity order, default=2",
						//"	-p [visc]	potential viscosity, default=0",
						//"	-P [visc]	potential hyperviscosity, default=0",
						"	-f [float]	f-parameter for f-plane or coriolis omega term, default=0",
						"	-b [float]	beta-parameter for beta-plane, default=0",
						"	            Use -1 to set f*sin(phi) with phi in [-pi/2;pi/2] in y",
						"	-g [float]	gravity, default=9.81",
						"	-a [float]	earth radius, default=1",
						"	-H [float]	average (initial) height of water, default=1000",
						"",
						"Simulation setup parameters",
						"	-s [scen]	scenario id, set to -1 for overview"
						"	-x [float]	x coordinate for setup \\in [0;1], default=0.5",
						"	-y [float]	y coordinate for setup \\in [0;1], default=0.5",
						"	-r [radius]	scale factor of radius for initial condition, default=1",
						"	--initial-freq-x-mul [float]	Frequency for the waves initial conditions in x, default=2",
						"	--initial-freq-y-mul [float]	Frequency for the waves initial conditions in y, default=1",
						"	--initial-coord-x [float]	Same as -x",
						"	--initial-coord-y [float]	Same as -y",
						"",
						"Discretization:",
						"  >Space:",
						"	-N [res]	resolution in x and y direction, default=0",
						"	-n [resx]	resolution in x direction, default=0",
						"	-m [resy]	resolution in y direction, default=0",
						"	-M [modes]	modes in x/y or lon/lat direction, default=0",
						"	-S [0/1]	Control Operator discretization for PlaneDatas",
						"               0: FD, 1: spectral derivatives, default:"
#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_SPHERE_SPECTRAL_SPACE
				"1"
#else
				"0"
#endif
						,
						"  >Time:",
						"	-W [0/1]	use up- and downwinding, default:0",
						"	-R [1-RKn]	order of Runge-Kutta method, default:4",
						"	-C [cfl]	CFL condition, use negative value for fixed time step size, default=0.05",
						"",
						"Control:",
						"	-t [time]	maximum simulation time, default=-1 (infinity)",
						"	-T [stepnr]	maximum number of time steps, default=-1 (infinity)",
						"	-o [time]	time interval at which output should be written",
						"",
						"Misc options",
						"	-v [int]	verbosity level",
						"	-V [double]	period of output",
						"	-G [0/1]	graphical user interface",
						"	-O [string]	string prefix for filename of output of simulation data",
						"	-d [int]	accuracy of floating point output",
						"	-i [file0][;file1][;file3]...	string with filenames for initial conditions",
						"	            specify BINARY; as first file name to read files as binary raw data",
						"	--use-robert-functions [bool]	Use Robert function formulation for velocities on the sphere",
						"	--nonlinear [bool]	Use non-linear (1) if available or linear (0) formulation, default: 1",
						"",
						"Rexi",
						"	--rexi [bool]	Time stepping method: 0: explicit, 1: REXI, -1: implicit",
						"	--rexi-h [float]	REXI parameter h",
						"	--rexi-m [int]	REXI parameter M",
						"	--rexi-l [int]	REXI parameter L",
						"	--rexi-half [bool]	Use half REXI poles, default:1",
						"	--rexi-normalization [bool]	Use REXI normalization around geostrophic balance, default:1",
						"	--rexi-ext-modes [int]	Use this number of extended modes in spherical harmonics",
						"",

				};

				std::cerr << "Usage information: " << std::endl;
				for (std::size_t i = 0; i < sizeof(help_strings)/sizeof(*help_strings); i++)
					std::cerr << help_strings[i] << std::endl;


#if SWEET_PARAREAL
				parareal.setup_printOptions();
#endif

				std::cerr << std::endl;
				std::cerr << "Unknown option '" << (char)opt << "'" << std::endl;
				return false;
			}
		}

		if (	(disc.res_physical[0] == 0 || disc.res_physical[1] == 0)	&&
				(disc.res_spectral[0] == 0 || disc.res_spectral[1] == 0)
			)
		{
			FatalError("Select physical resolution or spectral modes");
		}

		reset();



#if SWEET_PARAREAL
		// if max simulation time was not set for parareal, copy max simulation time from default parameters to parareal parameters.
		if (parareal.max_simulation_time <= 0)
			parareal.max_simulation_time = timecontrol.max_simulation_time;
#endif

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

		if (misc.output_file_name_prefix == "-")
			misc.output_file_name_prefix = "";

		std::cout << std::setprecision(misc.output_floating_point_precision);
		std::cerr << std::setprecision(misc.output_floating_point_precision);

		return true;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
