/*
 * SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
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
#include <sweet/StringSplit.hpp>


#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#ifndef SWEET_GUI
#	define SWEET_GUI 1
#endif

#ifndef SWEET_PFASST_CPP
#	define SWEET_PFASST_CPP 1
#endif

#ifndef SWEET_LIBPFASST
#	define SWEET_LIBPFASST 1
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal_SimulationVariables.hpp>
#endif

#if SWEET_LIBPFASST
#       include <libpfasst/LibPFASST_SimulationVariables.hpp>
#endif

#include <rexi/REXI_SimulationVariables.hpp>

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

public:
	REXI_SimulationVariables rexi;


public:
	/**
	 * Diagnostic variables
	 */
	struct Diagnostics
	{
		/// total mass
		double total_mass = 0;

		/// kinetic energy
		double kinetic_energy = 0;

		/// potential energy
		double potential_energy = 0;

		/// total energy
		double total_energy = 0;

		/// total potential enstropy
		double total_potential_enstrophy = 0;


		double ref_total_mass = -1;
		double ref_kinetic_energy = -1;
		double ref_potential_energy = -1;
		double ref_total_energy = -1;
		double ref_total_potential_enstrophy = -1;


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "DIAGNOSTICS:" << std::endl;
			std::cout << " + total_mass: " << total_mass << std::endl;
			std::cout << " + total_energy: " << total_energy << std::endl;
			std::cout << " + kinetic_energy: " << kinetic_energy << std::endl;
			std::cout << " + potential_energy: " << potential_energy << std::endl;
			std::cout << " + total_potential_enstrophy: " << total_potential_enstrophy << std::endl;
			std::cout << std::endl;
		}

		void backup_reference()
		{
			ref_total_mass = total_mass;
			ref_kinetic_energy = kinetic_energy;
			ref_potential_energy = potential_energy;
			ref_total_energy = total_energy;
			ref_total_potential_enstrophy = total_potential_enstrophy;
		}
	} diag;



public:
	/**
	 * Information on Partial Differential Equation
	 */
	struct PDE
	{
		/*
		 * ID of PDE to use
		 */
		int id = 0;

		/// Use non-linear equations for simulations
		int use_nonlinear_equations = 1;


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "PDE:" << std::endl;
			std::cout << " + id: " << id << std::endl;
			std::cout << " + use_nonlinear_equations: " << use_nonlinear_equations << std::endl;
			std::cout << std::endl;
		}


		void outputProgParams()
		{
			std::cout << std::endl;
			std::cout << "Partial differential equation:" << std::endl;
			std::cout << "	--pde-id [0/1]		PDE to solve (0: SWE, 1: advection)" << std::endl;
			std::cout << "	--nonlinear [int]		Use non-linear (>=1) if available or linear (0) formulation, default: 1" << std::endl;
			std::cout << "						0: Linear " << std::endl;
			std::cout << "						1: Nonlinear (default)" << std::endl;
			std::cout << "						2: Linear + nonlinear advection only (needs -H to be set)" << std::endl;
			std::cout << "" << std::endl;
		}
	} pde;



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


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "SETUP:" << std::endl;
			std::cout << " + benchmark_scenario_id: " << benchmark_scenario_id << std::endl;
			std::cout << " + radius_scale: " << radius_scale << std::endl;
			std::cout << " + setup_coord_x: " << setup_coord_x << std::endl;
			std::cout << " + setup_coord_y: " << setup_coord_y << std::endl;
			std::cout << " + advection_rotation_angle: " << advection_rotation_angle << std::endl;
			std::cout << " + input_data_filenames:" << std::endl;
			for (std::size_t i = 0; i < input_data_filenames.size(); i++)
				std::cout << "    - filename " << i << " " << input_data_filenames[i] << std::endl;
			std::cout << " + input_data_binary: " << input_data_binary << std::endl;
			std::cout << std::endl;
		}


		void outputProgParams()
		{
			std::cout << std::endl;
			std::cout << "SIMULATION SETUP PARAMETERS:" << std::endl;
			std::cout << "	-s [scen]				scenario id, set to -1 for overview" << std::endl;
			std::cout << "	-x [float]				x coordinate for setup \\in [0;1], default=0.5" << std::endl;
			std::cout << "	-y [float]				y coordinate for setup \\in [0;1], default=0.5" << std::endl;
			std::cout << "	-r [radius]				scale factor of radius for initial condition, default=1" << std::endl;
			std::cout << "	--initial-freq-x-mul [float]		Frequency for the waves initial conditions in x, default=2" << std::endl;
			std::cout << "	--initial-freq-y-mul [float]		Frequency for the waves initial conditions in y, default=1" << std::endl;
			std::cout << "	--initial-coord-x [float]		Same as -x" << std::endl;
			std::cout << "	--initial-coord-y [float]		Same as -y" << std::endl;
			std::cout << "	--advection-rotation-angle [float]	Rotation angle for e.g. advection test case" << std::endl;

			std::cout << "" << std::endl;
		}
	} setup;



	/**
	 * simulation coefficients
	 */
	struct SimulationCoefficients
	{
		/// average height for initialization
		/// use value similar to Galewski benchmark
		double h0 = 10000.0;


		/// For more information on viscosity,
		/// see 13.3.1 "Generic Form of the Explicit Diffusion Mechanism"
		/// in "Numerical Techniques for Global Atmospheric Models"

		/// viscosity-term on velocities with 2nd order diff operator
		double viscosity = 0.0;

		/// hyper viscosity-term on velocities with 4th order diff operator
		int viscosity_order = 2;

		/// CFL condition
		double CFL = 0.05;


		/**
		 * Coriolis effect
		 * 7.2921 × 10 −5
		 */
		double coriolis_omega = 0.000072921;


		/**
		 * Coriolis frequency f0
		 */
		double f0 = 0.00007292*2;

		// constants from Galwesky et al. paper

		/**
		 * Earth radius for simulations on the sphere
		 */
		double earth_radius = 6.37122e6;


		/**
		 * Simulation on f-sphere? (constant f0 term over entire sphere)
		 */
		bool f_sphere = false;

		/**
		 * Gravitational constant
		 */
		double gravitation = 9.80616;

		/// zero out the v-component at the top and bottom layer
		bool top_bottom_zero_v_velocity = false;

		/// domain size
		double domain_size[2] = {1.0, 1.0};


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "SIMULATION COEFFICIENTS:" << std::endl;
			std::cout << " + h0: " << h0 << std::endl;
			std::cout << " + viscosity: " << viscosity << std::endl;
			std::cout << " + viscosity_order: " << viscosity_order << std::endl;
			std::cout << " + CFL: " << CFL << std::endl;
			std::cout << " + f0: " << f0 << std::endl;
			std::cout << " + earth_radius: " << earth_radius << std::endl;
			std::cout << " + coriolis_omega: " << coriolis_omega << std::endl;
			std::cout << " + f_sphere: " << f_sphere << std::endl;
			std::cout << " + gravitation: " << gravitation << std::endl;
			std::cout << " + top_bottom_zero_v_velocity: " << top_bottom_zero_v_velocity << std::endl;
			std::cout << " + domain_size (2D): " << domain_size[0] << " x " << domain_size[1] << std::endl;
			std::cout << std::endl;
		}


		void outputProgParams()
		{
			std::cout << "Simulation parameters" << std::endl;
			std::cout << "	-X [length]	length of simulation domain in x direction, default=1" << std::endl;
			std::cout << "	-Y [width]	width of simulation domain in y direction, default=1" << std::endl;
			std::cout << "	-u [visc]	viscosity, , default=0" << std::endl;
			std::cout << "	-U [visc]	viscosity order, default=2" << std::endl;
			std::cout << "	-f [float]	f-parameter for f-plane or coriolis omega term, default=0" << std::endl;
			std::cout << "	-F [int]	Simulation on f-sphere, default=0" << std::endl;
			std::cout << "	-g [float]	gravity" << std::endl;
			std::cout << "	-a [float]	earth radius" << std::endl;
			std::cout << "	-H [float]	average (initial) height of water" << std::endl;
			std::cout << "" << std::endl;

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
//		double cell_size[2] = {0, 0};



		/// Leapfrog: Robert Asselin filter
		double leapfrog_robert_asselin_filter = 0;

		/// Crank-Nicolson filter
		double crank_nicolson_filter = 0.5;


		/// String of time stepping method
		/// See doc/swe/swe_plane_timesteppings
		std::string timestepping_method;

		/// Order of time stepping
		int timestepping_order = -1;

		/// Order of 2nd time stepping which might be used
		int timestepping_order2 = 0;


		/// use spectral differential operators
		bool use_spectral_basis_diffs =
#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_SPHERE_SPECTRAL_SPACE
				true;
#else
				false;
#endif

		bool use_staggering = false;

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		int normal_mode_analysis_generation = 0;


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "DISCRETIZATION:" << std::endl;
			std::cout << " + res_physical: " << res_physical[0] << " x " << res_physical[1] << std::endl;
			std::cout << " + res_spectral: " << res_spectral[0] << " x " << res_spectral[1] << std::endl;
//			std::cout << " + cell_size (2D): " << res_physical[0] << " x " << cell_size[1] << std::endl;
			std::cout << " + timestepping_method: " << timestepping_method << std::endl;
			std::cout << " + timestepping_order: " << timestepping_order << std::endl;
			std::cout << " + timestepping_order2: " << timestepping_order2 << std::endl;
			std::cout << " + leapfrog_robert_asselin_filter: " << leapfrog_robert_asselin_filter << std::endl;
			std::cout << " + crank_nicolson_filter: " << crank_nicolson_filter << std::endl;
			std::cout << " + use_spectral_basis_diffs: " << use_spectral_basis_diffs << std::endl;
			std::cout << " + use_staggering: " << use_staggering << std::endl;
			std::cout << " + normal_mode_analysis_generation: " << normal_mode_analysis_generation << std::endl;

			std::cout << " + dealiasing (compile time): " <<
		#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
					1
		#else
					0
		#endif
					<< std::endl;
			std::cout << std::endl;
		}

		void outputProgParams()
		{
			std::cout << "Discretization:" << std::endl;
			std::cout << "  >Space:" << std::endl;
			std::cout << "	--staggering [0/1]	Use staggering" << std::endl;
			std::cout << "	-N [res]		resolution in x and y direction, default=0" << std::endl;
			std::cout << "	-n [resx]		resolution in x direction, default=0" << std::endl;
			std::cout << "	-m [resy]		resolution in y direction, default=0" << std::endl;
			std::cout << "	-M [modes]		modes in x/y or lon/lat direction, default=0" << std::endl;
			std::cout << "	-S [0/1]		Control Operator discretization for PlaneData" << std::endl;
			std::cout << "					0: FD, 1: spectral derivatives, default: ";

#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_SPHERE_SPECTRAL_SPACE
	std::cout << "1" << std::endl;
#else
	std::cout << "0" << std::endl;
#endif

			std::cout << "  >Time:" << std::endl;
			std::cout << "	-W [0/1]					use up- and downwinding, default:0" << std::endl;
			std::cout << "	-R [1-RKn]					order of time stepping method, default:0" << std::endl;
			std::cout << "	-C [cfl]					CFL condition, use negative value for fixed time step size, default=0.05" << std::endl;
			std::cout << "	--timestepping-method [string]	String of time stepping method" << std::endl;
			std::cout << "	--timestepping-order [int]			Specify the order of the time stepping" << std::endl;
			std::cout << "	--timestepping-order2 [int]			Specify the order of the time stepping" << std::endl;
			std::cout << "	--leapfrog-robert-asselin-filter [0;1]		Damping parameter for Robert-Asselin filter" << std::endl;
			std::cout << "	--normal-mode-analysis-generation [0;1;2;3]	Generate output data for normal mode analysis" << std::endl;
			std::cout << "							0: don't generate" << std::endl;
			std::cout << "							1: generate in physical space" << std::endl;
			std::cout << "							2: generate in spectral space" << std::endl;
			std::cout << "							3: generate in spectral space with complex matrix)" << std::endl;

		}

	} disc;



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
		auto const nlevels   = pfasst::config::get_value<int>("nlevels" << std::endl; 1);
		auto const nnodes    = pfasst::config::get_value<int>("nnodes" << std::endl; 3);
		auto const nspace    = pfasst::config::get_value<int>("nspace" << std::endl; 8193);
		auto const nsteps    = pfasst::config::get_value<int>("nsteps" << std::endl; 16);
		auto const niters    = pfasst::config::get_value<int>("niters" << std::endl; 4);
		auto const dt        = pfasst::config::get_value<double>("dt" << std::endl; 0.1);
#endif
		void outputConfig()
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
		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "MISC:" << std::endl;
			std::cout << " + verbosity: " << verbosity << std::endl;
			std::cout << " + compute_errors " << compute_errors << std::endl;
			std::cout << " + stability_checks: " << stability_checks << std::endl;
			std::cout << " + output_floating_point_precision: " << output_floating_point_precision << std::endl;
			std::cout << " + gui_enabled: " << gui_enabled << std::endl;
			std::cout << " + be_verbose_after_this_simulation_time_period: " << be_verbose_after_this_simulation_time_period << std::endl;
			std::cout << " + output_file_name_prefix: " << output_file_name_prefix << std::endl;
			std::cout << " + output_each_sim_seconds: " << output_each_sim_seconds << std::endl;
			std::cout << " + output_next_sim_seconds: " << output_next_sim_seconds << std::endl;
			std::cout << " + vis_id: " << vis_id << std::endl;
			std::cout << " + sphere_use_robert_functions: " << sphere_use_robert_functions << std::endl;
			std::cout << " + output_time_scale: " << output_time_scale << std::endl;
			std::cout << std::endl;
		}


		/// set verbosity of simulation
		int verbosity = 0;

		/// compute errors
		int compute_errors = 0;

		/// do stability checks for simulation
		int stability_checks = 1;

		/// precision for floating point outputConfig to std::cout and std::endl
		int output_floating_point_precision = -1;

		/// activate GUI mode?
		bool gui_enabled = (SWEET_GUI == 0 ? false : true);

		/// outputConfig verbose information every given period of simulation time.
		double be_verbose_after_this_simulation_time_period = 0;

		/// prefix of filename for outputConfig of data
		std::string output_file_name_prefix = "output_%s_t%020.8f.csv";

		/// prefix of filename for outputConfig of data
		double output_each_sim_seconds = -1;

		/// Simulation seconds for next outputConfig
		double output_next_sim_seconds = 0;

		/// id for visualization
		int vis_id = 0;


		/// Use robert function formulation on the sphere
		bool sphere_use_robert_functions = true;

		/// time scaling for outputConfig
		/// e.g. use scaling by 1.0/(60*60) to output days instead of seconds
		double output_time_scale = 1.0;
	} misc;


	/**
	 * Timestepping
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
		int max_timesteps_nr = std::numeric_limits<int>::max();

		/// maximum simulation time to execute the simulation for
		double max_simulation_time = std::numeric_limits<double>::infinity();


		void outputConfig()
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


	void outputConfig()
	{
		sim.outputConfig();
		disc.outputConfig();
		setup.outputConfig();
		timecontrol.outputConfig();

		rexi.outputConfig();
		pde.outputConfig();
		misc.outputConfig();
		diag.outputConfig();

#if SWEET_PARAREAL
		parareal.outputConfig();
#endif

#if SWEET_LIBPFASST
		libpfasst.outputConfig();
#endif


#if SWEET_PFASST_CPP 
		pfasst.outputConfig();
#endif
	}


	/**
	 * update variables which are based on others
	 */
	void reset()
	{
		if (timecontrol.max_simulation_time < 0)
			FatalError("timecontrol.max_simulation_time < 0");

		if (timecontrol.max_timesteps_nr < 0)
			FatalError("timecontrol.max_timesteps_nr < 0");

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
			const char *bogus_var_names[] = nullptr,			///< list of strings of simulation-specific variables, has to be terminated by nullptr
			bool i_run_prog_parameter_validation = true
	)
	{
		const int max_options = 100;
        struct option long_options[max_options+1];

        for (std::size_t i = 0; i < max_options+1; i++)
        {
        	long_options[i].flag = 0;
        	long_options[i].has_arg = 0;
        	long_options[i].name = 0;
        	long_options[i].val = 0;
        }

		int next_free_program_option = 0;

		// SETUP
        long_options[next_free_program_option] = {"initial-coord-x", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"initial-coord-y", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"advection-rotation-angle", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;


        // MISC
        long_options[next_free_program_option] = {"compute-errors", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"stability-checks", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"use-robert-functions", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;


        // PDE
        long_options[next_free_program_option] = {"nonlinear", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"pde-id", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        // DISC
        long_options[next_free_program_option] = {"timestepping-method", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"timestepping-order", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"timestepping-order2", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"leapfrog-robert-asselin-filter", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"normal-mode-analysis-generation", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"crank-nicolson-filter", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"staggering", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

        long_options[next_free_program_option] = {"dt", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;



// leave this commented to avoid mismatch with following parameters!
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


        int rexi_start_option_index = next_free_program_option;
        rexi.setup_longOptionList(
        		long_options,
				next_free_program_option,	///< also updated (IO)
				max_options
			);

        // Test dummy object
        long_options[next_free_program_option] = {"dummy", required_argument, 0, 256+next_free_program_option};
        next_free_program_option++;

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
				FatalError("Max number of arguments reached. Reduce number of program arguments");
			}
        }

		// index into long_options for determined argument
		int option_index = 0;

		int opt;
		while (1)
		{
			opt = getopt_long(
							i_argc, i_argv,
							"N:M:n:m:C:u:U:s:X:Y:f:F:b:x:y:t:i:T:v:V:O:o:H:r:a:R:W:F:S:g:G:d:z",
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
					int c = 0;	if (i == c)	{	setup.setup_coord_x = atof(optarg);		continue;	}
					c++;		if (i == c)	{	setup.setup_coord_y = atof(optarg);		continue;	}
					c++;		if (i == c)	{	setup.advection_rotation_angle = atof(optarg);		continue;	}

					c++;		if (i == c)	{	misc.compute_errors = atoi(optarg);					continue;	}
					c++;		if (i == c)	{	misc.stability_checks = atoi(optarg);				continue;	}
					c++;		if (i == c)	{	misc.sphere_use_robert_functions = atoi(optarg);	continue;	}

					c++;		if (i == c)	{	pde.use_nonlinear_equations = atoi(optarg);			continue;	}
					c++;		if (i == c)	{	pde.id = atoi(optarg);								continue;	}

					c++;		if (i == c)	{	disc.timestepping_method = optarg;					continue;	}
					c++;		if (i == c)	{	disc.timestepping_order = atoi(optarg);				continue;	}
					c++;		if (i == c)	{	disc.timestepping_order2 = atoi(optarg);			continue;	}
					c++;		if (i == c)	{	disc.leapfrog_robert_asselin_filter = atof(optarg);	continue;	}
					c++;		if (i == c)	{	disc.normal_mode_analysis_generation = atoi(optarg);	continue;	}
					c++;		if (i == c)	{	disc.crank_nicolson_filter = atof(optarg);			continue;	}
					c++;		if (i == c)	{	disc.use_staggering = atof(optarg);					continue;	}

					c++;		if (i == c)	{	timecontrol.current_timestep_size = atof(optarg);		continue;	}

#if SWEET_PFASST_CPP 
					c++;		if (i == c)	{	pfasst.nlevels = atoi(optarg);	continue;	}
					c++;		if (i == c)	{	pfasst.nnodes = atoi(optarg);	continue;	}
					c++;		if (i == c)	{	pfasst.nspace = atoi(optarg);	continue;	}
					c++;		if (i == c)	{	pfasst.nsteps = atoi(optarg);	continue;	}
					c++;		if (i == c)	{	pfasst.niters = atoi(optarg);	continue;	}
					c++;		if (i == c)	{	pfasst.dt = atof(optarg);	continue;	}
#endif

#if SWEET_PARAREAL
					{
						int retval = parareal.setup_longOptionValue(i-parareal_start_option_index, optarg);
						if (retval == 0)
							continue;
						c += retval;
					}
#endif

#if SWEET_LIBPFASST
					{
					  int retval = libpfasst.setup_longOptionValue(i-libpfasst_start_option_index, optarg);
					  if (retval == 0)
					    continue;
					  c += retval;
					}
#endif

					{
						int retval = rexi.setup_longOptionValue(i-rexi_start_option_index, optarg);
						if (retval == 0)
							continue;
						c += retval;
					}

					c++;

					/*
					 * This can be tested with the --dummy parameter
					 */
					if (c != next_free_program_option-1)
						FatalError("Inconsistent processing of arguments");

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
				{
					std::vector<std::string> res = StringSplit::split(optarg, ",");

					int c = res.size();

					if (c == 0)
						FatalError("Invaild format for modes");

					if (c == 1)
					{
						disc.res_physical[0] = atoi(res[0].c_str());
						disc.res_physical[1] = disc.res_physical[0];
					}
					else if (c == 2)
					{
						disc.res_physical[0] = atoi(res[0].c_str());
						disc.res_physical[1] = atoi(res[1].c_str());
					}
					else
					{
						FatalError("More than 2 max resolutions given");
					}
				}
				break;

			case 'M':
				{
					std::vector<std::string> modes = StringSplit::split(optarg, ",");

					int c = modes.size();

					if (c == 0)
						FatalError("Invaild format for modes");

					if (c == 1)
					{
						disc.res_spectral[0] = atoi(modes[0].c_str());
						disc.res_spectral[1] = disc.res_spectral[0];
					}
					else if (c == 2)
					{
						disc.res_spectral[0] = atoi(modes[0].c_str());
						disc.res_spectral[1] = atoi(modes[1].c_str());
					}
					else
					{
						FatalError("More than 2 max modes given");
					}
				}
				break;

			case 'n':
				disc.res_physical[0] = atoi(optarg);
				break;

			case 'm':
				disc.res_physical[1] = atoi(optarg);
				break;

			case 'C':
				sim.CFL = atof(optarg);
				if (sim.CFL < 0)
					timecontrol.current_timestep_size = -sim.CFL;

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
				sim.coriolis_omega = atof(optarg);
				break;

			case 'F':
				sim.f_sphere = atoi(optarg);
				break;

			case 'a':
				sim.earth_radius = atof(optarg);
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
				disc.timestepping_order = atoi(optarg);
				break;

			case 'i':
				setup.setup_initial_condition_filenames(optarg);
				break;


			default:
				sim.outputProgParams();
				setup.outputProgParams();
				pde.outputProgParams();
				disc.outputProgParams();

				std::cout << "" << std::endl;
				std::cout << "Control:" << std::endl;
				std::cout << "	-t [time]	maximum simulation time, default=-1 (infinity)" << std::endl;
				std::cout << "	-T [stepnr]	maximum number of time steps, default=-1 (infinity)" << std::endl;
				std::cout << "	-o [time]	time interval at which output should be written, (set to 0 for output at every time step), default=-1 (no output) " << std::endl;
				std::cout << "" << std::endl;
				std::cout << "Misc options:" << std::endl;
				std::cout << "	-v [int]			verbosity level" << std::endl;
				std::cout << "	-V [double]			period of outputConfig" << std::endl;
				std::cout << "	-G [0/1]			graphical user interface" << std::endl;
				std::cout << "	-O [string]			string prefix for filename of output of simulation data" << std::endl;
				std::cout << "	-d [int]			accuracy of floating point output" << std::endl;
				std::cout << "	-i [file0][;file1][;file3]...	string with filenames for initial conditions" << std::endl;
				std::cout << "					specify BINARY; as first file name to read files as binary raw data" << std::endl;
				std::cout << "	--compute-errors [int]          Compute errors when possible [1], default=0	" << std::endl;
				std::cout << "	--use-robert-functions [bool]	Use Robert function formulation for velocities on the sphere" << std::endl;
				std::cout << "" << std::endl;
				rexi.outputProgParams();


#if SWEET_PARAREAL
				parareal.printOptions();
#endif

#if SWEET_LIBPFASST
				libpfasst.printOptions();
#endif

				std::cerr << std::endl;

				if ((char)opt != 'h')
					std::cerr << "Unknown option '" << (char)opt << "'" << std::endl;
				return false;
			}
		}


		if (i_run_prog_parameter_validation)
		{
			if (	(disc.res_physical[0] == 0 || disc.res_physical[1] == 0)	&&
					(disc.res_spectral[0] == 0 || disc.res_spectral[1] == 0)
				)
			{
				FatalError("Select physical resolution or spectral modes");
			}
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

		if (misc.output_floating_point_precision >= 0)
		{
			std::cout << std::setprecision(misc.output_floating_point_precision);
			std::cerr << std::setprecision(misc.output_floating_point_precision);
		}

		return true;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
