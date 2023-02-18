/*
 * SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
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
#include <sweet/dict/Dict.hpp>

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

#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#ifndef SWEET_GUI
#	define SWEET_GUI 1
#endif

#ifndef SWEET_LIBPFASST
#	define SWEET_LIBPFASST 1
#endif

#ifndef SWEET_XBRAID
#	define SWEET_XVRAID 1
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


#define SWEET_USE_PLANE_TEST 1


/*
 * REXI program parameters
 */
#include <rexi/REXI_SimulationVariables.hpp>

/*
 * Polvani program parameters
 */
#include <../programs/swe_plane_benchmarks/SWE_bench_Polvani_SimulationVariables.hpp>


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


	/**
	 * Input and output data
	 */
	struct IOData
	{
		/// filenames of input data for setup (this has to be setup by each application individually)
		std::vector<std::string> initial_condition_data_filenames;

		/// use "BINARY;filename1;filename2" to specify that the binary files should be read in binary format
		bool initial_condition_input_data_binary = false;


		/// prefix of filename for outputConfig of data
		std::string output_file_name = "X";

		/// output mode of variables
		std::string output_file_mode = "default";

		/// prefix of filename for outputConfig of data
		double output_each_sim_seconds = -1;

		/// Simulation seconds for next outputConfig
		double output_next_sim_seconds = 0;

		/// time scaling for outputConfig
		/// e.g. use scaling by 1.0/(60*60) to output days instead of seconds
		double output_time_scale = 1.0;

		/// precision for floating point outputConfig to std::cout and std::endl
		int output_floating_point_precision = std::numeric_limits<double>::digits10 + 1;



		void setup_initial_condition_filenames(
				const std::string i_string
		)
		{
			std::size_t last_pos = 0;
			for (std::size_t pos = 0; i_string[pos] != '\0'; pos++)
			{
				if (i_string[pos] != ';')
					continue;

				initial_condition_data_filenames.push_back(i_string.substr(last_pos, pos-last_pos));
				last_pos = pos+1;
			}

			initial_condition_data_filenames.push_back(i_string.substr(last_pos));

			if (initial_condition_data_filenames.size() > 0)
			{
				if (initial_condition_data_filenames[0] == "BINARY")
				{
					initial_condition_data_filenames.erase(initial_condition_data_filenames.begin());
					initial_condition_input_data_binary = true;
				}
			}
		}

		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "INPUT/OUTPUT:" << std::endl;
			for (std::size_t i = 0; i < initial_condition_data_filenames.size(); i++)
				std::cout << "    - filename " << i << " " << initial_condition_data_filenames[i] << std::endl;
			std::cout << " + input_data_binary: " << initial_condition_input_data_binary << std::endl;
			std::cout << " + output_file_name " << output_file_name << std::endl;
			std::cout << " + output_file_mode " << output_file_mode << std::endl;
			std::cout << " + output_each_sim_seconds: " << output_each_sim_seconds << std::endl;
			std::cout << " + output_next_sim_seconds: " << output_next_sim_seconds << std::endl;
			std::cout << " + output_time_scale: " << output_time_scale << std::endl;
			std::cout << " + output_floating_point_precision: " << output_floating_point_precision << std::endl;
			std::cout << std::endl;
		}



		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{
	        long_options[next_free_program_option] = {"output-file-name", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"output-file-mode", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;
		}

		void outputProgParams()
		{
			std::cout << std::endl;
			std::cout << "IOData:" << std::endl;
			std::cout << "	--output-file-name [string]		String specifying the name of the output file" << std::endl;
			std::cout << "	--output-file-mode [string]		Format of output file, default: default" << std::endl;

			std::cout << "" << std::endl;
		}



		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{
			switch(i_option_index)
			{
			case 0:
				output_file_name = i_value;
				return -1;

			case 1:
				output_file_mode = i_value;
				return -1;
			}

			return 2;
		}

	} iodata;



public:
	/**
	 * Values and parameters to setup benchmarks simulations
	 */
	struct Benchmark
	{
		/// seed for random number generator
		int random_seed = 0;

		/// benchmark scenario
		std::string benchmark_name = "";

		/// May the benchmark setup overwrite the simulation variables
		bool benchmark_override_simvars = true;

		/// Use 2/3 rule in physical space for dealiasing
		bool setup_dealiased = true;


		/// Galewsky-benchmark specific: velocity
		double benchmark_galewsky_umax = -1;

		/// Galewsky-benchmark specific: amplitude of bump
		double benchmark_galewsky_hamp = -1;

		/// Galewsky-benchmark specific: latitude coordinate
		double benchmark_galewsky_phi2 = -1;

		/// Normal modes benchmark scenario
		std::string benchmark_normal_modes_case = "";

		/// radius
		double object_scale = 1;

		/// setup coordinate of e.g. radial breaking dam, x-placement \in [0;1]
		double object_coord_x = 0.5;

		/// setup coordinate of e.g. radial breaking dam, y-placement \in [0;1]
		double object_coord_y = 0.5;

		/// rotation angle for advection equation
		double sphere_advection_rotation_angle = 0;


		/**
		 * Flag to indicate the presence of topography
		 */
		bool use_topography = false;


#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		/**
		 * Topography vector
		 */
		SphereData_Physical h_topo;
#endif



#if 0
		static void fun_no_forces(int, double, void*, void*)
		{
			SWEETError("External forces not available");
		};
#endif


		/// load external forces if available from benchmark scenario
		void (*getExternalForcesCallback)(int, double, void*, void*) = nullptr;// = &fun_no_forces;		/// SET TO NULLPTR
		void *getExternalForcesUserData = nullptr;



		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "BENCHMARK:" << std::endl;
			std::cout << " + random_seed: " << random_seed << std::endl;
			std::cout << " + benchmark_name: " << benchmark_name << std::endl;
			std::cout << " + benchmark_override_simvars: " << benchmark_override_simvars << std::endl;
			std::cout << " + setup_dealiased: " << setup_dealiased << std::endl;
			std::cout << " + benchmark_galewsky_umax: " << benchmark_galewsky_umax << std::endl;
			std::cout << " + benchmark_galewsky_hamp: " << benchmark_galewsky_hamp << std::endl;
			std::cout << " + benchmark_galewsky_phi2: " << benchmark_galewsky_phi2 << std::endl;
			std::cout << " + benchmark_normal_modes_case: " << benchmark_normal_modes_case << std::endl;
			std::cout << " + object_scale: " << object_scale << std::endl;
			std::cout << " + object_coord_x: " << object_coord_x << std::endl;
			std::cout << " + object_coord_y: " << object_coord_y << std::endl;
			std::cout << " + sphere_advection_rotation_angle: " << sphere_advection_rotation_angle << std::endl;
			std::cout << " + input_data_filenames:" << std::endl;
			std::cout << std::endl;
		}


		void outputProgParams()
		{
			std::cout << std::endl;
			std::cout << "SIMULATION SETUP PARAMETERS:" << std::endl;
			std::cout << "	--random-seed [int]		random seed for random number generator" << std::endl;
			std::cout << "	--benchmark-name [string]	benchmark name" << std::endl;
			std::cout << "	--benchmark-override-simvars [bool]	Allow overwriting simulation variables by benchmark (default: 1)" << std::endl;
			std::cout << "	--benchmark-setup-dealiased [bool]	Use dealiasing for setup (default: 1)" << std::endl;
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


		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{
	        long_options[next_free_program_option] = {"random-seed", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"initial-coord-x", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"initial-coord-y", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"advection-rotation-angle", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"benchmark-name", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"benchmark-override-simvars", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"benchmark-setup-dealiased", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;


	        long_options[next_free_program_option] = {"benchmark-galewsky-umax", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"benchmark-galewsky-hamp", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"benchmark-galewsky-phi2", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

			long_options[next_free_program_option] = {"benchmark-normal-modes-case", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;
		}



		/*
		 * This method is called to parse a particular
		 * long option related to some ID.
		 *
		 * \return: -1 if the option has been processed
		 */
		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{

			switch(i_option_index)
			{
			case 0:
				random_seed = atoi(i_value);
				return -1;

			case 1:
				object_coord_x = atof(i_value);
				return -1;

			case 2:
				object_coord_y = atof(i_value);
				return -1;

			case 3:
				sphere_advection_rotation_angle = atof(i_value);
				return -1;

			case 4:
				benchmark_name = i_value;
				return -1;

			case 5:
				benchmark_override_simvars = atoi(i_value);
				return -1;

			case 6:
				setup_dealiased = atof(i_value);
				return -1;

			case 7:
				benchmark_galewsky_umax = atof(i_value);
				return -1;

			case 8:
				benchmark_galewsky_hamp = atof(i_value);
				return -1;

			case 9:
				benchmark_galewsky_phi2 = atof(i_value);
				return -1;

			case 10:
				benchmark_normal_modes_case = i_value;
				return -1;
			}

			return 11;
		}
	} benchmark;



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


#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		/**
		 * Earth radius for simulations on the sphere
		 */
		double sphere_radius = 6.37122e6;

		/**
		 * Simulation on f-sphere? (constant f0 term over entire sphere)
		 */
		bool sphere_use_fsphere = false;

		/**
		 * Coriolis effect
		 * 7.2921 x 10^{-5}
		 */
		double sphere_rotating_coriolis_omega = 0.000072921;

		double sphere_fsphere_f0 = 0.00007292*2; //Sphere
#endif


		/**
		 * Plane with f-Coriolis rotation
		 */
#if SWEET_USE_PLANE_TEST
		double plane_rotating_f0 = 1.0; //Plane
#endif


		/**
		 * Gravitational constant
		 */
		double gravitation = 9.80616;


		/**
		 * domain size if running simulation on the plane
		 */
		double plane_domain_size[2] = {1.0, 1.0};


		/**
		 * Velocity and additional parameter for advection test cases
		 */
		double advection_velocity[3] = {0, 0, 0};


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "SIMULATION COEFFICIENTS:" << std::endl;
			std::cout << " + h0: " << h0 << std::endl;
			std::cout << " + viscosity: " << viscosity << std::endl;
			std::cout << " + viscosity_order: " << viscosity_order << std::endl;
			std::cout << " + gravitation: " << gravitation << std::endl;
			std::cout << " + domain_size (2D): " << plane_domain_size[0] << " x " << plane_domain_size[1] << std::endl;
			std::cout << " + advection_velocity (x, y, rotation speed): " << advection_velocity[0] << ", " << advection_velocity[1] << ", " << advection_velocity[2] << std::endl;

#if SWEET_USE_PLANE_TEST
			std::cout << " + plane_rotating_f0: " << plane_rotating_f0 << std::endl;
#endif

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
			std::cout << " + sphere_radius: " << sphere_radius << std::endl;
			std::cout << " + sphere_rotating_coriolis_omega: " << sphere_rotating_coriolis_omega << std::endl;
			std::cout << " + sphere_use_fsphere: " << sphere_use_fsphere << std::endl;
			std::cout << " + sphere_fsphere_f0: " << sphere_fsphere_f0 << std::endl;
#endif

			std::cout << std::endl;
		}


		void outputProgParams()
		{
			std::cout << "Simulation parameters:" << std::endl;
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


		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{
	        // sim
	        long_options[next_free_program_option] = {"advection-velocity", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;
		}


		/*
		 * This method is called to parse a particular
		 * long option related to some ID.
		 *
		 * \return: -1 if the option has been processed
		 */
		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{
			switch(i_option_index)
			{
			case 0:
				split3double(i_value, &advection_velocity[0], &advection_velocity[1], &advection_velocity[2]);
				return -1;
			}

			return 1;
		}
	} sim;



	/**
	 * This class stored the discretization-related parameters
	 *
	 * resolution / timestepping
	 */
	struct Discretization
	{
		/**
		 * resolution in physical space (grid cells)
		 */
		int space_res_physical[2] = {0, 0};


		/**
		 * resolution in spectral space (number of modes)
		 */
		int space_res_spectral[2] = {0, 0};


		/**
		 * use spectral differential operators
		 */
		bool space_use_spectral_basis_diffs =
#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_SPHERE_SPECTRAL_SPACE
				true;
#else
				false;
#endif

		/**
		 * Use C-grid staggering
		 */
		bool space_grid_use_c_staggering = false;

		/// Leapfrog: Robert Asselin filter
		double timestepping_leapfrog_robert_asselin_filter = 0;

		/// Crank-Nicolson filter
		double timestepping_crank_nicolson_filter = 0.5;

		/// Number of iterations for semi-Lagrangian methods
		int semi_lagrangian_max_iterations = 2;

		/// Method to use for computing depature points
		std::string semi_lagrangian_departure_point_method = "settls";

		/// Use limiter for higher-order SL interpolation (cubic)
		bool semi_lagrangian_interpolation_limiter = false;

		/// Create pseudo points at poles for interpolation
		bool semi_lagrangian_sampler_use_pole_pseudo_points = false;

		/// Convergence threshold for semi-Lagrangian methods (set to -1 to ignore error)
		double semi_lagrangian_convergence_threshold = -1;

		/// Use accurate spherical geometry (???) or approximation (Ritchie 1995)
		double semi_lagrangian_approximate_sphere_geometry = 0;


		/// String of time stepping method
		/// See doc/swe/swe_plane_timesteppings
		std::string timestepping_method;

		/// Order of time stepping
		int timestepping_order = -1;

		/// Order of 2nd time stepping which might be used
		int timestepping_order2 = -1;



		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "DISCRETIZATION:" << std::endl;
			std::cout << " + space_res_physical: " << space_res_physical[0] << " x " << space_res_physical[1] << std::endl;
			std::cout << " + space_res_spectral: " << space_res_spectral[0] << " x " << space_res_spectral[1] << std::endl;
			std::cout << " + space_use_spectral_basis_diffs: " << space_use_spectral_basis_diffs << std::endl;
			std::cout << " + space_grid_use_c_staggering: " << space_grid_use_c_staggering << std::endl;
			std::cout << " + timestepping_method: " << timestepping_method << std::endl;
			std::cout << " + timestepping_order: " << timestepping_order << std::endl;
			std::cout << " + timestepping_order2: " << timestepping_order2 << std::endl;
			std::cout << " + timestepping_leapfrog_robert_asselin_filter: " << timestepping_leapfrog_robert_asselin_filter << std::endl;
			std::cout << " + timestepping_crank_nicolson_filter: " << timestepping_crank_nicolson_filter << std::endl;
			std::cout << " + semi_lagrangian_max_iterations: " << semi_lagrangian_max_iterations << std::endl;
			std::cout << " + semi_lagrangian_departure_point_method: " << semi_lagrangian_departure_point_method << std::endl;
			std::cout << " + semi_lagrangian_interpolation_limiter: " << semi_lagrangian_interpolation_limiter << std::endl;
			std::cout << " + semi_lagrangian_sampler_use_pole_pseudo_points: " << semi_lagrangian_sampler_use_pole_pseudo_points << std::endl;
			std::cout << " + semi_lagrangian_convergence_threshold: " << semi_lagrangian_convergence_threshold << std::endl;
			std::cout << " + semi_lagrangian_approximate_sphere_geometry: " << semi_lagrangian_approximate_sphere_geometry << std::endl;
			std::cout << " + plane_dealiasing (compile time): " <<
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
			std::cout << "	--space-grid-use-c-staggering [0/1]	Use staggering" << std::endl;
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
			std::cout << "	--timestepping-method [string]	String of time stepping method" << std::endl;
			std::cout << "	--timestepping-order [int]			Specify the order of the time stepping" << std::endl;
			std::cout << "	--timestepping-order2 [int]			Specify the order of the time stepping" << std::endl;
			std::cout << "	--leapfrog-robert-asselin-filter [0;1]		Damping parameter for Robert-Asselin filter" << std::endl;
			std::cout << "	--normal-mode-analysis-generation [0;1;2;3]	Generate output data for normal mode analysis" << std::endl;
			std::cout << "							0: don't generate" << std::endl;
			std::cout << "							1: generate in physical space" << std::endl;
			std::cout << "							2: generate in spectral space" << std::endl;
			std::cout << "							3: generate in spectral space with complex matrix" << std::endl;
			std::cout << "	--semi-lagrangian-max-iterations [int]		Number of max. iterations during semi-Lagrangian time integration" << std::endl;
			std::cout << "	--semi-lagrangian-departure-point-method [str]		'settls' (default), 'midpoint', 'std'" << std::endl;
			std::cout << "	--semi-lagrangian-sampler-use-pole-pseudo-points [bool]" << std::endl;
			std::cout << "								Use pseudo poles for sampling" << std::endl;
			std::cout << "									false (default)" << std::endl;
			std::cout << "	--semi-lagrangian-interpolation-limiter [bool]	Use limiter for cubic interpolation" << std::endl;
			std::cout << "	--semi-lagrangian-convergence-threshold [float]	Threshold to stop iterating, Use -1 to disable" << std::endl;
			std::cout << "	--semi-lagrangian-approximate-sphere-geometry [int]	0: no approximation, 1: Richies approximation, default: 0" << std::endl;

		}



		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{

	        // DISC
	        long_options[next_free_program_option] = {"timestepping-method", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"timestepping-order", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"timestepping-order2", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"leapfrog-robert-asselin-filter", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"crank-nicolson-filter", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"semi-lagrangian-max-iterations", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"semi-lagrangian-departure-point-method", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"semi-lagrangian-sampler-use-pole-pseudo-points", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"semi-lagrangian-interpolation-limiter", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"semi-lagrangian-convergence-threshold", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"semi-lagrangian-approximate-sphere-geometry", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"space-grid-use-c-staggering", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;
		}



		/*
		 * This method is called to parse a particular
		 * long option related to some ID.
		 *
		 * \return: -1 if the option has been processed
		 */
		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{

			switch(i_option_index)
			{
			case 0:
				timestepping_method = i_value;
				return -1;

			case 1:
				timestepping_order = atoi(i_value);
				return -1;

			case 2:
				timestepping_order2 = atoi(i_value);
				return -1;

			case 3:
				timestepping_leapfrog_robert_asselin_filter = atof(i_value);
				return -1;

			case 4:
				timestepping_crank_nicolson_filter = atof(i_value);
				return -1;

			case 5:
				semi_lagrangian_max_iterations = atoi(i_value);
				return -1;

			case 6:
				semi_lagrangian_departure_point_method = i_value;
				return -1;

			case 7:
				semi_lagrangian_sampler_use_pole_pseudo_points = atoi(i_value);
				return -1;

			case 8:
				semi_lagrangian_interpolation_limiter = atoi(i_value);
				return -1;

			case 9:
				semi_lagrangian_convergence_threshold = atof(i_value);
				return -1;

			case 10:
				semi_lagrangian_approximate_sphere_geometry = atof(i_value);
				return -1;

			case 11:
				space_grid_use_c_staggering = atof(i_value);
				return -1;
			}

			return 12;
		}
	} disc;



	/**
	 * program parameters without specific association.
	 * These variables can be used different for each program
	 */
	struct UserDefined
	{
		std::string var[20];
	} user_defined;

	/**
	 * Miscellaneous variables
	 */
	struct Misc
	{
		/// set verbosity of simulation
		int verbosity = 0;

		/// compute errors
		int compute_errors = 0;

		/// do instability checks for simulation
		int instability_checks = 1;

		/// activate GUI mode?
		bool gui_enabled = (SWEET_GUI == 0 ? false : true);


		/// id for visualization
		int vis_id = 0;

		/// Diffusion applied only on nonlinear divergence
		int use_nonlinear_only_visc = 0;

		/// Load / Save plans for SHTNS (useful for reproducibility)
		TransformationPlans::TRANSFORMATION_PLAN_CACHE reuse_spectral_transformation_plans = TransformationPlans::QUICK;

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		int normal_mode_analysis_generation = 0;

		/*
		 * Some flexible variable where one can just add options like
		 * --comma-separated-tags=galewsky_analytical_geostrophic_setup
		 */
		std::string comma_separated_tags = "";

		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "MISC:" << std::endl;
			std::cout << " + verbosity: " << verbosity << std::endl;
			std::cout << " + compute_errors " << compute_errors << std::endl;
			std::cout << " + instability_checks: " << instability_checks << std::endl;
			std::cout << " + gui_enabled: " << gui_enabled << std::endl;
			std::cout << " + vis_id: " << vis_id << std::endl;
			std::cout << " + use_nonlinear_only_visc: " << use_nonlinear_only_visc << std::endl;
			std::cout << " + reuse_spectral_transformation_plans: " << TransformationPlans::getStringFromEnum(reuse_spectral_transformation_plans) << std::endl;
			std::cout << std::endl;
			std::cout << " + normal_mode_analysis_generation: " << normal_mode_analysis_generation << std::endl;
			std::cout << " + comma_separated_tags: " << comma_separated_tags << std::endl;
			std::cout << std::endl;
		}



		void outputProgParams()
		{
			std::cout << "" << std::endl;
			std::cout << "Misc options:" << std::endl;
			std::cout << "	-v [int]			verbosity level" << std::endl;
	//		std::cout << "	-V [double]			period of outputConfig" << std::endl;
			std::cout << "	-G [0/1]			graphical user interface" << std::endl;
			std::cout << "	-O [string]			string prefix for filename of output of simulation data (default output_%s_t%020.8f.csv)" << std::endl;
			std::cout << "	-d [int]			accuracy of floating point output" << std::endl;
			std::cout << "	-i [file0][;file1][;file3]...	string with filenames for initial conditions" << std::endl;
			std::cout << "					specify BINARY; as first file name to read files as binary raw data" << std::endl;
			std::cout << "	--compute-errors [int]          Compute errors when possible [1], default=0	" << std::endl;
			std::cout << "	--use-local-visc [0/1]	Viscosity will be applied only on nonlinear divergence, default:0" << std::endl;
			std::cout << "	--reuse-plans [0/1]	Save plans for fftw transformations and SH transformations" << std::endl;
			std::cout << "					-1: use only estimated plans (no wisdom)" << std::endl;
			std::cout << "					0: compute optimized plans (no wisdom)" << std::endl;
			std::cout << "					1: compute optimized plans, use wisdom if available and store wisdom" << std::endl;
			std::cout << "					2: use wisdom if available if not, trigger error if wisdom doesn't exist (not yet working for SHTNS)" << std::endl;
			std::cout << "					default: -1 (quick mode)" << std::endl;
			std::cout << "" << std::endl;
		}

		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{
	        // MISC
	        long_options[next_free_program_option] = {"compute-errors", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"instability-checks", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"use-nonlinear-only-visc", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"reuse-plans", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"normal-mode-analysis-generation", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;

	        long_options[next_free_program_option] = {"comma-separated-tags", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;
		}



		/*
		 * This method is called to parse a particular
		 * long option related to some ID.
		 *
		 * \return: -1 if the option has been processed
		 */
		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{

			switch(i_option_index)
			{
			case 0:
				compute_errors = atoi(i_value);
				return -1;

			case 1:
				instability_checks = atoi(i_value);
				return -1;

			case 2:
				use_nonlinear_only_visc = atoi(i_value);
				return -1;

			case 3:
				reuse_spectral_transformation_plans = TransformationPlans::getEnumFromString(i_value);
				return -1;

			case 4:
				normal_mode_analysis_generation = atoi(i_value);
				return -1;

			case 5:
				comma_separated_tags = i_value;
				return -1;
			}

			return 6;
		}


	} misc;


	/**
	 * Miscellaneous variables
	 */
	struct Parallelization
	{
		/// number of threads
		int num_threads_space = -1;

		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "PARALLELIZATION:" << std::endl;
			std::cout << " + num_threads: " << num_threads_space << std::endl;
			std::cout << std::endl;
		}



		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{
			long_options[next_free_program_option] = {"num-threads-space", required_argument, 0, 256+next_free_program_option};
	        next_free_program_option++;
		}

		void outputProgParams()
		{
			std::cout << "" << std::endl;
			std::cout << "Parallelization options:" << std::endl;
			std::cout << "	--num-threads-space [int]			Specify how many threads to use for spatial parallelization (very useful for nested parallel regions)" << std::endl;
		}

		/*
		 * This method is called to parse a particular
		 * long option related to some ID.
		 *
		 * \return: -1 if the option has been processed
		 */
		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{
			switch(i_option_index)
			{
			case 0:
				num_threads_space = atoi(i_value);
				return -1;

			}

			return 0;
		}


	} parallelization;


	/**
	 * Timestepping
	 */
	struct TimestepControl
	{
		/// Continue running simulation timestepping.
		/// This is beneficial to pause simulations if driven interactively.
		bool run_simulation_timesteps = true;

		/// Number of simulated time steps
		int current_timestep_nr = 0;

		/// Time step size used during setup
		double setup_timestep_size = -1;

		/// Current time step size
		double current_timestep_size = -1;

		/// Time in simulation
		double current_simulation_time = 0;

		/// Maximum number of time steps to simulate
		int max_timesteps_nr = std::numeric_limits<int>::max();

		/// Maximum simulation time to execute the simulation for
		double max_simulation_time = std::numeric_limits<double>::infinity();

		/// Maximum wallclock time to execute the simulation for
		double max_wallclock_time = -1;


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "TIMECONTROL:" << std::endl;
			std::cout << " + run_simulation_timesteps: " << run_simulation_timesteps << std::endl;
			std::cout << " + current_timestep_nr: " << current_timestep_nr << std::endl;
			std::cout << " + setup_timestep_size: " << setup_timestep_size << std::endl;
			std::cout << " + current_timestep_size: " << current_timestep_size << std::endl;
			std::cout << " + current_simulation_time: " << current_simulation_time << std::endl;
			std::cout << " + max_timesteps_nr: " << max_timesteps_nr << std::endl;
			std::cout << " + max_simulation_time: " << max_simulation_time << std::endl;
			std::cout << " + max_wallclock_time: " << max_wallclock_time << std::endl;
			std::cout << std::endl;
		}


		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{
			long_options[next_free_program_option] = {"dt", required_argument, 0, 256+next_free_program_option};
			next_free_program_option++;

			long_options[next_free_program_option] = {"max-wallclock-time", required_argument, 0, 256+next_free_program_option};
			next_free_program_option++;
		}



		/*
		 * This method is called to parse a particular
		 * long option related to some ID.
		 *
		 * \return: -1 if the option has been processed
		 */
		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{
			switch(i_option_index)
			{
			case 0:
				current_timestep_size = atof(i_value);
				setup_timestep_size = current_timestep_size;
				return -1;

			case 1:
				max_wallclock_time = atof(i_value);
				return -1;
			}

			return 2;
		}

	} timecontrol;

	/**
	 * SDC parameters
	 */
	struct SDC {

		std::string fileName = "";
		int nNodes = 3;	

		// Nodes values (between 0 and 1)
		sweet::DictArrayND<1, double> nodes;

		// Quadrature weights
		sweet::DictArrayND<1, double> weights;

		// Collocation matrix
		sweet::DictArrayND<2, double> qMatrix;

		// QDelta matrix for implicit sweep
		sweet::DictArrayND<2, double> qDeltaI;

		// QDelta matrix for explicit sweep
		sweet::DictArrayND<2, double> qDeltaE;

		// QDelta matrix for initial (implicit) sweep
		sweet::DictArrayND<2, double> qDelta0;

		// Number of iterations (sweeps)
		sweet::Dict::int64 nIter=0;

		// Type of initial sweep to use
		std::string initSweepType;

		// Wether or not use the diagonal implementation
		sweet::Dict::int64 diagonal=0;

		// Wether or not use collocation update for end point
		sweet::Dict::int64 useEndUpdate=0;

		// Unique string ID
		std::string idString;	

		void outputProgParams()
		{
			std::cout << std::endl;
			std::cout << "SDC option:" << std::endl;
			std::cout << "	--sdc-file [path]   SDC parameters in sweet::Dict format" << std::endl;
			std::cout << std::endl;
		}


		void outputConfig()
		{
			std::cout << std::endl;
			std::cout << "SDC:" << std::endl;
			std::cout << " + nNodes: " << nNodes << std::endl;
			std::cout << " + nodes: " << nodes << std::endl;
			std::cout << " + weights: " << weights << std::endl;
			std::cout << " + qMatrix: " << qMatrix << std::endl;
			std::cout << " + qDeltaI: " << qDeltaI << std::endl;
			std::cout << " + qDeltaE: " << qDeltaE << std::endl;
			std::cout << " + qDelta0: " << qDelta0 << std::endl;
			std::cout << " + nIter: " << nIter << std::endl;
			std::cout << " + initSweepType: " << initSweepType << std::endl;
			std::cout << " + diagonal: " << diagonal << std::endl;
			std::cout << " + useEndUpdate: " << useEndUpdate << std::endl;
			std::cout << std::endl;
		}


		void setup_longOptionsList(
				struct option *long_options,
				int &next_free_program_option
		)
		{
			long_options[next_free_program_option] = {"sdc-file", optional_argument, 0, 256+next_free_program_option};
			next_free_program_option++;
		}


		int setup_longOptionValue(
				int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
				const char *i_value		///< Value in string format
		)
		{
			switch(i_option_index)
			{
			case 0:
				fileName = i_value;
				sweet::Dict params(fileName);
				params.get("nodes", nodes);
				nNodes = nodes.size();
				params.get("weights", weights);
				params.get("qMatrix", qMatrix);
				params.get("qDeltaI", qDeltaI);
				params.get("qDeltaE", qDeltaE);
				params.get("qDelta0", qDelta0);
				params.get("nIter", nIter);
				params.get("diagonal", diagonal);
				params.get("initSweepType", initSweepType);
				params.get("useEndUpdate", useEndUpdate);
				params.get("idString", idString);
				return -1;
			}

			return 2;
		}
	} sdc;


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

	static
	int split2int(
			const char *i_str,
			int *o_int0,
			int *o_int1
	)
	{
		std::vector<std::string> res = StringSplit::split(i_str, ",");

		int c = res.size();

		if (c == 0)
			SWEETError("Invalid format for modes");

		if (c == 1)
		{
			*o_int0 = atoi(res[0].c_str());
			return 1;
		}
		else if (c == 2)
		{
			*o_int0 = atoi(res[0].c_str());
			*o_int1 = atoi(res[1].c_str());
			return 2;
		}

		SWEETError("More than 2 values given");
		return -1;
	}


	static
	int split2double(
			const char *i_str,
			double *o_int0,
			double *o_int1
	)
	{
		std::vector<std::string> res = StringSplit::split(i_str, ",");

		int c = res.size();

		if (c == 0)
			SWEETError("Invalid format for modes");

		if (c == 1)
		{
			*o_int0 = atof(res[0].c_str());
			return 1;
		}
		else if (c == 2)
		{
			*o_int0 = atof(res[0].c_str());
			*o_int1 = atof(res[1].c_str());
			return 2;
		}

		SWEETError("More than 2 values given");
		return -1;
	}


	static
	int split3double(
			const char *i_str,
			double *o_int0,
			double *o_int1,
			double *o_int2
	)
	{
		std::vector<std::string> res = StringSplit::split(i_str, ",");

		int c = res.size();

		if (c == 0)
			SWEETError("Invalid format for modes");

		if (c == 1)
		{
			*o_int0 = atof(res[0].c_str());
			return 1;
		}
		else if (c == 2)
		{
			*o_int0 = atof(res[0].c_str());
			*o_int1 = atof(res[1].c_str());
			return 2;
		}
		else if (c == 3)
		{
			*o_int0 = atof(res[0].c_str());
			*o_int1 = atof(res[1].c_str());
			*o_int2 = atof(res[2].c_str());
			return 3;
		}

		SWEETError("More than 3 values given");
		return -1;
	}


	void print_params(const char *i_user_defined_program_parameters[])
	{

		sim.outputProgParams();
		benchmark.outputProgParams();
		parallelization.outputProgParams();
		disc.outputProgParams();
		iodata.outputProgParams();

		// TODO : cleanup ??
		std::cout << "" << std::endl;
		std::cout << "Control:" << std::endl;
		std::cout << "	--dt [time]	timestep size, default=?" << std::endl;
		std::cout << "	--max-wallclock-time [time]	wallclock time limitation, default=-1" << std::endl;
		std::cout << "	-t [time]	maximum simulation time, default=-1 (infinity)" << std::endl;
		std::cout << "	-T [stepnr]	maximum number of time steps, default=-1 (infinity)" << std::endl;
		std::cout << "	-o [time]	time interval at which output should be written, (set to 0 for output at every time step), default=-1 (no output) " << std::endl;

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
					int c = split2int(optarg, &disc.space_res_physical[0], &disc.space_res_physical[1]);
					if (c == 1)
						disc.space_res_physical[1] = disc.space_res_physical[0];
				}
				break;

			case 'M':
				{
					int c = split2int(optarg, &disc.space_res_spectral[0], &disc.space_res_spectral[1]);
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
#if SWEET_USE_PLANE_TEST
				sim.plane_rotating_f0 = atof(optarg);
#endif
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
