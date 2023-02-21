/*
 * ShackBenchmark.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_

#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>


/**
 * Values and parameters to setup benchmarks simulations
 */
class Benchmark	:
		public sweet::ClassDictionaryInterface
{
public:
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


	/// load external forces if available from benchmark scenario
	void (*getExternalForcesCallback)(int, double, void*, void*) = nullptr;// = &fun_no_forces;		/// SET TO NULLPTR
	void *getExternalForcesUserData = nullptr;



	void outputConfig()
	{
		printClass();
	}


	void outputProgParams()
	{
		printProgramArguments();
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

	void printProgramArguments(const std::string& i_prefix = "")
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

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--random-seed", random_seed);
		i_pa.getArgumentValueBy2Keys("--initial-coord-x", "-x", object_coord_x);
		i_pa.getArgumentValueBy2Keys("--initial-coord-y", "-y", object_coord_y);
		i_pa.getArgumentValueByKey("--advection-rotation-angle", sphere_advection_rotation_angle);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);
		i_pa.getArgumentValueByKey("--benchmark-override-simvars", benchmark_override_simvars);

		i_pa.getArgumentValueByKey("--benchmark-setup-dealiased", setup_dealiased);
		i_pa.getArgumentValueByKey("--benchmark-galewsky-umax", benchmark_galewsky_umax);
		i_pa.getArgumentValueByKey("--benchmark-galewsky-hamp", benchmark_galewsky_hamp);
		i_pa.getArgumentValueByKey("--benchmark-galewsky-phi2", benchmark_galewsky_phi2);
		i_pa.getArgumentValueByKey("--benchmark-override-simvars", benchmark_override_simvars);
		i_pa.getArgumentValueByKey("--benchmark-normal-modes-case", benchmark_normal_modes_case);
		i_pa.getArgumentValueByKey("-r", object_scale);


		if (error.exists())
			return error.forwardFromWithPositiveReturn(i_pa.error);

		if (random_seed >= 0)
			srandom(random_seed);

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
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
};



#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_ */
