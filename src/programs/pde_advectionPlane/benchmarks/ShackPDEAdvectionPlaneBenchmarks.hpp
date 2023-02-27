/*
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_

#include <string>
#include <iostream>
#include <sweet/core/shacks/ShackInterface.hpp>

/**
 * Values and parameters to setup benchmarks simulations
 */
class ShackPDEAdvectionPlaneBenchmarks	:
		public sweet::ShackInterface
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


	/// radius
	double object_scale = 1;

	/// setup coordinate of e.g. radial breaking dam, x-placement \in [0;1]
	double object_coord_x = 0.5;

	/// setup coordinate of e.g. radial breaking dam, y-placement \in [0;1]
	double object_coord_y = 0.5;


	/**
	 * Flag to indicate the presence of topography
	 */
	bool use_topography = false;


	/// load external forces if available from benchmark scenario
	void (*getExternalForcesCallback)(int, double, sweet::PlaneData_Spectral*, ShackPDEAdvectionPlaneBenchmarks*) = nullptr;
	void *getExternalForcesUserData = nullptr;



	/**
	 * Velocity and additional parameter for advection test cases
	 */
	double advection_velocity[3] = {0, 0, 0};

	bool validateNonzeroAdvection()
	{
		if (advection_velocity[0] == 0 && advection_velocity[1] == 0)
			return error.set("Both advection velocities are 0, use --advection-velocity=...");

		return true;
	}

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << i_prefix << "	--random-seed [int]		random seed for random number generator" << std::endl;
		std::cout << i_prefix << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << i_prefix << "	--advection-velocity=[float],[float],[float]	advection velocity components (x, y, rotational)" << std::endl;
		std::cout << i_prefix << "	--benchmark-override-simvars [bool]	Allow overwriting simulation variables by benchmark (default: 1)" << std::endl;
		std::cout << i_prefix << "	--benchmark-setup-dealiased [bool]	Use dealiasing for setup (default: 1)" << std::endl;
		std::cout << i_prefix << "	-x [float]				x coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-y [float]				y coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-r [radius]				scale factor of radius for initial condition, default=1" << std::endl;
		std::cout << i_prefix << "	--initial-freq-x-mul [float]		Frequency for the waves initial conditions in x, default=2" << std::endl;
		std::cout << i_prefix << "	--initial-freq-y-mul [float]		Frequency for the waves initial conditions in y, default=1" << std::endl;
		std::cout << i_prefix << "	--initial-coord-x [float]		Same as -x" << std::endl;
		std::cout << i_prefix << "	--initial-coord-y [float]		Same as -y" << std::endl;
		std::cout << i_prefix << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--random-seed", random_seed);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

		std::string tmp;
		if (i_pa.getArgumentValueByKey("--advection-velocity", tmp))
			StringSplit::split3double(tmp, &advection_velocity[0], &advection_velocity[1], &advection_velocity[2]);

		i_pa.getArgumentValueBy2Keys("--initial-coord-x", "-x", object_coord_x);
		i_pa.getArgumentValueBy2Keys("--initial-coord-y", "-y", object_coord_y);

		i_pa.getArgumentValueByKey("--benchmark-setup-dealiased", setup_dealiased);
		i_pa.getArgumentValueByKey("--benchmark-override-simvars", benchmark_override_simvars);
		i_pa.getArgumentValueByKey("-r", object_scale);


		if (error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		if (random_seed >= 0)
			srandom(random_seed);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "BENCHMARK:" << std::endl;
		std::cout << i_prefix << " + random_seed: " << random_seed << std::endl;
		std::cout << i_prefix << " + benchmark_name: " << benchmark_name << std::endl;
		std::cout << i_prefix << " + advection_velocity (x, y, rotation speed): " << advection_velocity[0] << ", " << advection_velocity[1] << ", " << advection_velocity[2] << std::endl;
		std::cout << i_prefix << " + benchmark_override_simvars: " << benchmark_override_simvars << std::endl;
		std::cout << i_prefix << " + setup_dealiased: " << setup_dealiased << std::endl;
		std::cout << i_prefix << " + object_scale: " << object_scale << std::endl;
		std::cout << i_prefix << " + object_coord_x: " << object_coord_x << std::endl;
		std::cout << i_prefix << " + object_coord_y: " << object_coord_y << std::endl;
		std::cout << i_prefix << " + input_data_filenames:" << std::endl;
		std::cout << i_prefix << std::endl;
	}
};



#endif
