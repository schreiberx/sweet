/*
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_

#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
	#include <sweet/sphere/SphereData_Spectral.hpp>
#endif

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


#if SWEET_USE_SPHERE_SPECTRAL_SPACE
	/**
	 * Topography vector
	 */
	SphereData_Physical h_topo;
#endif


	/// load external forces if available from benchmark scenario
	void (*getExternalForcesCallback)(int, double, PlaneData_Spectral*, ShackPDEAdvectionPlaneBenchmarks*) = nullptr;
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
		std::cout << std::endl;
		std::cout << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << "	--random-seed [int]		random seed for random number generator" << std::endl;
		std::cout << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << "	--advection-velocity=[float],[float],[float]	advection velocity components (x, y, rotational)" << std::endl;
		std::cout << "	--benchmark-override-simvars [bool]	Allow overwriting simulation variables by benchmark (default: 1)" << std::endl;
		std::cout << "	--benchmark-setup-dealiased [bool]	Use dealiasing for setup (default: 1)" << std::endl;
		std::cout << "	-x [float]				x coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << "	-y [float]				y coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << "	-r [radius]				scale factor of radius for initial condition, default=1" << std::endl;
		std::cout << "	--initial-freq-x-mul [float]		Frequency for the waves initial conditions in x, default=2" << std::endl;
		std::cout << "	--initial-freq-y-mul [float]		Frequency for the waves initial conditions in y, default=1" << std::endl;
		std::cout << "	--initial-coord-x [float]		Same as -x" << std::endl;
		std::cout << "	--initial-coord-y [float]		Same as -y" << std::endl;

		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--random-seed", random_seed);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);


		std::string tmp;
		if (i_pa.getArgumentValueByKey("--advection-velocity", tmp))
		{
			StringSplit::split3double(tmp, &advection_velocity[0], &advection_velocity[1], &advection_velocity[2]);
		}

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
		std::cout << std::endl;
		std::cout << "BENCHMARK:" << std::endl;
		std::cout << " + random_seed: " << random_seed << std::endl;
		std::cout << " + benchmark_name: " << benchmark_name << std::endl;
		std::cout << " + advection_velocity (x, y, rotation speed): " << advection_velocity[0] << ", " << advection_velocity[1] << ", " << advection_velocity[2] << std::endl;
		std::cout << " + benchmark_override_simvars: " << benchmark_override_simvars << std::endl;
		std::cout << " + setup_dealiased: " << setup_dealiased << std::endl;
		std::cout << " + object_scale: " << object_scale << std::endl;
		std::cout << " + object_coord_x: " << object_coord_x << std::endl;
		std::cout << " + object_coord_y: " << object_coord_y << std::endl;
		std::cout << " + input_data_filenames:" << std::endl;
		std::cout << std::endl;
	}
};



#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKBENCHMARK_HPP_ */
