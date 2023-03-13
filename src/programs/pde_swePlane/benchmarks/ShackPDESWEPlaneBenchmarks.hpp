/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACK_PDE_SWE_PLANE_BENCHMARKS_HPP_
#define SRC_INCLUDE_SWEET_SHACK_PDE_SWE_PLANE_BENCHMARKS_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>



/**
 * simulation coefficients
 */
class ShackPDESWEPlaneBenchmarks	:
		public sweet::ShackInterface
{
public:

	/// seed for random number generator
	int random_seed = 0;

	/// benchmark scenario
	std::string benchmark_name = "";

	/// Normal modes benchmark scenario
	std::string benchmark_normal_modes_case = "";

	/// radius
	double object_scale = 1;

	/// setup coordinate of e.g. radial breaking dam, x-placement \in [0;1]
	double object_coord_x = 0.5;

	/// setup coordinate of e.g. radial breaking dam, y-placement \in [0;1]
	double object_coord_y = 0.5;

	/// load external forces if available from benchmark scenario
	void (*getExternalForcesCallback)(int, double, sweet::PlaneData_Spectral*, ShackPDESWEPlaneBenchmarks*) = nullptr;// = &fun_no_forces;		/// SET TO NULLPTR
	ShackPDESWEPlaneBenchmarks *getExternalForcesUserData = nullptr;


	/**
	 * Velocity and additional parameter for advection test cases
	 */
	double advection_velocity[2] = {0, 0};


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "SIMULATION SETUP PARAMETERS:" << std::endl;
		std::cout << i_prefix << "	--random-seed [int]		random seed for random number generator" << std::endl;
		std::cout << i_prefix << "	--benchmark-name [string]	benchmark name" << std::endl;
		std::cout << i_prefix << "	-x [float]				x coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-y [float]				y coordinate for setup \\in [0;1], default=0.5" << std::endl;
		std::cout << i_prefix << "	-r [radius]				scale factor of radius for initial condition, default=1" << std::endl;
		std::cout << i_prefix << "	--advection-velocity=[float],[float],[float]	advection velocity components (x, y, rotational)" << std::endl;

		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--random-seed", random_seed);
		i_pa.getArgumentValueBy2Keys("--initial-coord-x", "-x", object_coord_x);
		i_pa.getArgumentValueBy2Keys("--initial-coord-y", "-y", object_coord_y);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

		i_pa.getArgumentValueByKey("--benchmark-normal-modes-case", benchmark_normal_modes_case);
		i_pa.getArgumentValueByKey("-r", object_scale);

		std::string tmp;
		if (i_pa.getArgumentValueByKey("--advection-velocity", tmp))
			StringSplit::split2double(tmp, &advection_velocity[0], &advection_velocity[1]);


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
		std::cout << " + benchmark_normal_modes_case: " << benchmark_normal_modes_case << std::endl;
		std::cout << " + object_scale: " << object_scale << std::endl;
		std::cout << " + object_coord_x: " << object_coord_x << std::endl;
		std::cout << " + object_coord_y: " << object_coord_y << std::endl;
		std::cout << i_prefix << " + advection_velocity (x, y, rotation speed): " << advection_velocity[0] << ", " << advection_velocity[1] << std::endl;
		std::cout << std::endl;
	}
};




#endif
