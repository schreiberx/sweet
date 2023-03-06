/*
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACK_PDE_ADVECTION_SPHERE_BENCHMARK_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACK_PDE_ADVECTION_SPHERE_BENCHMARK_HPP_

#include <string>
#include <iostream>
#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/ScalarDataArray.hpp>

/**
 * Values and parameters to setup benchmarks simulations
 */
class ShackPDEAdvectionSphereBenchmarks	:
		public sweet::ShackInterface
{
public:
	/// seed for random number generator
	int random_seed = 0;

	/// benchmark scenario
	std::string benchmark_name = "";

	/// rotation angle for advection equation
	double sphere_advection_rotation_angle = 0;

	/**
	 * Flag to indicate the presence of topography
	 */
	bool use_topography = false;


	/*
	 * Get updated velocities for particular point in time
	 */
	void (*getVelocities)(
			sweet::SphereData_Physical&,
			sweet::SphereData_Physical&,
			double i_time,
			ShackPDEAdvectionSphereBenchmarks* user_ptr
	) = nullptr;
	ShackPDEAdvectionSphereBenchmarks *getVelocitiesUserData = nullptr;

	/**
	 * Callback for special benchmark
	 */
	void (*callback_slComputeDeparture3rdOrder)(
			void *i_this,
			const sweet::ScalarDataArray &i_pos_lon_A,	///< longitude coordinate to compute the velocity for
			const sweet::ScalarDataArray &i_pos_lat_A,	///< latitude coordinate to compute the velocity for
			sweet::ScalarDataArray &o_pos_lon_D,		///< velocity along longitude
			sweet::ScalarDataArray &o_pos_lat_D,		///< velocity along latitude
			double i_dt,
			double i_timestamp_arrival			///< timestamp at arrival point
	);
	void *slComputeDeparture3rdOrderUserData = nullptr;


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
		std::cout << i_prefix << "	--initial-coord-x [float]		Same as -x" << std::endl;
		std::cout << i_prefix << "	--initial-coord-y [float]		Same as -y" << std::endl;
		std::cout << i_prefix << "	--advection-rotation-angle [float]	Rotation angle for e.g. advection test case" << std::endl;
		std::cout << i_prefix << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--random-seed", random_seed);
		i_pa.getArgumentValueByKey("--benchmark-name", benchmark_name);

		std::string tmp;
		if (i_pa.getArgumentValueByKey("--advection-velocity", tmp))
			StringSplit::split3double(tmp, &advection_velocity[0], &advection_velocity[1], &advection_velocity[2]);

		i_pa.getArgumentValueByKey("--advection-rotation-angle", sphere_advection_rotation_angle);

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
		std::cout << i_prefix << " + sphere_advection_rotation_angle: " << sphere_advection_rotation_angle << std::endl;
		std::cout << i_prefix << " + input_data_filenames:" << std::endl;
		std::cout << i_prefix << std::endl;
	}
};


#endif

