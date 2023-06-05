/*
 * ShackTimestepControl.hpp
 *
 *  Created on: Feb 21, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKTIMESTEPCONTROL_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKTIMESTEPCONTROL_HPP_


#include <string>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>


namespace sweet
{

#define SHACK_TIMESTEP_CONTROL_EPSILON 1e-13

/**
 * Timestep Control
 */
class ShackTimestepControl	:
		public ShackInterface
{
public:
	/// Continue running simulation timestepping.
	/// This is beneficial to pause simulations if driven interactively.
	bool run_simulation_timesteps = true;

	/// Number of simulated time steps
	int current_timestep_nr = 0;

	/// Time step size used during setup
	//double setup_timestepSize = -1;

	/// Current time step size
	double current_timestepSize = -1;

	/// Time in simulation
	double current_simulation_time = 0;

	/// Maximum number of time steps to simulate
	int max_timesteps_nr = -1;

	/// Maximum simulation time to execute the simulation for
	double max_simulation_time = -1;

	/// Maximum wallclock time to execute the simulation for
	double max_wallclock_time = -1;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Timecontrol:" << std::endl;
		std::cout << "	--dt [float]	timestep size, default=?" << std::endl;
		std::cout << "	--max-wallclock-time [float]	wallclock time limitation, default=-1" << std::endl;
		std::cout << "	-t [float]	maximum simulation time, default=-1 (infinity)" << std::endl;
		std::cout << "	-T [int]	maximum number of time steps, default=-1 (infinity)" << std::endl;
		std::cout << "	-o [float]	time interval at which output should be written, (set to 0 for output at every time step), default=-1 (no output) " << std::endl;
	}

	/*
	 * Check arguments
	 */
	bool validateMaxSimulationTimeOrTimestepNr()
	{
		if (!(max_simulation_time >= 0) && !(max_timesteps_nr >= 0))
			return error.set("You need to set the maximum simulation time using -t [float] or time step numbers using -T [int]");

		return true;
	}

	/*
	 * Check arguments
	 */
	bool validateMaxTimestepNr()
	{
		if (!(max_timesteps_nr >= 0))
			return error.set("You need to set the maximal number of time steps");

		return true;
	}
	/*
	 * Check arguments
	 */
	bool validateMaxSimulationTime()
	{
		if (!(max_simulation_time >= 0))
			return error.set("You need to set the maximum simulation time using -t [float]");

		return true;
	}

	bool validateTimestepSize()
	{
		if (!(current_timestepSize > 0))
			return error.set("Timestep size not set, use --dt=[float]");

		return true;
	}

	/*
	 * Code which we require again and again to be executed before each time step.
	 *
	 * We have a special function for this since we need to keep round-off errors in mind.
	 *
	 * \return true if the current time step size was modified
	 */
	bool timestepHelperStart()
	{
		// If we didn't set max_simulation_time we just continue
		if (max_simulation_time == -1)
			return false;

		// Check whether there might be some numerical issues of round-off errors in the last time step
		double diff = max_simulation_time - (current_simulation_time + current_timestepSize);

		/*
		 * Check if we're not close to the maximum simulation time and return
		 */
		if (diff > SHACK_TIMESTEP_CONTROL_EPSILON*max_simulation_time)
			return false;

		/*
		 *
		 * This might also change the time step size if it's not necessary,
		 * but we avoid yet another if condition.
		 */
		current_timestepSize = max_simulation_time - current_simulation_time;

		return true;
	}

	bool timestepHelperEnd()
	{
		// advance in time
		current_simulation_time += current_timestepSize;
		current_timestep_nr++;

		return false;
	}

	bool isFinalTimestepReached()
	{
		if (max_timesteps_nr >= 0)
		{
			assert(max_timesteps_nr >= 0);
			assert(current_timestep_nr <= max_timesteps_nr);

			if (max_timesteps_nr == current_timestep_nr)
				return true;
		}

		if (max_simulation_time >= 0)
		{
			double diff = max_simulation_time - current_simulation_time;

			if (diff < 0)
				SWEETError("Internal error: This should never happen (diff < 0)");

			if (diff == 0)
			{
				assert(max_simulation_time == current_simulation_time);
				return true;
			}

#if SWEET_DEBUG
			if (diff < SHACK_TIMESTEP_CONTROL_EPSILON*max_simulation_time)
				SWEETError("Internal error: This should never happen (diff > epsilon)");
#endif
		}

		return false;
	}

	bool processProgramArguments(ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--dt", current_timestepSize);
		i_pa.getArgumentValueByKey("--max-wallclock-time", max_wallclock_time);
		i_pa.getArgumentValueByKey("-t", max_simulation_time);
		i_pa.getArgumentValueByKey("-T", max_timesteps_nr);

		//setup_timestepSize = current_timestepSize;

		current_timestep_nr = 0;
		current_simulation_time = 0;
		//current_timestepSize = setup_timestepSize;

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(i_pa);
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "TIMECONTROL:" << std::endl;
		std::cout << " + run_simulation_timesteps: " << run_simulation_timesteps << std::endl;
		std::cout << " + current_timestep_nr: " << current_timestep_nr << std::endl;
		//std::cout << " + setup_timestepSize: " << setup_timestepSize << std::endl;
		std::cout << " + current_timestepSize: " << current_timestepSize << std::endl;
		std::cout << " + current_simulation_time: " << current_simulation_time << std::endl;
		std::cout << " + max_timesteps_nr: " << max_timesteps_nr << std::endl;
		std::cout << " + max_simulation_time: " << max_simulation_time << std::endl;
		std::cout << " + max_wallclock_time: " << max_wallclock_time << std::endl;
		std::cout << std::endl;
	}
};

}

#endif
