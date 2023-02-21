/*
 * ShackTimestepControl.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKTIMESTEPCONTROL_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKTIMESTEPCONTROL_HPP_


#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include "../shacks/ShackInterface.hpp"


/**
 * Timestepping
 */
class TimestepControl	:
		public sweet::ClassDictionaryInterface
{
public:
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
		printClass();
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

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "" << std::endl;
		std::cout << "Timecontrol:" << std::endl;
		std::cout << "	--dt [time]	timestep size, default=?" << std::endl;
		std::cout << "	--max-wallclock-time [time]	wallclock time limitation, default=-1" << std::endl;
		std::cout << "	-t [time]	maximum simulation time, default=-1 (infinity)" << std::endl;
		std::cout << "	-T [stepnr]	maximum number of time steps, default=-1 (infinity)" << std::endl;
		std::cout << "	-o [time]	time interval at which output should be written, (set to 0 for output at every time step), default=-1 (no output) " << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--dt", current_timestep_size);
		i_pa.getArgumentValueBy2Keys("--max-wallclock-time", "-t", max_wallclock_time);
		i_pa.getArgumentValueByKey("--T", max_timesteps_nr);

		setup_timestep_size = current_timestep_size;

		if (max_simulation_time < 0)
			return error.set("timecontrol.max_simulation_time < 0");

		if (max_timesteps_nr < 0)
			return error.set("timecontrol.max_timesteps_nr < 0");

		current_timestep_nr = 0;
		current_simulation_time = 0;
		current_timestep_size = setup_timestep_size;

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
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
};




#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKTIMESTEPCONTROL_HPP_ */
