/*
 * PararealController.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREALCONTROLLER_SERIAL_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREALCONTROLLER_SERIAL_HPP_

#include <parareal/PararealSimulation_Base.hpp>

/**
 * This class takes over the control and
 * calls methods offered via PararealSimulation.
 */
class PararealController_Serial
{
	PararealSimulation_Base *base = nullptr;
	int coarse_time_steps = -1;
	double timeframe_start = -1;
	double timeframe_end = -1;

	void set_simulation_base_class(
			PararealSimulation_Base *i_base
	)
	{
		base = i_base;
	}

	void set_coarse_time_steps(
			int i_coarse_time_steps
	)
	{
		coarse_time_steps = i_coarse_time_steps;
	}

	void set_timeframe(
			double i_start,		///< start of total simulation time
			double i_end		///< end of total simulation time
	)
	{
		timeframe_start = i_start;
		timeframe_end = i_end;
	}

	void run(
	)
	{

	}
};





#endif /* SRC_INCLUDE_PARAREAL_PARAREALCONTROLLER_SERIAL_HPP_ */
