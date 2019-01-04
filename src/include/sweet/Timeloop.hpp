/*
 * Timeloop.hpp
 *
 *  Created on: Jan 3, 2019
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_TIMELOOP_HPP_
#define SRC_INCLUDE_SWEET_TIMELOOP_HPP_



#define SWEET_TIMELOOP		\
	for (	simVars.timecontrol.current_simulation_time = 0;	\
			simVars.timecontrol.current_simulation_time < simVars.timecontrol.max_simulation_time*(1.0-1e-12);	\
			simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size	\
	)



#endif /* SRC_INCLUDE_SWEET_TIMELOOP_HPP_ */
