/*
 * Copyright 2010 Martin Schreiber
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */



#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

#include <cstddef>
#include <cassert>
#include <iostream>

#define SWEET_TIMER_CHRONO	1

#if SWEET_TIMER_CHRONO
#	include <chrono>
#else
#	include <sys/time.h>
#endif

/**
 * \brief start, stop, continue and restart a virtual stopwatch
 */
class Stopwatch
{
	/**
	 * some storage for the time values at start of stopwatch and at stop of stopwatch
	 */
private:

#if SWEET_TIMER_CHRONO
	std::chrono::time_point<std::chrono::system_clock> timevalue_start;	///< time value of last start
	std::chrono::time_point<std::chrono::system_clock> timevalue_stop;	///< time value of last stop
#else

	struct timeval timevalue_start;	///< time value of last start
	struct timeval timevalue_stop;	///< time value of last stop
#endif

	int recursive_counter;	/// count recursions to support nested calls

public:
	double time;		///< stopped time (difference between start and stop time)

	/**
	 * reset time with 0
	 */
	void reset()
	{
		time = 0.0f;

		recursive_counter = 0;
	}

	/**
	 * initialize the stopwatch
	 */
	Stopwatch(bool i_start = false)
	{
		reset();

		if (i_start)
			start();
	}


	/**
	 * start the stop watch:
	 *
	 * this function only reads out the current time value and stored it to the stopwatch
	 * class for a later computation of the delta time between start() and stop() functions
	 */
	inline void start()
	{
		if (recursive_counter == 0)
		{
#if SWEET_TIMER_CHRONO
			timevalue_start = std::chrono::system_clock::now();
#else
			gettimeofday(&timevalue_start, NULL);
#endif
		}

		recursive_counter++;
	}


	/**
	 * stop the stop watch
	 *
	 * add the delta value to the 'time' counter
	 */
	inline void stop()
	{
		recursive_counter--;

		if (recursive_counter == 0)
		{
#if SWEET_TIMER_CHRONO
			timevalue_stop = std::chrono::system_clock::now();

			time += ((std::chrono::duration<double>)(timevalue_stop-timevalue_start)).count();
#else
			gettimeofday(&timevalue_stop, NULL);

			time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
			suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

			time += (double)dsec + (double)dsusec/1000000.0;
#endif
		}
	}


	/**
	 * return the time in seconds
	 */
	inline double getTimeSinceStart()
	{
#if SWEET_TIMER_CHRONO
		timevalue_stop = std::chrono::system_clock::now();

		return ((std::chrono::duration<double>)(timevalue_stop-timevalue_start)).count();
#else
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		return time + (double)dsec + (double)dsusec/1000000.0;
#endif
	}


	/**
	 * return the time in seconds
	 */
	inline double operator()()
	{
#if SWEET_DEBUG
		if (recursive_counter != 0)
		{
			std::cout << "Recursion counter should be 0" << std::endl;
			assert(false);
		}
#endif

		return time;
	}
};

#endif
