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
#include <sys/time.h>


/**
 * \brief start, stop, continue and restart a virtual stopwatch
 */
class Stopwatch
{
	/**
	 * some storage for the time values at start of stopwatch and at stop of stopwatch
	 */
private:
	struct timeval timevalue_start;	///< time value of last start
	struct timeval timevalue_stop;	///< time value of last stop

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
			gettimeofday(&timevalue_start, NULL);

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
			gettimeofday(&timevalue_stop, NULL);

			time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
			suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

			time += (double)dsec + (double)dsusec/1000000.0;
		}
	}


	/**
	 * return the time in seconds
	 */
	inline double getTimeSinceStart()
	{
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		return time + (double)dsec + (double)dsusec/1000000.0;
	}


	/**
	 * return the current time value in seconds
	 *
	 * This has nothing to do with a stop watch and is only a convenient feature to access the current time value
	 */
	inline static double getCurrentClockSeconds()
	{
		struct timeval timevalue;		///< time value of last stop
		gettimeofday(&timevalue, NULL);
		return (double)timevalue.tv_sec + (double)timevalue.tv_usec/1000000.0f;
	}

	/**
	 * return the time in seconds
	 */
	inline double operator()()
	{
#if SWEET_DEBUG
		assert(recursive_counter == 0);
#endif

		return time;
	}
};

#endif
