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


#ifndef CSTOPWATCH_HPP
#define CSTOPWATCH_HPP

#define OS_UNIX	1
//#include "lib/os_preprocessor_defines.h"
#include <cstddef>

#if OS_UNIX
extern "C"
{
	#include <sys/time.h>
}
#endif

#if OS_WINDOWS
extern "C"
{
	#include <windows.h>
	#include <mmsystem.h>
}
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
#if OS_WINDOWS
	DWORD timevalue_start;
	DWORD timevalue_stop;
#endif

#if OS_UNIX
	struct timeval timevalue_start;	///< time value of last start
	struct timeval timevalue_stop;	///< time value of last stop
#endif

public:
	double time;		///< stopped time (difference between start and stop time)

	/**
	 * reset time with 0
	 */
	void reset()
	{
		time = 0.0f;
	}

	/**
	 * initialize the stopwatch
	 */
	Stopwatch()
	{
		reset();
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
#if OS_UNIX
		gettimeofday(&timevalue_start, NULL);
#endif

#if OS_WINDOWS
		timevalue_start = timeGetTime();
#endif
	}

	/**
	 * stop the stop watch
	 *
	 * add the delta value to the 'time' counter
	 */
	inline void stop()
	{
#if OS_UNIX
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		time += (double)dsec + (double)dsusec/1000000.0;
#endif

#if OS_WINDOWS
		DWORD timevalue = timeGetTime();
		time += (double)timevalue/1000.0;
#endif
	}


	/**
	 * return the time in seconds
	 */
	inline double getTimeSinceStart()
	{
#if OS_UNIX
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		return time + (double)dsec + (double)dsusec/1000000.0;
#endif

#if OS_WINDOWS
		DWORD timestamp = timeGetTime();
		return time + (double)timestamp/1000.0;
#endif
	}


	/**
	 * return the current time value in seconds
	 *
	 * This has nothing to do with a stop watch and is only a convenient feature to access the current time value
	 */
	inline static double getCurrentClockSeconds()
	{
#if OS_UNIX
		struct timeval timevalue;		///< time value of last stop
		gettimeofday(&timevalue, NULL);
		return (double)timevalue.tv_sec + (double)timevalue.tv_usec/1000000.0f;
#endif

#if OS_WINDOWS
		DWORD timestamp = timeGetTime();
		return (double)timestamp/1000.0;
#endif
	}

	/**
	 * return the time in seconds
	 */
	inline double operator()()
	{
		return time;
	}
};

#endif
