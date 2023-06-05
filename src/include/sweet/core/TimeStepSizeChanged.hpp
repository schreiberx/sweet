/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_TIMESTEPSIZECHANGED_HPP_
#define SRC_INCLUDE_SWEET_TIMESTEPSIZECHANGED_HPP_

#include <algorithm>
#include <iostream>


class TimeStepSizeChanged
{
	/*
	 * Determine whether a time step size has been sufficiently changed
	 * so that a new configuration of e.g. implicit solvers is required.
	 *
	 * This can happen due to round-off errors. Hence, with 52 bit mantissa/fraction,
	 * this is of  2^52 = 4,50359×10^15 accuracy. We choose 10^14 as a threshold.
	 */
public:
	static
	bool is_changed(
		double i_old_timestepSize,		///< old timestep size with which the timestepper has been configured
		double i_new_timestepSize,		///< new timestep size
		bool i_output = true
	)
	{
		if (std::abs(i_old_timestepSize - i_new_timestepSize)/std::max(i_old_timestepSize, i_new_timestepSize) > 1e-14)
		{
//#if SWEET_DEBUG
			if (i_output)
				std::cout << "Warning: Reducing time step size from " << i_old_timestepSize << " to " << i_new_timestepSize << std::endl;
//#endif

			return true;
		}

		return false;
	}

};



#endif
