/*
 * SWE_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereOperators.hpp>


class SWE_Sphere_TS_interface
{
public:
	virtual void run_timestep(
			SphereDataSpectral &io_h,	///< prognostic variables
			SphereDataSpectral &io_u,	///< prognostic variables
			SphereDataSpectral &io_v,	///< prognostic variables

			double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp
	) = 0;
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
