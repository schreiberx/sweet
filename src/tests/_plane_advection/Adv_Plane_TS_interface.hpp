/*
 * ADV_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_ADV_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/core/SimulationVariables.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>


class Adv_Plane_TS_interface
{
public:
	virtual void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp
	) = 0;
};

#endif /* SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_LN_ERK_HPP_ */
