/*
 * ADV_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_ADV_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>


class Adv_Plane_TS_interface
{
public:
	virtual void run_timestep(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	) = 0;
};

#endif /* SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_LN_ERK_HPP_ */
