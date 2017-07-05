/*
 * Burgers_Plane_TS_ln_erk.hpp
 *
 *  Created on: 15 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_Burgers_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_Burgers_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>


class Burgers_Plane_TS_interface
{
public:
	virtual void run_timestep(
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables
			PlaneData &io_u_prev,	///< prognostic variables
			PlaneData &io_v_prev,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	) = 0;
};

#endif /* SRC_PROGRAMS_Burgers_PLANE_TS_LN_ERK_HPP_ */
