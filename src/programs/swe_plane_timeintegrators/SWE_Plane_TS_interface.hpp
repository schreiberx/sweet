/*
 * SWE_Plane_TS_interface.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>


class SWE_Plane_TS_interface
{
public:
	virtual void run_timestep(
			PlaneData_Spectral &io_h_pert,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,		///< time step size
			double i_sim_timestamp
	) = 0;

	// for parareal SL
	virtual void set_previous_solution(
				PlaneData &i_h_prev,
				PlaneData &i_u_prev,
				PlaneData &i_v_prev
	)
	{
	};

};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_INTERFACE_HPP_ */
