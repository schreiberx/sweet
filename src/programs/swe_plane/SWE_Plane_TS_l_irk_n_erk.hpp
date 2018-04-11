/*
 * SWE_Plane_TS_l_irk_n_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_SWE_PLANE_TS_L_IRK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_SWE_PLANE_TS_L_IRK_N_ERK_HPP_

#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "SWE_Plane_TS_l_irk.hpp"

#include "../swe_plane/SWE_Plane_TS_interface.hpp"
#include "../swe_plane/SWE_Plane_TS_l_irk.hpp"


class SWE_Plane_TS_l_irk_n_erk	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order_linear;
	int timestepping_order_nonlinear;

	PlaneDataTimesteppingExplicitRK timestepping_rk;
	SWE_Plane_TS_l_irk ts_l_irk;


public:
	void euler_timestep_update_nonlinear(
			const PlaneData &i_h,	///< prognostic variables
			const PlaneData &i_u,	///< prognostic variables
			const PlaneData &i_v,	///< prognostic variables

			PlaneData &o_h_t,	///< time updates
			PlaneData &o_u_t,	///< time updates
			PlaneData &o_v_t	///< time updates
	);


public:
	SWE_Plane_TS_l_irk_n_erk(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_l_order,
			int i_n_order
	);

	void run_timestep(
			PlaneData &io_h,	///< prognostic variables
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);

        SWE_Plane_TS_l_irk& get_implicit_timestepper() {return ts_l_irk;}

	virtual ~SWE_Plane_TS_l_irk_n_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
