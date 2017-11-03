/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_settls.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include "SWE_Plane_TS_interface.hpp"
#include "SWE_Plane_TS_l_direct.hpp"
#include "SWE_Plane_TS_l_rexi.hpp"



class SWE_Plane_TS_l_rexi_na_sl_nd_settls	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	SWE_Plane_TS_l_rexi ts_l_rexi;

	int with_linear_div_only;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	//Previous values (t_n-1)
	PlaneData h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

public:
	SWE_Plane_TS_l_rexi_na_sl_nd_settls(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			REXI_SimulationVariables &i_rexi,

			int i_with_nonlinear
	);

	void run_timestep(
			PlaneData &io_h,	///< prognostic variables
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Plane_TS_l_rexi_na_sl_nd_settls();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_ */
