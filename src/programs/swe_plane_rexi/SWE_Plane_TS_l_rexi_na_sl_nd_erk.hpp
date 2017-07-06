/*
 * SWE_Plane_TS_l_rexi_ns_sl_nd_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_ERK_HPP_

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



class SWE_Plane_TS_l_rexi_na_sl_nd_erk	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	SWE_Plane_TS_l_rexi ts_l_rexi;

	int with_nonlinear;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	PlaneData h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

public:
	SWE_Plane_TS_l_rexi_na_sl_nd_erk(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			double i_h,						///< sampling size
			int i_M,						///< number of sampling points
			int i_L,						///< number of sampling points for Gaussian approximation
											///< set to 0 for auto detection
			bool i_rexi_half,				///< use half-pole reduction
			bool i_rexi_normalization,		///< REXI normalization

			int i_with_nonlinear
	);

	void run_timestep(
			PlaneData &io_h,	///< prognostic variables
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	);



	virtual ~SWE_Plane_TS_l_rexi_na_sl_nd_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
