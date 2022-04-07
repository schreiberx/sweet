/*
 * SWE_Plane_TS_l_cn.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "SWE_Plane_TS_l_erk.hpp"

#include "../swe_plane_timeintegrators/SWE_Plane_TS_interface.hpp"
#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_erk.hpp"
#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_irk.hpp"


class SWE_Plane_TS_l_cn	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	double crank_nicolson_damping_factor = 0.5;
	//int timestepping_order_linear = 1;

	SWE_Plane_TS_l_erk ts_l_erk;
	SWE_Plane_TS_l_irk ts_l_irk;


private:
	void backward_euler_timestep_linear(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables
			double i_dt
	);



public:
	SWE_Plane_TS_l_cn(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			//int i_l_order,
			double i_crank_nicolson_damping_factor
	);

	void run_timestep(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Plane_TS_l_cn();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_CN_HPP_ */
