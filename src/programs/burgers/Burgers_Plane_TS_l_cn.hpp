/*
 * Burgers_Plane_TS_l_cn.hpp
 *
 *  Created on: 02 October 2017
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_HPP_

#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "Burgers_Plane_TS_interface.hpp"
#include <benchmarks_plane/BurgersValidationBenchmarks.hpp>



class Burgers_Plane_TS_l_cn	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;


public:
	Burgers_Plane_TS_l_cn(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup();

	void run_timestep(
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables
			PlaneData &io_u_prev,	///< prognostic variables
			PlaneData &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~Burgers_Plane_TS_l_cn();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_HPP_ */
