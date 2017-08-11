/*
 * Burgers_Plane_TS_ln_adomian.hpp
 *
 *  Created on: 03 August 2017
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_Burgers_PLANE_TS_LN_ADOMIAN_HPP_
#define SRC_PROGRAMS_Burgers_PLANE_TS_LN_ADOMIAN_HPP_

#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "Burgers_Plane_TS_interface.hpp"
#include <benchmarks_plane/BurgersValidationBenchmarks.hpp>



class Burgers_Plane_TS_ln_adomian	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;

public:
	Burgers_Plane_TS_ln_adomian(
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



	virtual ~Burgers_Plane_TS_ln_adomian();
};

#endif /* SRC_PROGRAMS_Burgers_PLANE_TS_LN_ADOMIAN_HPP_ */
