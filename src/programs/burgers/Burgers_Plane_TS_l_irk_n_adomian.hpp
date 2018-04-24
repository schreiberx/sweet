/*
 * Burgers_Plane_TS_l_irk_n_adomian.hpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_N_ADOMIAN_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_N_ADOMIAN_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "Burgers_Plane_TS_interface.hpp"

#include "Burgers_Plane_TS_l_irk.hpp"
#include "Burgers_Plane_TS_n_adomian.hpp"



class Burgers_Plane_TS_l_irk_n_adomian	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	Burgers_Plane_TS_l_irk ts_l_irk;
	Burgers_Plane_TS_n_adomian ts_n_adomian;

public:
	Burgers_Plane_TS_l_irk_n_adomian(
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



	virtual ~Burgers_Plane_TS_l_irk_n_adomian();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_N_ADOMIAN_HPP_ */
