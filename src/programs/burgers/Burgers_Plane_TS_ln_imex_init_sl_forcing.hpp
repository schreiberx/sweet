/*
 * Burgers_Plane_TS_ln_imex_init_sl_forcing.hpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_INIT_SL_FORCING_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_INIT_SL_FORCING_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "Burgers_Plane_TS_interface.hpp"

#include "Burgers_Plane_TS_ln_imex_forcing.hpp"
#include "Burgers_Plane_TS_l_irk_n_sl_forcing.hpp"



class Burgers_Plane_TS_ln_imex_init_sl_forcing	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	Burgers_Plane_TS_ln_imex_forcing ts_ln_imex_forcing;
	Burgers_Plane_TS_l_irk_n_sl_forcing ts_l_irk_n_sl_forcing;

public:
	Burgers_Plane_TS_ln_imex_init_sl_forcing(
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



	virtual ~Burgers_Plane_TS_ln_imex_init_sl_forcing();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_INIT_SL_FORCING_HPP_ */
