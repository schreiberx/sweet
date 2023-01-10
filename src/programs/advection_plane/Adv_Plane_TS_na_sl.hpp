/*
 * Adv_Plane_TS_na_sl.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_SL_HPP_
#define SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_SL_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>

#include "Adv_Plane_TS_interface.hpp"


class Adv_Plane_TS_na_sl	: public Adv_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;

	PlaneDataSampler sampler2D;
	PlaneDataSemiLagrangian semiLagrangian;

	PlaneData_Spectral prog_u_prev, prog_v_prev;


	/**
	 * Position in physical space given in longitude/latitude angles
	 */
	ScalarDataArray posx_a, posy_a;

public:
	Adv_Plane_TS_na_sl(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			PlaneData_Spectral &io_phi,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~Adv_Plane_TS_na_sl();
};

#endif
