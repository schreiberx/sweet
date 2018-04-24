/*
 * Burgers_Plane_TS_l_irk_mms.hpp
 *
 *  Created on: 16 April 2018
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_MMS_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_MMS_HPP_

#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "Burgers_Plane_TS_interface.hpp"
#include <benchmarks_plane/BurgersValidationBenchmarks.hpp>



class Burgers_Plane_TS_l_irk_mms	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;
	double **table;
	bool table_crated;
	bool second_time;

public:
	Burgers_Plane_TS_l_irk_mms(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables
			PlaneData &io_u_prev,	///< prognostic variables
			PlaneData &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);

	void return_initial(
			PlaneData &init
	);

	virtual ~Burgers_Plane_TS_l_irk_mms();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_MMS_HPP_ */
