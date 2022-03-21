/*
 * Burgers_Plane_TS_ln_cole_hopf.hpp
 *
 *  Created on: 08 August 2017
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_LN_COLE_HOPF_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_LN_COLE_HOPF_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#include "Burgers_Plane_TS_l_direct.hpp"
#include "Burgers_Plane_TS_interface.hpp"
#include "../burgers_benchmarks/BurgersValidationBenchmarks.hpp"



class Burgers_Plane_TS_ln_cole_hopf	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;

public:
	Burgers_Plane_TS_l_direct ts_l_direct;

public:
	Burgers_Plane_TS_ln_cole_hopf(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup();

	void run_timestep(
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables
			PlaneData_Spectral &io_u_prev,	///< prognostic variables
			PlaneData_Spectral &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~Burgers_Plane_TS_ln_cole_hopf();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_LN_COLE_HOPF_HPP_ */
