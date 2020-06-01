/*
 * Burgers_Plane_TS_ln_erk_forcing.hpp
 *
 *  Created on: 17 June 2017
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_LN_ERK_FORCING_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_LN_ERK_FORCING_HPP_


#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include "Burgers_Plane_TS_interface.hpp"
#include "../burgers_benchmarks/BurgersValidationBenchmarks.hpp"


class Burgers_Plane_TS_ln_erk_forcing	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;
	PlaneDataTimesteppingExplicitRK timestepping_rk;

private:
	void euler_timestep_update(
			const PlaneData &i_tmp,	///< prognostic variables
			const PlaneData &i_u,	///< prognostic variables
			const PlaneData &i_v,	///< prognostic variables

			PlaneData &o_tmp_t,	///< time updates
			PlaneData &o_u_t,	///< time updates
			PlaneData &o_v_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	Burgers_Plane_TS_ln_erk_forcing(
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



	virtual ~Burgers_Plane_TS_ln_erk_forcing();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_LN_ERK_FORCING_HPP_ */
