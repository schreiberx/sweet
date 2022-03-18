/*
 * Burgers_Plane_TS_l_erk.hpp
 *
 *  Created on: 17 June 2017
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_L_ERK_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_L_ERK_HPP_


#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include "../burgers_timeintegrators/Burgers_Plane_TS_interface.hpp"

#include "Burgers_Plane_TS_interface.hpp"
#include "../burgers_benchmarks/BurgersValidationBenchmarks.hpp"


class Burgers_Plane_TS_l_erk	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;
	PlaneDataTimesteppingExplicitRK timestepping_rk;

private:
	void euler_timestep_update(
			const PlaneData_Spectral &i_tmp,	///< prognostic variables
			const PlaneData_Spectral &i_u,	///< prognostic variables
			const PlaneData_Spectral &i_v,	///< prognostic variables

			PlaneData_Spectral &o_tmp_t,	///< time updates
			PlaneData_Spectral &o_u_t,	///< time updates
			PlaneData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	Burgers_Plane_TS_l_erk(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables
			PlaneData_Spectral &io_u_prev,	///< prognostic variables
			PlaneData &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~Burgers_Plane_TS_l_erk();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_L_ERK_HPP_ */
