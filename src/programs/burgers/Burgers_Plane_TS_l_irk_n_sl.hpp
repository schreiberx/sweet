/*
 * Burgers_Plane_TS_l_irk_n_sl.hpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_Burgers_PLANE_TS_L_IRK_N_SL_HPP_
#define SRC_PROGRAMS_Burgers_PLANE_TS_L_IRK_N_SL_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include "Burgers_Plane_TS_interface.hpp"
#include <benchmarks_plane/BurgersValidationBenchmarks.hpp>



class Burgers_Plane_TS_l_irk_n_sl	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

public:
	Burgers_Plane_TS_l_irk_n_sl(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup();

	void run_timestep(
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables
			PlaneData &io_u_prev,	///< prognostic variables
			PlaneData &io_v_prev,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	);



	virtual ~Burgers_Plane_TS_l_irk_n_sl();
};

#endif /* SRC_PROGRAMS_Burgers_PLANE_TS_L_IRK_N_SL_HPP_ */
