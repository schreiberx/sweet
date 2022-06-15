/*
 * Burgers_Plane_TS_l_cn_n_sl.hpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_N_SL_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_N_SL_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include "../burgers_timeintegrators/Burgers_Plane_TS_interface.hpp"
#include "../burgers_timeintegrators/Burgers_Plane_TS_l_cn.hpp"



class Burgers_Plane_TS_l_cn_n_sl	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	PlaneData_Spectral u_prev, v_prev;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

	Burgers_Plane_TS_l_cn ts_l_cn;

public:
	Burgers_Plane_TS_l_cn_n_sl(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup();

	void run_timestep(
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables
			///PlaneData_Spectral &io_u_prev,	///< prognostic variables
			///PlaneData_Spectral &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


#if SWEET_PARAREAL
	void set_previous_solution(
				PlaneData_Spectral &i_u_prev,
				PlaneData_Spectral &i_v_prev
	) override
	{
		if (simVars.misc.verbosity > 5)
			std::cout << "set_previous_solution()" << std::endl;
		u_prev = i_u_prev;
		v_prev = i_v_prev;
	}
#endif



	virtual ~Burgers_Plane_TS_l_cn_n_sl();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_N_SL_HPP_ */
