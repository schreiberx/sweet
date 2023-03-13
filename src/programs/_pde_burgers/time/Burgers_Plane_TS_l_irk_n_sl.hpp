/*
 * Burgers_Plane_TS_l_irk_n_sl.hpp
 *
 *  Created on: 14 June 2017
 * Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_N_SL_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_L_IRK_N_SL_HPP_

#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/sweet::PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include "../burgers_timeintegrators/Burgers_Plane_TS_interface.hpp"
#include "../burgers_timeintegrators/Burgers_Plane_TS_l_irk.hpp"



class Burgers_Plane_TS_l_irk_n_sl	: public Burgers_Plane_TS_interface
{
	sweet::ShackDictionary &shackDict;
	PlaneOperators &op;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

	Burgers_Plane_TS_l_irk ts_l_irk;

	sweet::PlaneData_Spectral u_prev, v_prev;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;


public:
	Burgers_Plane_TS_l_irk_n_sl(
			sweet::ShackDictionary &i_shackDict,
			PlaneOperators &i_op
		);

	void setup();

	void runTimestep(
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables
			///sweet::PlaneData_Spectral &io_u_prev,	///< prognostic variables
			///sweet::PlaneData_Spectral &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

#if SWEET_PARAREAL || SWEET_XBRAID
	void set_previous_solution(
				sweet::PlaneData_Spectral &i_u_prev,
				sweet::PlaneData_Spectral &i_v_prev
	) override
	{
		if (shackDict.misc.verbosity > 5)
			std::cout << "set_previous_solution()" << std::endl;
		u_prev = i_u_prev;
		v_prev = i_v_prev;
	}
#endif


	virtual ~Burgers_Plane_TS_l_irk_n_sl();
};

#endif
