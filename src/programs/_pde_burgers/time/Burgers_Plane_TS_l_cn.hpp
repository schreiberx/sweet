/*
 * Burgers_Plane_TS_l_cn.hpp
 *
 *  Created on: 02 October 2017
 * Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_L_CN_HPP_

#include <limits>
#include <sweet/core/plane/sweet::PlaneData_Spectral.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>

#include "../burgers_benchmarks/BurgersValidationBenchmarks.hpp"
#include "../burgers_timeintegrators/Burgers_Plane_TS_interface.hpp"



class Burgers_Plane_TS_l_cn	: public Burgers_Plane_TS_interface
{
	sweet::ShackDictionary &shackDict;
	PlaneOperators &op;


public:
	Burgers_Plane_TS_l_cn(
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



	virtual ~Burgers_Plane_TS_l_cn();
};

#endif
