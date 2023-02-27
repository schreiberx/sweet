/*
 * SWE_Plane_TS_l_cn.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include "SWE_Plane_TS_l_erk.hpp"

#include "SWE_Plane_TS_interface.hpp"
#include "SWE_Plane_TS_l_erk.hpp"
#include "SWE_Plane_TS_l_irk.hpp"


class SWE_Plane_TS_l_cn	: public SWE_Plane_TS_interface
{
	sweet::ShackDictionary *shackDict;
	sweet::PlaneOperators &op;

	double crank_nicolson_damping_factor = 0.5;
	//int timestepping_order_linear = 1;

	SWE_Plane_TS_l_erk ts_l_erk;
	SWE_Plane_TS_l_irk ts_l_irk;


private:
	void backward_euler_timestep_linear(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables
			double i_dt
	);



public:
	SWE_Plane_TS_l_cn(
			sweet::ShackDictionary *io_shackDict,
			sweet::PlaneOperators &i_op
		);

	void setup(
			//int i_l_order,
			double i_crank_nicolson_damping_factor
	);

	void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Plane_TS_l_cn();
};

#endif
