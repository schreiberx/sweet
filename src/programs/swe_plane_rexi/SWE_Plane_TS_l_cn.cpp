/*
 * SWE_Plane_TS_l_cn.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_cn.hpp"
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>








void SWE_Plane_TS_l_cn::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt <= 0)
		FatalError("SWE_Plane_TS_l_cn: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time - i_simulation_timestamp;


	PlaneData h_linear_t1 = io_h;
	PlaneData u_linear_t1 = io_u;
	PlaneData v_linear_t1 = io_v;

	double o_dummy;
	ts_l_erk.run_timestep(
			io_h,
			io_u,
			io_v,

			o_dummy,
			i_fixed_dt*(1.0-crank_nicolson_damping_factor),
			i_simulation_timestamp,
			i_max_simulation_time
		);

	ts_l_irk.run_timestep(
			io_h, io_u, io_v,
			o_dummy,
			i_fixed_dt*crank_nicolson_damping_factor,
			i_simulation_timestamp,
			i_max_simulation_time
		);

	o_dt = i_fixed_dt;
}



/*
 * Setup
 */
void SWE_Plane_TS_l_cn::setup(
		int i_l_order,
		double i_crank_nicolson_damping_factor
)
{
	timestepping_order_linear = i_l_order;

	if (timestepping_order_linear != 2)
		FatalError("SWE_Plane_TS_l_cn: Only 2nd order TS (Because of Crank Nicolson) supported with this implementation");

	crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;
}


SWE_Plane_TS_l_cn::SWE_Plane_TS_l_cn(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		ts_l_erk(simVars, op),
		ts_l_irk(simVars, op)
{
	ts_l_irk.setup(1);
	ts_l_erk.setup(1);
}



SWE_Plane_TS_l_cn::~SWE_Plane_TS_l_cn()
{
}

