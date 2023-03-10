/*
 * Burgers_Plane_TS_l_cn.cpp
 *
 *  Created on: 02 October 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_l_cn.hpp"


void Burgers_Plane_TS_l_cn::runTimestep(
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_u_prev,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_v_prev,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Burgers_Plane_TS_l_cn: Only constant time step size allowed");

	double dt = i_fixed_dt;

	sweet::PlaneData_Spectral rhs_u = io_u + 0.5*dt*shackDict.sim.viscosity*(op.diff2_c_x(io_u) + op.diff2_c_y(io_u));
	sweet::PlaneData_Spectral rhs_v = io_v + 0.5*dt*shackDict.sim.viscosity*(op.diff2_c_x(io_v) + op.diff2_c_y(io_v));
	sweet::PlaneData_Spectral lhs = ((-0.5*dt)*shackDict.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);

	io_u = rhs_u.spectral_div_element_wise(lhs);
	io_v = rhs_v.spectral_div_element_wise(lhs);

}



/*
 * Setup
 */
void Burgers_Plane_TS_l_cn::setup()
{
}


Burgers_Plane_TS_l_cn::Burgers_Plane_TS_l_cn(
		sweet::ShackDictionary &i_shackDict,
		PlaneOperators &i_op
)	:
		shackDict(i_shackDict),
		op(i_op)
{
}



Burgers_Plane_TS_l_cn::~Burgers_Plane_TS_l_cn()
{
}

