/*
 * SWE_Plane_TS_l_cn.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 *
 */

#include "SWE_Plane_TS_l_cn.hpp"



bool SWE_Plane_TS_l_cn::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	ts_l_erk.shackRegistration(io_shackDict);
	ts_l_irk.shackRegistration(io_shackDict);
	PDESWEPlaneTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}



/*
 * This is a Crank-Nicolson scheme for linear equation
 *
 * 1) Takes an explicit 1/2 step (or whatever controlled by crank_nicolson_damping_factor) of explicit euler
 * 2) Then, takes an implicit 1/2 step (or whatever controlled by crank_nicolson_damping_factor) with and implicit euler scheme
 *
 * With explicit Euler, may be viewed as classic CN:
 *
 * U_t = L U(n)
 *
 *
 * (U(n+1) - U(n)) / dt = 0.5*(L U(n+1) + L U(n))
 *
 * <=> U(n+1) - U(n) =  dt/2 *(L U(n+1) + L U(n))
 *
 * <=> (1-dt/2L) U(n+1)= (1 + dt/2 L) U(n)
 *     ---------------   -----------------
 *           |                  |
 *  (1/2 implicit euler)  (1/2 explicit euler)
 *
 * Comment added by P. Peixoto on 4 Sept 2017
 */

void SWE_Plane_TS_l_cn::runTimestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_cn: Only constant time step size allowed (Please set --dt)");


	sweet::PlaneData_Spectral h_linear_t1 = io_h;
	sweet::PlaneData_Spectral u_linear_t1 = io_u;
	sweet::PlaneData_Spectral v_linear_t1 = io_v;

	ts_l_erk.runTimestep(
			io_h,
			io_u,
			io_v,
			i_dt*(1.0-crank_nicolson_damping_factor),
			i_simulation_timestamp
		);

	ts_l_irk.runTimestep(
			io_h, io_u, io_v,
			i_dt*crank_nicolson_damping_factor,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
bool SWE_Plane_TS_l_cn::setup(
		sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_cn");

	crank_nicolson_damping_factor = 0.5;

	ts_l_irk.setup(io_ops, 1);
	ts_l_erk.setup(io_ops, 1);

	return true;
}
