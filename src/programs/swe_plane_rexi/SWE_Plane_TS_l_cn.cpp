/*
 * SWE_Plane_TS_l_cn.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 *
 */

#include "SWE_Plane_TS_l_cn.hpp"
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>

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

void SWE_Plane_TS_l_cn::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Plane_TS_l_cn: Only constant time step size allowed (Please set --dt)");


	PlaneData h_linear_t1 = io_h;
	PlaneData u_linear_t1 = io_u;
	PlaneData v_linear_t1 = io_v;

	ts_l_erk.run_timestep(
			io_h,
			io_u,
			io_v,
			i_dt*(1.0-crank_nicolson_damping_factor),
			i_simulation_timestamp
		);

	ts_l_irk.run_timestep(
			io_h, io_u, io_v,
			i_dt*crank_nicolson_damping_factor,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Plane_TS_l_cn::setup(
		//int i_l_order,
		double i_crank_nicolson_damping_factor
)
{
	//Force 1st order implicit and explicit schemes to achieve the 2nd order CN
	//timestepping_order_linear = 1; //i_l_order;

	//if (timestepping_order_linear != 1)
	//{
	//	std::cout << "SWE_Plane_TS_l_cn Warning: Using half explicit/implicit euler to achieve Crank-Nicolson" << std::endl;
	//   std::cout << "                             even though you wanted to have a " << timestepping_order_linear << " order explicit scheme." << std::endl;
	//	//FatalError("SWE_Plane_TS_l_cn: Only 2nd order TS (Because of Crank Nicolson) supported with this implementation");
	//}

	if (simVars.disc.use_staggering)
		FatalError("Staggering not supported for l_cn");

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
	//Force 1st order implicit and explicit schemes to achieve the 2nd order CN
	ts_l_irk.setup(1);
	ts_l_erk.setup(1);
}



SWE_Plane_TS_l_cn::~SWE_Plane_TS_l_cn()
{
}

