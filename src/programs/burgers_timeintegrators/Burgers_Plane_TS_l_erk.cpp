/*
 * Burgers_Plane_TS_l_erk.cpp
 *
 *  Created on: 17 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "../burgers_timeintegrators/Burgers_Plane_TS_l_erk.hpp"


/*
 * Main routine for method to be used in case of finite differences
 */
void Burgers_Plane_TS_l_erk::euler_timestep_update(
		const PlaneData_Spectral &i_tmp,	///< prognostic variables (perturbed part of height)
		const PlaneData_Spectral &i_u,	///< prognostic variables
		const PlaneData_Spectral &i_v,	///< prognostic variables

		PlaneData_Spectral &o_tmp_t,	///< time updates
		PlaneData_Spectral &o_u_t,	///< time updates
		PlaneData_Spectral &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (simVars.misc.verbosity > 2)
		std::cout << "p_run_euler_timestep_update()" << std::endl;

	//TODO: staggering vs. non staggering

///#if SWEET_USE_PLANE_SPECTRAL_SPACE
   o_tmp_t.spectral_set_zero();
   o_u_t.spectral_set_zero();
   o_v_t.spectral_set_zero();
////#endif
////   o_tmp_t.physical_set_all(0);
////   o_u_t.physical_set_all(0);
////   o_v_t.physical_set_all(0);

	// u and v updates
	o_u_t = simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));

	o_v_t = simVars.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));

	o_tmp_t.spectral_set_zero();
}



void Burgers_Plane_TS_l_erk::run_timestep(
		PlaneData_Spectral &io_u,	///< prognostic variables
		PlaneData_Spectral &io_v,	///< prognostic variables
		///PlaneData_Spectral &io_u_prev,	///< prognostic variables
		///PlaneData_Spectral &io_v_prev,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Burgers_Plane_TS_l_erk: Only constant time step size allowed");


	// setup dummy data
	PlaneData_Spectral tmp(io_u.planeDataConfig);
//#if SWEET_USE_PLANE_SPECTRAL_SPACE
	tmp.spectral_set_zero();
//#endif
//	tmp.physical_set_all(0);

	// run standard Runge Kutta
	timestepping_rk.run_timestep(
		this,
		&Burgers_Plane_TS_l_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
		tmp, io_u, io_v,
		i_fixed_dt,
		timestepping_order,
		i_simulation_timestamp
	);

}



/*
 * Setup
 */
void Burgers_Plane_TS_l_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Burgers_Plane_TS_l_erk::Burgers_Plane_TS_l_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Burgers_Plane_TS_l_erk::~Burgers_Plane_TS_l_erk()
{
}

