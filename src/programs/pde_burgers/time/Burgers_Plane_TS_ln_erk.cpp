/*
 * Burgers_Plane_TS_ln_erk.cpp
 *
 *  Created on: 17 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "../burgers_timeintegrators/Burgers_Plane_TS_ln_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void Burgers_Plane_TS_ln_erk::euler_timestep_update(
		const sweet::PlaneData_Spectral &i_tmp,	///< prognostic variables (perturbed part of height)
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_tmp_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (shackDict.misc.verbosity > 2)
		std::cout << "p_run_euler_timestep_update()" << std::endl;

	//TODO: staggering vs. non staggering

///#if SWEET_USE_PLANE_SPECTRAL_SPACE
   o_tmp_t.spectral_set_zero();
   o_u_t.spectral_set_zero();
   o_v_t.spectral_set_zero();
///#endif
///	o_tmp_t.physical_set_all(0);
///	o_u_t.physical_set_all(0);
///	o_v_t.physical_set_all(0);

	// u and v updates
	o_u_t = -(i_u*op.diff_c_x(i_u) + i_v*op.diff_c_y(i_u));
	o_u_t += shackDict.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));

	o_v_t = -(i_u*op.diff_c_x(i_v) + i_v*op.diff_c_y(i_v));
	o_v_t += shackDict.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));

	o_tmp_t.spectral_set_zero();
}



void Burgers_Plane_TS_ln_erk::run_timestep(
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_u_prev,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_v_prev,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Burgers_Plane_TS_ln_erk: Only constant time step size allowed");


	// setup dummy data
	sweet::PlaneData_Spectral tmp(io_u.planeDataConfig);
//#if SWEET_USE_PLANE_SPECTRAL_SPACE
	tmp.spectral_set_zero();
//#endif
//	tmp.physical_set_all(0);

	if (timestepping_order != 21)
	{
		// setup dummy data
		tmp.spectral_set_zero();

		// run standard Runge Kutta
		timestepping_rk.run_timestep(
			this,
			&Burgers_Plane_TS_ln_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			tmp, io_u, io_v,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
	}
	// Explicit Runge-Kutta with order 1 in diffusion and order 2 in advection
	else
	{
		if (shackDict.misc.verbosity > 2)
			std::cout << "run_timestep_erk()" << std::endl;

		sweet::PlaneData_Spectral u=io_u;
		sweet::PlaneData_Spectral v=io_v;

		double t = i_fixed_dt;
/*
		// Modify timestep to final time if necessary
		double t = i_fixed_dt;
		if (shackDict.timecontrol.current_simulation_time+i_fixed_dt < shackDict.timecontrol.max_simulation_time)
			t = i_fixed_dt;
		else
			t = shackDict.timecontrol.max_simulation_time-shackDict.timecontrol.current_simulation_time;
*/

		if (shackDict.disc.space_use_spectral_basis_diffs) //spectral
		{
			sweet::PlaneData_Spectral u1 = u + t*shackDict.sim.viscosity*(op.diff2_c_x(u)+op.diff2_c_y(u))
						   - 0.5*t*(u*op.diff_c_x(u)+v*op.diff_c_y(u));
			sweet::PlaneData_Spectral v1 = v + t*shackDict.sim.viscosity*(op.diff2_c_x(v)+op.diff2_c_y(v))
						   - 0.5*t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));

			io_u = u + t*shackDict.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
				  - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1));
			io_v = v + t*shackDict.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
				  - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));

		} else { //Jacobi
			SWEETError("NOT available");
		}
	}
}



/*
 * Setup
 */
void Burgers_Plane_TS_ln_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Burgers_Plane_TS_ln_erk::Burgers_Plane_TS_ln_erk(
		sweet::ShackDictionary &i_shackDict,
		PlaneOperators &i_op
)	:
		shackDict(i_shackDict),
		op(i_op)
{
	setup(shackDict.disc.timestepping_order);
}



Burgers_Plane_TS_ln_erk::~Burgers_Plane_TS_ln_erk()
{
}

