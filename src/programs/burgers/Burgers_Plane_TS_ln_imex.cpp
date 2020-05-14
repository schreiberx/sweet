/*
 * Burgers_Plane_TS_ln_imex.cpp
 *
 *  Created on: 17 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_ln_imex.hpp"



void Burgers_Plane_TS_ln_imex::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (simVars.misc.verbosity > 2)
		std::cout << "Burgers_Plane::run_timestep_imex()" << std::endl;

	PlaneData u=io_u;
	PlaneData v=io_v;
	double t = i_fixed_dt;

	// Setting explicit right hand side and operator of the left hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;

	if (timestepping_order == 1)
	{
		rhs_u += - t*(u*op.diff_c_x(u)+v*op.diff_c_y(u));
		rhs_v += - t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));
	}
	else if (timestepping_order == 2)
	{
		rhs_u += - 0.5*t*(u*op.diff_c_x(u)+v*op.diff_c_y(u));
		rhs_v += - 0.5*t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));
	}
	else
		SWEETError("The chosen timestepping-order is not possible with IMEX");

	if (simVars.disc.space_use_spectral_basis_diffs) //spectral
	{

		PlaneData lhs = u;
		if (timestepping_order == 1)
		{
			lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		}
		else
		{
			lhs = ((-t*0.5)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		}
		PlaneData u1 = rhs_u.spectral_div_element_wise(lhs);
		PlaneData v1 = rhs_v.spectral_div_element_wise(lhs);

		io_u = u + t*simVars.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
			  - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1));
		io_v = v + t*simVars.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
			  - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));

	} else { //Jacobi
		SWEETError("NOT available");
	}
}



/*
 * Setup
 */
void Burgers_Plane_TS_ln_imex::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Burgers_Plane_TS_ln_imex::Burgers_Plane_TS_ln_imex(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Burgers_Plane_TS_ln_imex::~Burgers_Plane_TS_ln_imex()
{
}

