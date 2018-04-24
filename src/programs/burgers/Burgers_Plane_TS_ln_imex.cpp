/*
 * Burgers_Plane_TS_ln_imex.cpp
 *
 *  Created on: 17 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_ln_imex.hpp"


/*
 * Implementation of the IMEX method applied to Burgers' equation according to the algorithm provided in
 * Ascher et al. (1997) Implicit-Explicit Runge-Kutta Methods for Time-Dependent Partial Differential Equations
 * for a timestep from u^n to u^{n+1}
 *
 * First order implementation is IMEX(1,1,1)
 *     (The IMEX(1,2,1) implementation is available as-well, see below)
 * Second order implementation is IMEX(1,2,2)
 */
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

	// Use local variables for calculation
	PlaneData u=io_u;
	PlaneData v=io_v;
	double dt = i_fixed_dt;

	// Initialize variables for right-hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;
#if 1
	// Calculate \hat{k}_1 = f(u^n)
	PlaneData explK_u = -(u*op.diff_c_x(u)+v*op.diff_c_y(u));
	PlaneData explK_v = -(u*op.diff_c_x(v)+v*op.diff_c_y(v));

	// Calculate k_1
	if (timestepping_order == 1)
	{
		rhs_u = simVars.sim.viscosity*(op.diff2_c_x(u+dt*explK_u)+op.diff2_c_y(u+dt*explK_u));
		rhs_v = simVars.sim.viscosity*(op.diff2_c_x(v+dt*explK_v)+op.diff2_c_y(v+dt*explK_v));
	}
	else if (timestepping_order == 2)
	{
		rhs_u = simVars.sim.viscosity*(op.diff2_c_x(u+dt/2*explK_u)+op.diff2_c_y(u+dt/2*explK_u));
		rhs_v = simVars.sim.viscosity*(op.diff2_c_x(v+dt/2*explK_v)+op.diff2_c_y(v+dt/2*explK_v));
	}
	else
		FatalError("The chosen timestepping-order is not possible with IMEX");

	if (simVars.disc.use_spectral_basis_diffs) //spectral
	{

		PlaneData lhs = u;
		if (timestepping_order == 1)
		{
			lhs = ((-dt)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		}
		else
		{
			lhs = ((-dt*0.5)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		}
		PlaneData implK_u = rhs_u.spectral_div_element_wise(lhs);
		PlaneData implK_v = rhs_v.spectral_div_element_wise(lhs);

		// Calculate \hat{k}_2
		if (timestepping_order == 2)
		{
			explK_u = - ((u+dt/2*explK_u+dt/2*implK_u)*op.diff_c_x(u+dt/2*explK_u+dt/2*implK_u)
					+ (v+dt/2*explK_v+dt/2*implK_v)*op.diff_c_y(u+dt/2*explK_u+dt/2*implK_u));
			explK_v = - ((u+dt/2*explK_u+dt/2*implK_u)*op.diff_c_x(v+dt/2*explK_v+dt/2*implK_v)
					+ (v+dt/2*explK_v+dt/2*implK_v)*op.diff_c_y(v+dt/2*explK_v+dt/2*implK_v));
		}
#if 0 //used for Ascher(1,2,1) from Ascher et al. (1997)
		else
		{
			explK_u = - ((u+dt*explK_u+dt*implK_u)*op.diff_c_x(u+dt*explK_u+dt*implK_u)
					+ (v+dt*explK_v+dt*implK_v)*op.diff_c_y(u+dt*explK_u+dt*implK_u));
			explK_v = - ((u+dt*explK_u+dt*implK_u)*op.diff_c_x(v+dt*explK_v+dt*implK_v)
					+ (v+dt*explK_v+dt*implK_v)*op.diff_c_y(v+dt*explK_v+dt*implK_v));
		}
#endif

		// Calculate u^{n+1} = u^n + dt*k_1 + dt*\hat{k}_2
		io_u = u + dt*explK_u + dt*implK_u;
		io_v = v + dt*explK_v + dt*implK_v;

	} else { //Jacobi
		FatalError("NOT available");
}
#else //Wrong Implementation
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
		FatalError("The chosen timestepping-order is not possible with IMEX");

	if (simVars.disc.use_spectral_basis_diffs) //spectral
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
		FatalError("NOT available");
	}
#endif
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

