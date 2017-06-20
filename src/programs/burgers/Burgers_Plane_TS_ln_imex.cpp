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

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (simVars.misc.verbosity > 2)
		std::cout << "Burgers_Plane::run_timestep_imex()" << std::endl;

	PlaneData u=io_u;
	PlaneData v=io_v;

	// Modify timestep to final time if necessary
	double& t = o_dt;
	if (simVars.timecontrol.current_simulation_time+i_fixed_dt < i_max_simulation_time)
		t = i_fixed_dt;
	else
		t = i_max_simulation_time-simVars.timecontrol.current_simulation_time;

	// Initialize and set timestep dependent source for manufactured solution
	PlaneData f(io_u.planeDataConfig);
	PlaneData ff(io_u.planeDataConfig);
	BurgersValidationBenchmarks::set_source(simVars.timecontrol.current_simulation_time,simVars,simVars.disc.use_staggering,f);
	BurgersValidationBenchmarks::set_source(simVars.timecontrol.current_simulation_time+0.5*t,simVars,simVars.disc.use_staggering,ff);
	f.request_data_spectral();
	ff.request_data_spectral();

	// Setting explicit right hand side and operator of the left hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;

	rhs_u += - 0.5*t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + 0.5*t*f;
	rhs_v += - 0.5*t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));

	if (simVars.disc.use_spectral_basis_diffs) //spectral
	{

		PlaneData lhs = u;
		lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
        PlaneData u1 = rhs_u.spectral_div_element_wise(lhs);
        PlaneData v1 = rhs_v.spectral_div_element_wise(lhs);

        io_u = u + t*simVars.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
              - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1)) +ff*t;
        io_v = v + t*simVars.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
              - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));

	} else { //Jacobi
		FatalError("NOT available");
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

