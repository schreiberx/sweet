/*
 * Burgers_Plane_TS_ln_imex_forcing.cpp
 *
 *  Created on: 17 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_ln_imex_forcing.hpp"



void Burgers_Plane_TS_ln_imex_forcing::run_timestep(
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

	//std::cout << std::endl << std::endl << "rhs_u" << std::endl;
	//rhs_u.print_physicalArrayData();

	if (simVars.disc.use_spectral_basis_diffs) //spectral
	{

		PlaneData lhs = u;
		lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
        PlaneData u1 = rhs_u.spectral_div_element_wise(lhs);
        PlaneData v1 = rhs_v.spectral_div_element_wise(lhs);

        //std::cout << std::endl << std::endl << "u1" << std::endl;
        //u1.print_physicalArrayData();

        io_u = u + t*simVars.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
              - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1)) +ff*t;
        io_v = v + t*simVars.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
              - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));

        //std::cout << std::endl << std::endl << "io_u" << std::endl;
        //io_u.print_physicalArrayData();

	} else { //Jacobi
		FatalError("NOT available");
	}
}



/*
 * Setup
 */
void Burgers_Plane_TS_ln_imex_forcing::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Burgers_Plane_TS_ln_imex_forcing::Burgers_Plane_TS_ln_imex_forcing(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Burgers_Plane_TS_ln_imex_forcing::~Burgers_Plane_TS_ln_imex_forcing()
{
}

