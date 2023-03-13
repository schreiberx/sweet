/*
 * Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "../burgers_timeintegrators/Burgers_Plane_TS_ln_imex_forcing.hpp"



void Burgers_Plane_TS_ln_imex_forcing::runTimestep(
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_u_prev,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_v_prev,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (shackDict.misc.verbosity > 2)
		std::cout << "Burgers_Plane::run_timestep_imex()" << std::endl;

	sweet::PlaneData_Spectral u=io_u;
	sweet::PlaneData_Spectral v=io_v;
	double dt = i_fixed_dt;

	// Initialize and set timestep dependent source for manufactured solution
	sweet::PlaneData_Spectral f(io_u.planeDataConfig);
	sweet::PlaneData_Spectral ff(io_u.planeDataConfig);
	BurgersValidationBenchmarks::set_source(shackDict.timecontrol.current_simulation_time,shackDict,shackDict.disc.space_grid_use_c_staggering,f);
	BurgersValidationBenchmarks::set_source(shackDict.timecontrol.current_simulation_time+0.5*dt,shackDict,shackDict.disc.space_grid_use_c_staggering,ff);
//	f.request_data_spectral();
//	ff.request_data_spectral();

	// Setting explicit right hand side and operator of the left hand side
	sweet::PlaneData_Spectral rhs_u = u;
	sweet::PlaneData_Spectral rhs_v = v;

	rhs_u += - 0.5*dt*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + 0.5*dt*f;
	rhs_v += - 0.5*dt*(u*op.diff_c_x(v)+v*op.diff_c_y(v));

	//std::cout << std::endl << std::endl << "rhs_u" << std::endl;
	//rhs_u.print_physicalArrayData();

	if (shackDict.disc.space_use_spectral_basis_diffs) //spectral
	{

		sweet::PlaneData_Spectral lhs = u;
		lhs = ((-dt)*shackDict.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
        	sweet::PlaneData_Spectral u1 = rhs_u.spectral_div_element_wise(lhs);
        	sweet::PlaneData_Spectral v1 = rhs_v.spectral_div_element_wise(lhs);

        	//std::cout << std::endl << std::endl << "u1" << std::endl;
        	//u1.print_physicalArrayData();

        	io_u = u + dt*shackDict.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
        	      - dt*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1)) +ff*dt;
        	io_v = v + dt*shackDict.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
        	      - dt*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));

        	//std::cout << std::endl << std::endl << "io_u" << std::endl;
        	//io_u.print_physicalArrayData();

	} else { //Jacobi
		SWEETError("NOT available");
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
		sweet::ShackDictionary &i_shackDict,
		PlaneOperators &i_op
)	:
		shackDict(i_shackDict),
		op(i_op)
{
	setup(shackDict.disc.timestepping_order);
}



Burgers_Plane_TS_ln_imex_forcing::~Burgers_Plane_TS_ln_imex_forcing()
{
}

