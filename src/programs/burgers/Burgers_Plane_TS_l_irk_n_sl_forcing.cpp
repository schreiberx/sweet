/*
 * Burgers_Plane_TS_l_irk_n_sl_forcing.cpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_l_irk_n_sl_forcing.hpp"


void Burgers_Plane_TS_l_irk_n_sl_forcing::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Burgers_Plane_TS_l_irk_n_sl_forcing: Only constant time step size allowed");

	//Departure points and arrival points
	ScalarDataArray posx_d(io_u.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_u.planeDataConfig->physical_array_data_number_of_elements);

	Staggering staggering;

	//Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			io_u_prev, io_v_prev,
			io_u, io_v,
			posx_a, posy_a,
			i_fixed_dt,
			posx_d, posy_d,
			simVars.sim.domain_size,
			&staggering,
			2
			);

	// Save old velocities
	io_u_prev = io_u;
	io_v_prev = io_v;

	//Now interpolate to the the departure points
	//Departure points are set for physical space
	io_u = sampler2D.bicubic_scalar(
			io_u,
			posx_d,
			posy_d,
			staggering.u[0],
			staggering.u[1]
	);

	io_v = sampler2D.bicubic_scalar(
			io_v,
			posx_d,
			posy_d,
			staggering.v[0],
			staggering.v[1]
	);

	PlaneData u=io_u;
	PlaneData v=io_v;

	// Initialize and set timestep dependent source for manufactured solution
	PlaneData f(io_u.planeDataConfig);

	BurgersValidationBenchmarks::set_source(
			simVars.timecontrol.current_simulation_time+i_fixed_dt,
			simVars,
			simVars.disc.use_staggering,
			f
		);

	f.request_data_spectral();

	// Setting explicit right hand side and operator of the left hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;

	rhs_u += i_fixed_dt*f;

	if (simVars.disc.use_spectral_basis_diffs) //spectral
	{
		PlaneData lhs = u;

		lhs = ((-i_fixed_dt)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		io_u = rhs_u.spectral_div_element_wise(lhs);
		io_v = rhs_v.spectral_div_element_wise(lhs);

	} else { //Jacobi
		FatalError("NOT available");
	}

}



/*
 * Setup
 */
void Burgers_Plane_TS_l_irk_n_sl_forcing::setup()
{

	// Setup sampler for future interpolations
	sampler2D.setup(simVars.sim.domain_size, op.planeDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(simVars.sim.domain_size, op.planeDataConfig);


	PlaneData tmp_x(op.planeDataConfig);
	tmp_x.physical_update_lambda_array_indices(
		[&](int i, int j, double &io_data)
		{
			io_data = ((double)i)*simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
		},
		false
	);

	PlaneData tmp_y(op.planeDataConfig);
	tmp_y.physical_update_lambda_array_indices(
		[&](int i, int j, double &io_data)
		{
			io_data = ((double)j)*simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];
		},
		false
	);

	// Initialize arrival points with h position
	ScalarDataArray pos_x = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_x);
	ScalarDataArray pos_y = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_y);

	double cell_size_x = simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
	double cell_size_y = simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*cell_size_x;
	posy_a = pos_y+0.5*cell_size_y;


}


Burgers_Plane_TS_l_irk_n_sl_forcing::Burgers_Plane_TS_l_irk_n_sl_forcing(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements)
{
}



Burgers_Plane_TS_l_irk_n_sl_forcing::~Burgers_Plane_TS_l_irk_n_sl_forcing()
{
}

