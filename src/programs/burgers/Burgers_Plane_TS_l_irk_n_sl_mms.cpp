/*
 * Burgers_Plane_TS_l_irk_n_sl_mms.cpp
 *
 *  Created on: 16 April 2018
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_l_irk_n_sl_mms.hpp"


void Burgers_Plane_TS_l_irk_n_sl_mms::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Burgers_Plane_TS_l_irk_n_sl_mms: Only constant time step size allowed");

	if (timestepping_order<1 || timestepping_order>2)
		FatalError("Burgers_Plane_TS_l_irk_n_sl_mms: Only orders 1 and 2 possible");

	if (op.diff_c_y(io_u).reduce_maxAbs()>1e-11)
		FatalError("Burgers_Plane_TS_l_irk_mms: Implemented in 1D only!");

	//Departure points and arrival points
	ScalarDataArray posx_d = posx_a;
	ScalarDataArray posy_d = posy_a;

	double dt = i_fixed_dt;

	Staggering staggering;
	assert(staggering.staggering_type == 'a');

	//Calculate departure points
	if (timestepping_order == 1)
	{
		semiLagrangian.semi_lag_departure_points_first_order(
			io_u, io_v,
			posx_a, posy_a,
			dt,
			posx_d, posy_d,
			staggering
			);
	}
	else
	{
		semiLagrangian.semi_lag_departure_points_settls(
				io_u_prev, io_v_prev,
				io_u, io_v,
				posx_a, posy_a,
				dt,
				posx_d, posy_d,
				staggering
				);
	}

	// Save old velocities
	io_u_prev = io_u;
	io_v_prev = io_v;

	if (timestepping_order == 2)
	{
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_l_irk_mms.run_timestep(
				io_u, io_v,
				io_u_prev, io_v_prev,
				0.5*dt,
				i_simulation_timestamp
		);
	}

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

	if (timestepping_order == 2)
	{
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_l_irk_mms.run_timestep(
				io_u, io_v,
				io_u_prev, io_v_prev,
				0.5*dt,
				i_simulation_timestamp
		);
	}
	else
	{
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_l_irk_mms.run_timestep(
				io_u, io_v,
				io_u_prev, io_v_prev,
				dt,
				i_simulation_timestamp
		);
	}

}


void Burgers_Plane_TS_l_irk_n_sl_mms::return_initial(PlaneData& init)
{
	ts_l_irk_mms.return_initial(init);
}

void Burgers_Plane_TS_l_irk_n_sl_mms::setup_look_up_table(double start, double end, double step_size)
{
	ts_l_irk_mms.setup_look_up_table(start, end, step_size);
}


/*
 * Setup
 */
void Burgers_Plane_TS_l_irk_n_sl_mms::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	ts_l_irk_mms.setup(timestepping_order);

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


Burgers_Plane_TS_l_irk_n_sl_mms::Burgers_Plane_TS_l_irk_n_sl_mms(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements),

		ts_l_irk_mms(simVars, op)
{
	setup(simVars.disc.timestepping_order);
}



Burgers_Plane_TS_l_irk_n_sl_mms::~Burgers_Plane_TS_l_irk_n_sl_mms()
{
}

