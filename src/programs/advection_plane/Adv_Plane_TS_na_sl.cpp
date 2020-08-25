/*
 * Adv_Plane_TS_na_sl.cpp
 *
 *  Created on: 29 Mar 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "Adv_Plane_TS_na_sl.hpp"




void Adv_Plane_TS_na_sl::run_timestep(
		PlaneData &io_phi,		///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Only constant time step size allowed");

	double dt = simVars.timecontrol.current_timestep_size;

	if (simVars.benchmark.getExternalForcesCallback != nullptr)
	{
		simVars.benchmark.getExternalForcesCallback(1, i_simulation_timestamp, &io_u, simVars.benchmark.getExternalForcesUserData);
		simVars.benchmark.getExternalForcesCallback(2, i_simulation_timestamp, &io_v, simVars.benchmark.getExternalForcesUserData);
	}

	if (i_simulation_timestamp == 0)
	{
		prog_u_prev = io_u;
		prog_v_prev = io_v;
	}


	// OUTPUT: position of departure points at t
	ScalarDataArray posx_d(io_phi.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_phi.planeDataConfig->physical_array_data_number_of_elements);

	semiLagrangian.semi_lag_departure_points_settls(
			prog_u_prev, prog_v_prev,
			io_u, io_v,
			posx_a, posy_a,
			dt,
			posx_d, posy_d,
			simVars.sim.plane_domain_size,
			nullptr,
			timestepping_order,

			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold
	);

	prog_u_prev = io_u;
	prog_v_prev = io_v;

	PlaneData new_prog_phi(io_phi.planeDataConfig);

	if (timestepping_order == 1)
	{
		sampler2D.bilinear_scalar(
				io_phi,
				posx_d,
				posy_d,
				new_prog_phi
		);
	}
	else if (timestepping_order == 2)
	{
		sampler2D.bicubic_scalar(
				io_phi,
				posx_d,
				posy_d,
				new_prog_phi
		);
	}
	else
	{
		SWEETError("Timestepping order not available");
	}

	io_phi = new_prog_phi;
}



/*
 * Setup
 */
void Adv_Plane_TS_na_sl::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	const PlaneDataConfig *planeDataConfig = op.planeDataConfig;

	posx_a.setup(planeDataConfig->physical_array_data_number_of_elements);
	posy_a.setup(planeDataConfig->physical_array_data_number_of_elements);

	// setup some test sampling points
	// we use 2 arrays - one for each sampling position
	posx_a.update_lambda_array_indices(
		[&](int idx, double &io_data)
		{
			int i = idx % planeDataConfig->physical_res[0];

			io_data = (double)i*(double)simVars.sim.plane_domain_size[0]/(double)planeDataConfig->physical_res[0];

			assert(io_data >= 0.0);
			assert(io_data <= simVars.sim.plane_domain_size[0]);
		}
	);
	posy_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
		{
			//int i = idx % planeDataConfig->physical_data_size[0];
			int j = idx / (double)planeDataConfig->physical_res[0];

			io_data = (double)j*(double)simVars.sim.plane_domain_size[1]/(double)planeDataConfig->physical_res[1];

			assert(io_data >= 0.0);
			assert(io_data <= simVars.sim.plane_domain_size[1]);
		}
	);

	// TODO: Use semiLagrangian.sampler2D
	sampler2D.setup(simVars.sim.plane_domain_size, planeDataConfig);

	//PXT- This just calls sampler2D.setup, so any reason for having it?
	semiLagrangian.setup(simVars.sim.plane_domain_size, planeDataConfig);
}


Adv_Plane_TS_na_sl::Adv_Plane_TS_na_sl(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		prog_u_prev(op.planeDataConfig),
		prog_v_prev(op.planeDataConfig)
{
	setup(simVars.disc.timestepping_order);
}



Adv_Plane_TS_na_sl::~Adv_Plane_TS_na_sl()
{
}

