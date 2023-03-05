/*
 * Burgers_Plane_TS_l_irk_n_sl_forcing.cpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "../burgers_timeintegrators/Burgers_Plane_TS_l_irk_n_sl_forcing.hpp"


void Burgers_Plane_TS_l_irk_n_sl_forcing::run_timestep(
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_u_prev,	///< prognostic variables
		///sweet::PlaneData_Spectral &io_v_prev,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Burgers_Plane_TS_l_irk_n_sl_forcing: Only constant time step size allowed");

	if (i_simulation_timestamp == 0)
	{
#if !SWEET_PARAREAL
		/*
		 * First time step
		 */
		u_prev = io_u;
		v_prev = io_v;
#endif
	}


	//Departure points and arrival points
	ScalarDataArray posx_d(io_u.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_u.planeDataConfig->physical_array_data_number_of_elements);

	Staggering staggering;

	//Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev.toPhys(), v_prev.toPhys(),
			io_u.toPhys(), io_v.toPhys(),
			posx_a, posy_a,
			i_fixed_dt,
			posx_d, posy_d,
			shackDict.sim.plane_domain_size,
			&staggering,
			2,

			shackDict.disc.semi_lagrangian_max_iterations,
			shackDict.disc.semi_lagrangian_convergence_threshold
			);

	// Save old velocities
	u_prev = io_u;
	v_prev = io_v;

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

	sweet::PlaneData_Spectral u=io_u;
	sweet::PlaneData_Spectral v=io_v;

	// Initialize and set timestep dependent source for manufactured solution
	sweet::PlaneData_Spectral f(io_u.planeDataConfig);

	BurgersValidationBenchmarks::set_source(
			shackDict.timecontrol.current_simulation_time+i_fixed_dt,
			shackDict,
			shackDict.disc.space_grid_use_c_staggering,
			f
		);

	// Setting explicit right hand side and operator of the left hand side
	sweet::PlaneData_Spectral rhs_u = u;
	sweet::PlaneData_Spectral rhs_v = v;

	rhs_u += i_fixed_dt*f;

	if (shackDict.disc.space_use_spectral_basis_diffs) //spectral
	{
		sweet::PlaneData_Spectral lhs = u;

		lhs = ((-i_fixed_dt)*shackDict.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		io_u = rhs_u.spectral_div_element_wise(lhs);
		io_v = rhs_v.spectral_div_element_wise(lhs);

	} else { //Jacobi
		SWEETError("NOT available");
	}

}



/*
 * Setup
 */
void Burgers_Plane_TS_l_irk_n_sl_forcing::setup()
{

	// Setup sampler for future interpolations
	sampler2D.setup(shackDict.sim.plane_domain_size, op.planeDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(shackDict.sim.plane_domain_size, op.planeDataConfig);


	sweet::PlaneData_Physical tmp_x(op.planeDataConfig);
	tmp_x.physical_update_lambda_array_indices(
		[&](int i, int j, double &io_data)
		{
			io_data = ((double)i)*shackDict.sim.plane_domain_size[0]/(double)shackDict.disc.space_res_physical[0];
		},
		false
	);

	sweet::PlaneData_Physical tmp_y(op.planeDataConfig);
	tmp_y.physical_update_lambda_array_indices(
		[&](int i, int j, double &io_data)
		{
			io_data = ((double)j)*shackDict.sim.plane_domain_size[1]/(double)shackDict.disc.space_res_physical[1];
		},
		false
	);

	// Initialize arrival points with h position
	ScalarDataArray pos_x = Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_x);
	ScalarDataArray pos_y = Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_y);

	double cell_size_x = shackDict.sim.plane_domain_size[0]/(double)shackDict.disc.space_res_physical[0];
	double cell_size_y = shackDict.sim.plane_domain_size[1]/(double)shackDict.disc.space_res_physical[1];

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*cell_size_x;
	posy_a = pos_y+0.5*cell_size_y;


}


Burgers_Plane_TS_l_irk_n_sl_forcing::Burgers_Plane_TS_l_irk_n_sl_forcing(
		sweet::ShackDictionary &i_shackDict,
		PlaneOperators &i_op
)	:
		shackDict(i_shackDict),
		op(i_op),

		posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements),

		u_prev(i_op.planeDataConfig),
		v_prev(i_op.planeDataConfig)
{
}



Burgers_Plane_TS_l_irk_n_sl_forcing::~Burgers_Plane_TS_l_irk_n_sl_forcing()
{
}

