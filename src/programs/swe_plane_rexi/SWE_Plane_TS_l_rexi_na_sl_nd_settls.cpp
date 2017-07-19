/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_rexi_na_sl_nd_settls.hpp"


void SWE_Plane_TS_l_rexi_na_sl_nd_settls::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt <= 0)
		FatalError("SWE_Plane_TS_l_rexi_na_sl_nd_erk: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time - i_simulation_timestamp;

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step
		 */

		h_prev = io_h;
		u_prev = io_u;
		v_prev = io_v;
	}


	//Out vars
	PlaneData h(io_h.planeDataConfig);
	PlaneData u(io_h.planeDataConfig);
	PlaneData v(io_h.planeDataConfig);
	PlaneData N_h(io_h.planeDataConfig);
	PlaneData N_u(io_h.planeDataConfig);
	PlaneData N_v(io_h.planeDataConfig);
	PlaneData hdiv(io_h.planeDataConfig);

	//Departure points and arrival points
	ScalarDataArray posx_d(io_h.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_h.planeDataConfig->physical_array_data_number_of_elements);

	//Parameters
	double dt = i_fixed_dt;

	Staggering staggering;
	assert(staggering.staggering_type == 'a');

	if (with_nonlinear > 0)
	{
		// Calculate departure points
		semiLagrangian.semi_lag_departure_points_settls(
				u_prev,	v_prev,
				io_u,		io_v,
				posx_a,		posy_a,
				dt,
				posx_d,	posy_d,			// output
				staggering
		);
	}

	u = io_u;
	v = io_v;
	h = io_h;

	N_u.physical_set_all(0);
	N_v.physical_set_all(0);
	N_h.physical_set_all(0);
	hdiv.physical_set_all(0);

	//Calculate nonlinear terms
	if (with_nonlinear > 0)
	{
		//Calculate Divergence and vorticity spectrally
		hdiv =  - io_h * (op.diff_c_x(io_u) + op.diff_c_y(io_v));

		// Calculate nonlinear term for the previous time step
		// h*div is calculate in cartesian space (pseudo-spectrally)
		N_h = -h_prev * (op.diff_c_x(u_prev) + op.diff_c_y(v_prev));

		//Calculate exp(Ldt)N(n-1), relative to previous timestep
		//Calculate the V{n-1} term as in documentation, with the exponential integrator
		double o_dt;
		ts_l_rexi.run_timestep(N_h, N_u, N_v, o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time);

		//Use N_h to store now the nonlinearity of the current time (prev will not be required anymore)
		//Update the nonlinear terms with the constants relative to dt
		N_u = -0.5 * dt * N_u; // N^n of u term is zero
		N_v = -0.5 * dt * N_v; // N^n of v term is zero
		N_h = dt * hdiv - 0.5 * dt * N_h ; //N^n of h has the nonlin term
	}



	if (with_nonlinear > 0)
	{
		//Build variables to be interpolated to dep. points
		// This is the W^n term in documentation
		u = u + N_u;
		v = v + N_v;
		h = h + N_h;

		// Interpolate W to departure points
		u = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
		v = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);

		h = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);

	}

	/*
	 * Calculate the exp(Ldt) W{n}_* term as in documentation, with the exponential integrator?
	 */
	ts_l_rexi.run_timestep(h, u, v, o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time);

	if (with_nonlinear > 1)
	{
		// Add nonlinearity in h
		h = h + 0.5 * dt * hdiv;
	}

	// Set time (n) as time (n-1)
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	// output data
	io_h = h;
	io_u = u;
	io_v = v;

	o_dt = i_fixed_dt;
}



/*
 * Setup
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_settls::setup(
		REXI_SimulationVariables &i_rexi,

		int i_with_nonlinear
)
{
	ts_l_rexi.setup(simVars.rexi);

	with_nonlinear = i_with_nonlinear;

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


SWE_Plane_TS_l_rexi_na_sl_nd_settls::SWE_Plane_TS_l_rexi_na_sl_nd_settls(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		ts_l_rexi(i_simVars, i_op),

		h_prev(i_op.planeDataConfig),
		u_prev(i_op.planeDataConfig),
		v_prev(i_op.planeDataConfig),

		posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements)
{
}



SWE_Plane_TS_l_rexi_na_sl_nd_settls::~SWE_Plane_TS_l_rexi_na_sl_nd_settls()
{
}

