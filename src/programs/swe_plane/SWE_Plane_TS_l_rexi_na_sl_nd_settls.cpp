/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_settls.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */


#include "SWE_Plane_TS_l_rexi_na_sl_nd_settls.hpp"


void SWE_Plane_TS_l_rexi_na_sl_nd_settls::run_timestep(
		PlaneData &io_h,	///< prognostic variables - perturbation of height!
		PlaneData &io_u,	///< prognostic variables - zonal velocity
		PlaneData &io_v,	///< prognostic variables - meridional velocity

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Plane_TS_l_rexi_na_sl_nd_settls: Only constant time step size allowed (Please set --dt)");

	const PlaneDataConfig *planeDataConfig = io_h.planeDataConfig;

	//Out vars
	PlaneData h(io_h.planeDataConfig);
	PlaneData u(io_h.planeDataConfig);
	PlaneData v(io_h.planeDataConfig);

	//Temporary
	PlaneData N_h(io_h.planeDataConfig);
	PlaneData N_u(io_h.planeDataConfig);
	PlaneData N_v(io_h.planeDataConfig);
	PlaneData hdiv(io_h.planeDataConfig);
	PlaneData N_h_prev(io_h.planeDataConfig);
	PlaneData N_h_ext(io_h.planeDataConfig);

	//Departure points and arrival points
	ScalarDataArray posx_d(io_h.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_h.planeDataConfig->physical_array_data_number_of_elements);

	//Parameters
	double dt = i_dt;

	Staggering staggering;
	assert(staggering.staggering_type == 'a');

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step
		 */
		h_prev = io_h;
		u_prev = io_u;
		v_prev = io_v;
	}

	//Preserve io unmodified
	u = io_u;
	v = io_v;
	h = io_h;

	// Calculate departure points
	//Calculate departure points - always force to be second order accurate!
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev,	v_prev,
			u,		v,
			posx_a,		posy_a,
			dt,
			posx_d,	posy_d,			// output
			simVars.sim.domain_size,
			&staggering,
			2 //simVars.disc.timestepping_order
	);

	N_u.physical_set_all(0);
	N_v.physical_set_all(0);
	N_h.physical_set_all(0);
	N_h_prev.physical_set_all(0);
	N_h_ext.physical_set_all(0);
	hdiv.physical_set_all(0);

	//Original more stable scheme

	//Calculate nonlinear terms - not done in case of only linear divergence (linear div is already in linear part)
	if (simVars.pde.use_linear_div == 0) // Full nonlinear case
	{
		//FatalError("Not yet tested");

		// Calculate nonlinear term for the previous time step
		N_h = -h_prev * (op.diff_c_x(u_prev) + op.diff_c_y(v_prev));

		if(simVars.misc.use_local_visc != 0)
			N_h= op.implicit_diffusion(N_h, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);

		//Calculate exp(Ldt)N(n-1), relative to previous timestep
		ts_l_rexi.run_timestep(N_h, N_u, N_v, i_dt, i_simulation_timestamp);

		//Use N_h to store now the nonlinearity of the current time (prev will not be required anymore)
		//Update the nonlinear terms with the constants relative to dt
		// N=dtN^n-0.5dt exp(dtL)N^n-1 from paper
		// N=-h*div is calculate in cartesian space (pseudo-spectrally)
		hdiv =  - h * (op.diff_c_x(io_u) + op.diff_c_y(io_v));
		if(simVars.misc.use_local_visc != 0)
			hdiv= op.implicit_diffusion(hdiv, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);

		N_u = -0.5 * dt * N_u; // N^n of u term is zero
		N_v = -0.5 * dt * N_v; // N^n of v term is zero
		N_h = dt * hdiv - 0.5 * dt * N_h ; //N^n of h has the nonlin term

		//Build variables to be interpolated to dep. points
		// This is the W^n term in documentation
		u = u + N_u;
		v = v + N_v;
		h = h + N_h;
	}

	// Interpolate W to departure points
	h = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
	u = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
	v = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);


	// Add nonlinearity in h
	if (simVars.pde.use_linear_div == 0) // Full nonlinear case
	{
		h = h + 0.5 * dt * hdiv;
	}

	/*
	 * Calculate the exp(Ldt) of the resulting u,v,h
	 */
	//Calculate phi_0 of interpolated U
	PlaneData phi0_Un_h(planeDataConfig);
	PlaneData phi0_Un_u(planeDataConfig);
	PlaneData phi0_Un_v(planeDataConfig);
	//ts_l_rexi.run_timestep(h, u, v, i_dt, i_simulation_timestamp);
	ts_l_rexi.run_timestep(
			h, u, v,
			phi0_Un_h, phi0_Un_u, phi0_Un_v,
			i_dt,
			i_simulation_timestamp
	);

	h = phi0_Un_h;
	u = phi0_Un_u;
	v = phi0_Un_v;

	// Set time (n) as time (n-1)
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	// output data
	io_h = h;
	io_u = u;
	io_v = v;
}



/*
 * Setup
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_settls::setup(
		//REXI_SimulationVariables &i_rexi,
		//int i_use_linear_div
)
{
	ts_l_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);

	if (simVars.disc.use_staggering)
		FatalError("SWE_Plane_TS_l_rexi_na_sl_nd_settls: Staggering not supported for l_rexi_na_sl_nd_settls");

	//with_linear_div_only = i_use_linear_div;

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

