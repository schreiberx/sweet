/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_etdrk.hpp
 *
 *  Created on: 09 Oct 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *  Changelog:
 *      based on Martin Schreiber ETD timestepper
 */

#include "SWE_Plane_TS_l_rexi_na_sl_nd_etdrk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::euler_timestep_update_nonlinear(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	/*
	 * o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	 * o_u_t = -i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	 * o_v_t = -i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);
	 */
	// In lagrangian form, the only nonlinearity is the nonlinear divergence
	o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	o_u_t = 0.0; //-i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t = 0.0; //-i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

}



void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Plane_TS_l_rexi_na_sl_nd_etdrk: Only constant time step size allowed (Please set --dt)");


	const PlaneDataConfig *planeDataConfig = io_h.planeDataConfig;

	//Departure points and arrival points
	ScalarDataArray posx_d(io_h.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_h.planeDataConfig->physical_array_data_number_of_elements);

	// Tmp vars
	PlaneData h(io_h.planeDataConfig);
	PlaneData u(io_h.planeDataConfig);
	PlaneData v(io_h.planeDataConfig);

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

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev,	v_prev,
			io_u,		io_v,
			posx_a,		posy_a,
			i_dt,
			posx_d,	posy_d,			// output
			&staggering
	);

	u = io_u;
	v = io_v;
	h = io_h;

	// Interpolate U to departure points
	u = sampler2D.bicubic_scalar(io_u, posx_d, posy_d, -0.5, -0.5);
	v = sampler2D.bicubic_scalar(io_v, posx_d, posy_d, -0.5, -0.5);
	h = sampler2D.bicubic_scalar(io_h, posx_d, posy_d, -0.5, -0.5);

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		PlaneData phi0_Un_h(planeDataConfig);
		PlaneData phi0_Un_u(planeDataConfig);
		PlaneData phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
			);

		PlaneData FUn_h(planeDataConfig);
		PlaneData FUn_u(planeDataConfig);
		PlaneData FUn_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				h, u, v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		PlaneData phi1_FUn_h(planeDataConfig);
		PlaneData phi1_FUn_u(planeDataConfig);
		PlaneData phi1_FUn_v(planeDataConfig);

		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_dt,
				i_simulation_timestamp
			);

		io_h = phi0_Un_h + i_dt*phi1_FUn_h;
		io_u = phi0_Un_u + i_dt*phi1_FUn_u;
		io_v = phi0_Un_v + i_dt*phi1_FUn_v;
	}
	else if (timestepping_order == 2)
	{

		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		PlaneData phi0_Un_h(planeDataConfig);
		PlaneData phi0_Un_u(planeDataConfig);
		PlaneData phi0_Un_v(planeDataConfig);

		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
			);

		PlaneData FUn_h(planeDataConfig);
		PlaneData FUn_u(planeDataConfig);
		PlaneData FUn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				h, u, v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		PlaneData phi1_FUn_h(planeDataConfig);
		PlaneData phi1_FUn_u(planeDataConfig);
		PlaneData phi1_FUn_v(planeDataConfig);

		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_dt,
				i_simulation_timestamp
			);

		PlaneData A_h = phi0_Un_h + i_dt*phi1_FUn_h;
		PlaneData A_u = phi0_Un_u + i_dt*phi1_FUn_u;
		PlaneData A_v = phi0_Un_v + i_dt*phi1_FUn_v;

		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		PlaneData FAn_h(planeDataConfig);
		PlaneData FAn_u(planeDataConfig);
		PlaneData FAn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


		PlaneData phi2_X_h(planeDataConfig);
		PlaneData phi2_X_u(planeDataConfig);
		PlaneData phi2_X_v(planeDataConfig);

		ts_phi2_rexi.run_timestep(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h,
				phi2_X_u,
				phi2_X_v,

				i_dt,
				i_simulation_timestamp
			);

		io_h = A_h + i_dt*phi2_X_h;
		io_u = A_u + i_dt*phi2_X_u;
		io_v = A_v + i_dt*phi2_X_v;
	}
	else if (timestepping_order == 4)
	{
		double dt = i_dt;
		double dt_half = dt*0.5;



		/*
		 * Precompute commonly used terms
		 */
		PlaneData phi0_Un_h(planeDataConfig);
		PlaneData phi0_Un_u(planeDataConfig);
		PlaneData phi0_Un_v(planeDataConfig);

		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				dt_half,
				i_simulation_timestamp
			);

		PlaneData FUn_h(planeDataConfig);
		PlaneData FUn_u(planeDataConfig);
		PlaneData FUn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				h, u, v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);



		/*
		 * Some commonly shared buffers
		 */

		PlaneData phi1_h(planeDataConfig);
		PlaneData phi1_u(planeDataConfig);
		PlaneData phi1_v(planeDataConfig);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 */
		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		PlaneData A_h = phi0_Un_h + dt_half*phi1_h;
		PlaneData A_u = phi0_Un_u + dt_half*phi1_u;
		PlaneData A_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 */

		PlaneData FAn_h(planeDataConfig);
		PlaneData FAn_u(planeDataConfig);
		PlaneData FAn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.run_timestep(
				FAn_h, FAn_u, FAn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		PlaneData B_h = phi0_Un_h + dt_half*phi1_h;
		PlaneData B_u = phi0_Un_u + dt_half*phi1_u;
		PlaneData B_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 */

		PlaneData phi0_An_h(planeDataConfig);
		PlaneData phi0_An_u(planeDataConfig);
		PlaneData phi0_An_v(planeDataConfig);

		ts_phi0_rexi.run_timestep(
				A_h, A_u, A_v,
				phi0_An_h, phi0_An_u, phi0_An_v,
				dt_half,
				i_simulation_timestamp
			);


		PlaneData FBn_h(planeDataConfig);
		PlaneData FBn_u(planeDataConfig);
		PlaneData FBn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				B_h, B_u, B_v,
				FBn_h, FBn_u, FBn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.run_timestep(
				2.0*FBn_h - FUn_h,
				2.0*FBn_u - FUn_u,
				2.0*FBn_v - FUn_v,
				phi1_h,	phi1_u,	phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		PlaneData C_h = phi0_An_h + dt_half*phi1_h;
		PlaneData C_u = phi0_An_u + dt_half*phi1_u;
		PlaneData C_v = phi0_An_v + dt_half*phi1_v;



		/*
		 * R0 - R3
		 */
		PlaneData FCn_h(planeDataConfig);
		PlaneData FCn_u(planeDataConfig);
		PlaneData FCn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				C_h, C_u, C_v,
				FCn_h, FCn_u, FCn_v,
				i_simulation_timestamp + dt
		);

		PlaneData R0_h = h;
		PlaneData R0_u = u;
		PlaneData R0_v = v;

		PlaneData &R1_h = FUn_h;
		PlaneData &R1_u = FUn_u;
		PlaneData &R1_v = FUn_v;

		PlaneData R2_h = FAn_h + FBn_h;
		PlaneData R2_u = FAn_u + FBn_u;
		PlaneData R2_v = FAn_v + FBn_v;

		PlaneData &R3_h = FCn_h;
		PlaneData &R3_u = FCn_u;
		PlaneData &R3_v = FCn_v;


		/*
		 * U_{n+1} =
		 * 		\psi_{0}(\Delta tL)R_{0}
		 * 			+ \Delta t
		 * 			(
		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
		 * 				  \upsilon_{3}(\Delta tL) R_{3}
		 * 			)
		 */
		ts_ups0_rexi.run_timestep(
				R0_h, R0_u, R0_v,
				dt,		i_simulation_timestamp
			);

		ts_ups1_rexi.run_timestep(
				R1_h, R1_u, R1_v,
				dt,		i_simulation_timestamp
			);

		ts_ups2_rexi.run_timestep(
				R2_h, R2_u, R2_v,
				dt,		i_simulation_timestamp
			);

		ts_ups3_rexi.run_timestep(
				R3_h, R3_u, R3_v,
				dt,		i_simulation_timestamp
			);

		io_h = R0_h + dt*(R1_h + 2.0*R2_h + R3_h);
		io_u = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
		io_v = R0_v + dt*(R1_v + 2.0*R2_v + R3_v);
	}
	else
	{
		FatalError("TODO: This order is not implemented, yet!");
	}
	// Set time (n) as time (n-1)
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

}



/*
 * Setup
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::setup(
		REXI_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order
)
{
	timestepping_order = i_timestepping_order;


	if (timestepping_order == 1)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", simVars.timecontrol.current_timestep_size);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", simVars.timecontrol.current_timestep_size);
		ts_phi1_rexi.setup(i_rexiSimVars, "phi1", simVars.timecontrol.current_timestep_size);
		ts_phi2_rexi.setup(i_rexiSimVars, "phi2", simVars.timecontrol.current_timestep_size);
	}
	else if (timestepping_order == 4)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", simVars.timecontrol.current_timestep_size*0.5);
		ts_phi1_rexi.setup(i_rexiSimVars, "phi1", simVars.timecontrol.current_timestep_size*0.5);
		ts_phi2_rexi.setup(i_rexiSimVars, "phi2", simVars.timecontrol.current_timestep_size*0.5);

		ts_ups0_rexi.setup(i_rexiSimVars, "phi0", simVars.timecontrol.current_timestep_size);
		ts_ups1_rexi.setup(i_rexiSimVars, "ups1", simVars.timecontrol.current_timestep_size);
		ts_ups2_rexi.setup(i_rexiSimVars, "ups2", simVars.timecontrol.current_timestep_size);
		ts_ups3_rexi.setup(i_rexiSimVars, "ups3", simVars.timecontrol.current_timestep_size);
	}

	if (simVars.disc.use_staggering)
		FatalError("SWE_Plane_TS_l_rexi_na_sl_nd_etdrk: Staggering not supported for l_rexi_na_sl_nd_etdrk");

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


SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::SWE_Plane_TS_l_rexi_na_sl_nd_etdrk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		ts_phi0_rexi(simVars, op),
		ts_phi1_rexi(simVars, op),
		ts_phi2_rexi(simVars, op),

		ts_ups0_rexi(simVars, op),
		ts_ups1_rexi(simVars, op),
		ts_ups2_rexi(simVars, op),
		ts_ups3_rexi(simVars, op),

		h_prev(i_op.planeDataConfig),
		u_prev(i_op.planeDataConfig),
		v_prev(i_op.planeDataConfig),

		posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

		posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements)
{
}



SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::~SWE_Plane_TS_l_rexi_na_sl_nd_etdrk()
{
}

