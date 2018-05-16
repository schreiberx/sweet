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
	o_u_t.physical_set_zero(); //-i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t.physical_set_zero(); // = 0.0; //-i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

	if (simVars.pde.use_linear_div == 1) // linear div only
		o_h_t.physical_set_zero(); // = 0.0; //-op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	else //nonlinear div
		o_h_t = -i_h*(op.diff_c_x(i_u) + op.diff_c_y(i_v));

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

	// Tmp vars
	PlaneData h(io_h.planeDataConfig);
	PlaneData u(io_h.planeDataConfig);
	PlaneData v(io_h.planeDataConfig);

	//Departure points and arrival points
	ScalarDataArray posx_d(io_h.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_h.planeDataConfig->physical_array_data_number_of_elements);


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

	std::cout << "input: time = " << i_simulation_timestamp  << std::endl;
	std::cout <<  io_h.reduce_sum()  << std::endl;
	std::cout <<  io_u.reduce_sum()  << std::endl;
	std::cout <<  io_v.reduce_sum()  << std::endl;
	std::cout <<  io_h.file_physical_saveData_ascii("h_in.csv")  << std::endl;
	std::cout <<  io_u.file_physical_saveData_ascii("u_in.csv")  << std::endl;
	std::cout <<  io_v.file_physical_saveData_ascii("v_in.csv")  << std::endl;
	std::cout << "-------------------------------"   << std::endl;


	//Preserve io unmodified
	u = io_u;
	v = io_v;
	h = io_h;

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev,	v_prev,
			u,		v,
			posx_a,		posy_a,
			i_dt,
			posx_d,	posy_d,			// output
			simVars.sim.domain_size,
			&staggering,
			simVars.disc.timestepping_order
	);


	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \phi_{0}( \Delta t L ) [
		 * 			U_{0}_dep + \Delta t  (\phi_{1}(-\Delta tL) N(U_{0}))_dep.
		 *
		 *\phi_{1}(-\Delta tL)=psi_{1}(\Delta tL)
		 *
		 *F(U)=N(U)
		 *
		 */

		// Calculate term to be interpolated: u+dt*psi_1(dt L)N(U_{0})
		//Calculate N(U_{0})

		if (simVars.pde.use_linear_div == 0) //Full nonlinear case
		{
			PlaneData FUn_h(planeDataConfig);
			PlaneData FUn_u(planeDataConfig);
			PlaneData FUn_v(planeDataConfig);
			euler_timestep_update_nonlinear(
					h, u, v,
					FUn_h, FUn_u, FUn_v,
					i_simulation_timestamp
			);

			//Apply psi_1 to N(U_{0})
			PlaneData psi1_FUn_h(planeDataConfig);
			PlaneData psi1_FUn_u(planeDataConfig);
			PlaneData psi1_FUn_v(planeDataConfig);

			ts_psi1_rexi.run_timestep(
					FUn_h, FUn_u, FUn_v,
					psi1_FUn_h, psi1_FUn_u, psi1_FUn_v,
					i_dt,
					i_simulation_timestamp
			);

			//Add this to U and interpolate to departure points
			h = h + i_dt*psi1_FUn_h;
			u = u + i_dt*psi1_FUn_u;
			v = v + i_dt*psi1_FUn_v;
		}

		/*
		std::cout << "after calculation of departure points: time = " << i_simulation_timestamp  << std::endl;
		std::cout <<  h.reduce_sum()  << std::endl;
		std::cout <<  u.reduce_sum()  << std::endl;
		std::cout <<  v.reduce_sum()  << std::endl;
		std::cout <<  h.file_physical_saveData_ascii("h_after_dep.csv")  << std::endl;
		std::cout <<  u.file_physical_saveData_ascii("u_after_dep.csv")  << std::endl;
		std::cout <<  v.file_physical_saveData_ascii("v_after_dep.csv")  << std::endl;
		 */

		h = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
		u = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
		v = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);


		std::cout << "after interpolation to departure points: time = " << i_simulation_timestamp  << std::endl;
		std::cout <<  h.file_physical_saveData_ascii("h_after_int.csv")  << std::endl;
		std::cout <<  u.file_physical_saveData_ascii("u_after_int.csv")  << std::endl;
		std::cout <<  v.file_physical_saveData_ascii("v_after_int.csv")  << std::endl;
		std::cout <<  h.reduce_sum()  << std::endl;
		std::cout <<  u.reduce_sum()  << std::endl;
		std::cout <<  v.reduce_sum()  << std::endl;


		//Calculate phi_0 of interpolated U
		PlaneData phi0_Un_h(planeDataConfig);
		PlaneData phi0_Un_u(planeDataConfig);
		PlaneData phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
		);

		h = phi0_Un_h;
		u = phi0_Un_u;
		v = phi0_Un_v;

		//io_h = h;
		//io_u = u;
		//io_v = v;

		std::cout << "after interpolation phi0: time = " << i_simulation_timestamp  << std::endl;
		std::cout <<  h.file_physical_saveData_ascii("h_after_phi0.csv")  << std::endl;
		std::cout <<  u.file_physical_saveData_ascii("u_after_phi0.csv")  << std::endl;
		std::cout <<  v.file_physical_saveData_ascii("v_after_phi0.csv")  << std::endl;
		std::cout <<  h.reduce_sum()  << std::endl;
		std::cout <<  u.reduce_sum()  << std::endl;
		std::cout <<  v.reduce_sum()  << std::endl;


	}
	else if (timestepping_order == 2)
	{

		/*
		 * U_2^{n+1} =  {U}_1^{n+1}
		 *     + \Delta t\, \varphi_0(\Delta t L)\left[ \psi_2(\Delta t L) N({U}_1^{n+1})
		 *     - \left(\psi_2(\Delta t L)N(U^n) \right)_*^n\right],
		 *
		 *     F(U)=N(U)
		 */

		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1}
		//-----------------------------------------------------

		/*
		 * U_{1} = \phi_{0}( \Delta t L ) [
		 * 			U_{0}_dep + \Delta t  (\phi_{1}(-\Delta tL) N(U_{0}))_dep.
		 *
		 *\phi_{1}(-\Delta tL)=psi_{1}(\Delta tL)
		 *
		 *F(U)=N(U)
		 *
		 */

		// Calculate term to be interpolated: u+dt*psi_1(dt L)N(U_{0})

		//Calculate N(U_{0})
		PlaneData FUn_h(planeDataConfig);
		PlaneData FUn_u(planeDataConfig);
		PlaneData FUn_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				h, u, v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		//Apply psi_1 to N(U_{0})
		PlaneData psi1_FUn_h(planeDataConfig);
		PlaneData psi1_FUn_u(planeDataConfig);
		PlaneData psi1_FUn_v(planeDataConfig);

		ts_psi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				psi1_FUn_h, psi1_FUn_u, psi1_FUn_v,
				i_dt,
				i_simulation_timestamp
		);

		//Add this to U and interpolate to departure points
		PlaneData h_adv(planeDataConfig);
		PlaneData u_adv(planeDataConfig);
		PlaneData v_adv(planeDataConfig);
		h_adv = h + i_dt*psi1_FUn_h;
		u_adv = u + i_dt*psi1_FUn_u;
		v_adv = v + i_dt*psi1_FUn_v;

		PlaneData h_dep(io_h.planeDataConfig);
		PlaneData u_dep(io_h.planeDataConfig);
		PlaneData v_dep(io_h.planeDataConfig);

		h_dep = sampler2D.bicubic_scalar(h_adv, posx_d, posy_d, -0.5, -0.5);
		u_dep = sampler2D.bicubic_scalar(u_adv, posx_d, posy_d, -0.5, -0.5);
		v_dep = sampler2D.bicubic_scalar(v_adv, posx_d, posy_d, -0.5, -0.5);

		//Calculate phi_0 of interpolated U
		PlaneData phi0_Un_h(planeDataConfig);
		PlaneData phi0_Un_u(planeDataConfig);
		PlaneData phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				h_dep, u_dep, v_dep,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
		);

		//Save sl-etdrk1 solution
		PlaneData h_1(io_h.planeDataConfig);
		PlaneData u_1(io_h.planeDataConfig);
		PlaneData v_1(io_h.planeDataConfig);

		h_1 = phi0_Un_h;
		u_1 = phi0_Un_u;
		v_1 = phi0_Un_v;

		//Calculate psi2NU_1
		//-----------------------

		//NU_1
		PlaneData FU1_h(planeDataConfig);
		PlaneData FU1_u(planeDataConfig);
		PlaneData FU1_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				h_1, u_1, v_1,
				FU1_h, FU1_u, FU1_v,
				i_simulation_timestamp
		);
		//Apply psi2
		PlaneData psi2_FU1_h(planeDataConfig);
		PlaneData psi2_FU1_u(planeDataConfig);
		PlaneData psi2_FU1_v(planeDataConfig);
		ts_psi2_rexi.run_timestep(
				FU1_h, FU1_u, FU1_v,
				psi2_FU1_h, psi2_FU1_u, psi2_FU1_v,
				i_dt,
				i_simulation_timestamp
		);

		//Calculate psi2NUn

		PlaneData psi2_FUn_h(planeDataConfig);
		PlaneData psi2_FUn_u(planeDataConfig);
		PlaneData psi2_FUn_v(planeDataConfig);

		ts_psi2_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				psi2_FUn_h, psi2_FUn_u, psi2_FUn_v,
				i_dt,
				i_simulation_timestamp
		);

		//Interpolate psi2NUn to departure points ( )_*

		PlaneData psi2FUn_h_dep(io_h.planeDataConfig);
		PlaneData psi2FUn_u_dep(io_u.planeDataConfig);
		PlaneData psi2FUn_v_dep(io_v.planeDataConfig);

		psi2FUn_h_dep = sampler2D.bicubic_scalar(psi2_FUn_h, posx_d, posy_d, -0.5, -0.5);
		psi2FUn_u_dep = sampler2D.bicubic_scalar(psi2_FUn_u, posx_d, posy_d, -0.5, -0.5);
		psi2FUn_v_dep = sampler2D.bicubic_scalar(psi2_FUn_v, posx_d, posy_d, -0.5, -0.5);


		//psi2NU_1-psi2NUn_dep
		PlaneData dif2_h(io_h.planeDataConfig);
		PlaneData dif2_u(io_u.planeDataConfig);
		PlaneData dif2_v(io_v.planeDataConfig);

		dif2_h = psi2_FU1_h-psi2FUn_h_dep;
		dif2_u = psi2_FU1_u-psi2FUn_u_dep;
		dif2_v = psi2_FU1_v-psi2FUn_v_dep;

		//Apply phi0
		PlaneData phi0_dif2_h(planeDataConfig);
		PlaneData phi0_dif2_u(planeDataConfig);
		PlaneData phi0_dif2_v(planeDataConfig);

		ts_phi0_rexi.run_timestep(
				dif2_h, dif2_u, dif2_v,
				phi0_dif2_h, phi0_dif2_u, phi0_dif2_v,
				i_dt,
				i_simulation_timestamp
			);
		std::cout <<  h_1.reduce_maxAbs()  << std::endl;
		std::cout <<  u_1.reduce_maxAbs()  << std::endl;
		std::cout <<  v_1.reduce_maxAbs()  << std::endl;

		io_h = h_1 + i_dt*phi0_dif2_h;
		io_u = u_1 + i_dt*phi0_dif2_u;
		io_v = v_1 + i_dt*phi0_dif2_v;
	}
	else if (timestepping_order == 4)
	{
		double dt = i_dt;
		double dt_half = dt*0.5;

		FatalError("TODO: Order 4 is not implemented, yet!");


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

	// Save current time step for next step
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	io_h = h;
	io_u = u;
	io_v = v;

}



/*
 * Setup
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::setup(
		//REXI_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order
)
{
	timestepping_order = i_timestepping_order;


	if (timestepping_order == 1)
	{
		ts_phi0_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size);
		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size);
		ts_phi2_rexi.setup(simVars.rexi, "phi2", simVars.timecontrol.current_timestep_size);

		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size);
		ts_psi2_rexi.setup(simVars.rexi, "psi2", simVars.timecontrol.current_timestep_size);

	}
	else if (timestepping_order == 4)
	{
		ts_phi0_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size*0.5);
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size*0.5);
		ts_phi2_rexi.setup(simVars.rexi, "phi2", simVars.timecontrol.current_timestep_size*0.5);

		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size*0.5);
		ts_psi2_rexi.setup(simVars.rexi, "psi2", simVars.timecontrol.current_timestep_size*0.5);
		ts_psi3_rexi.setup(simVars.rexi, "psi3", simVars.timecontrol.current_timestep_size*0.5);

		ts_ups0_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);
		ts_ups1_rexi.setup(simVars.rexi, "ups1", simVars.timecontrol.current_timestep_size);
		ts_ups2_rexi.setup(simVars.rexi, "ups2", simVars.timecontrol.current_timestep_size);
		ts_ups3_rexi.setup(simVars.rexi, "ups3", simVars.timecontrol.current_timestep_size);
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

		ts_psi1_rexi(simVars, op),
		ts_psi2_rexi(simVars, op),
		ts_psi3_rexi(simVars, op),

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

