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

	// linear div only
	if (use_only_linear_divergence)
	{
		o_h_t.physical_set_zero(); // = 0.0; //-op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	}
	else
	{
		// nonlinear div
		o_h_t = -i_h*(op.diff_c_x(i_u) + op.diff_c_y(i_v));
		// Smooth spectrum to avoid instability
		if (simVars.misc.use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
			FatalError("Implicit diffusion only supported with spectral space activated");
#else
			o_h_t= op.implicit_diffusion(o_h_t, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
#endif
		}

	}
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
	//h, u, v tmp
	PlaneData h(planeDataConfig);
	PlaneData u(planeDataConfig);
	PlaneData v(planeDataConfig);
	//Nonlinear calculation of u,v,h from input
	PlaneData FUn_h(planeDataConfig);
	PlaneData FUn_u(planeDataConfig);
	PlaneData FUn_v(planeDataConfig);


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


	//Preserve io unmodified
	u = io_u;
	v = io_v;
	h = io_h;

	// Calculate departure points - always force to be second order accurate!
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev,	v_prev,
			u,		v,
			posx_a,		posy_a,
			i_dt,
			posx_d,	posy_d,			// output
			simVars.sim.plane_domain_size,
			&staggering,
			2, //simVars.disc.timestepping_order,

			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold
	);

	if (timestepping_order == 1 || timestepping_order == 2)
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

		if (!use_only_linear_divergence) //Full nonlinear case
		{

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



		h = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
		u = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
		v = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);


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
	}
	else
	{
		FatalError("TODO: This order is not implemented, yet!");
	}



	//Aditional steps for SL-ETD2RK (it depends on SL-ETD1RK) - only for full nonlinear case
	if (timestepping_order == 2 && !use_only_linear_divergence)
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
		// Save SL-ETD1RK from above
		PlaneData h1(planeDataConfig);
		PlaneData u1(planeDataConfig);
		PlaneData v1(planeDataConfig);
		h1 = h;
		u1 = u;
		v1 = v;

		//Calculate psi2NU_1
		//-----------------------

		//NU_1
		PlaneData FU1_h(planeDataConfig);
		PlaneData FU1_u(planeDataConfig);
		PlaneData FU1_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				h1, u1, v1,
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

		PlaneData psi2FUn_h_dep(planeDataConfig);
		PlaneData psi2FUn_u_dep(planeDataConfig);
		PlaneData psi2FUn_v_dep(planeDataConfig);
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
		//second order final forcing
		PlaneData phi0_dif2_h(planeDataConfig);
		PlaneData phi0_dif2_u(planeDataConfig);
		PlaneData phi0_dif2_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				dif2_h, dif2_u, dif2_v,
				phi0_dif2_h, phi0_dif2_u, phi0_dif2_v,
				i_dt,
				i_simulation_timestamp
		);

		h = h1 + i_dt*phi0_dif2_h;
		u = u1 + i_dt*phi0_dif2_u;
		v = v1 + i_dt*phi0_dif2_v;

	}



	// Save current time step for next step
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	io_h = h;
	io_u = u;
	io_v = v;


}



/**
 * Setup
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::setup(
	int i_timestepping_order,
	bool i_use_only_linear_divergence
)
{
	timestepping_order = i_timestepping_order;
	use_only_linear_divergence = i_use_only_linear_divergence;

	ts_phi0_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);

	if (timestepping_order == 1 && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size);
		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size);
	}
	else if (timestepping_order == 2 && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size);
		ts_phi2_rexi.setup(simVars.rexi, "phi2", simVars.timecontrol.current_timestep_size);

		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size);
		ts_psi2_rexi.setup(simVars.rexi, "psi2", simVars.timecontrol.current_timestep_size);

	}
	else if (timestepping_order == 4 && !use_only_linear_divergence)
	{
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

	if (simVars.disc.space_grid_use_c_staggering)
		FatalError("SWE_Plane_TS_l_rexi_na_sl_nd_etdrk: Staggering not supported for l_rexi_na_sl_nd_etdrk");

	//with_linear_div_only = i_use_linear_div;

	// Setup sampler for future interpolations
	sampler2D.setup(simVars.sim.plane_domain_size, op.planeDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(simVars.sim.plane_domain_size, op.planeDataConfig);


	PlaneData tmp_x(op.planeDataConfig);
	tmp_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)i)*simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0];
			},
			false
	);

	PlaneData tmp_y(op.planeDataConfig);
	tmp_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)j)*simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1];
			},
			false
	);

	// Initialize arrival points with h position
	ScalarDataArray pos_x = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_x);
	ScalarDataArray pos_y = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_y);


	double cell_size_x = simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0];
	double cell_size_y = simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1];

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

