/*
 *  Created on: 09 Oct 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *  Changelog:
 *      based on Martin Schreiber ETD timestepper
 */

#include "SWE_Plane_TS_l_rexi_na_sl_nr_etdrk.hpp"

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_rexi_na_sl_nr_etdrk::euler_timestep_update_nonlinear(
		const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

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
	 * o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	 * o_u_t = -i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
	 * o_v_t = -i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);
	 */
	// In Lagrangian form, the only nonlinearity is the nonlinear divergence
	o_u_t.spectral_set_zero(); //-i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
	o_v_t.spectral_set_zero(); // = 0.0; //-i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);

	// linear div only
	if (use_only_linear_divergence)
	{
		o_h_t.spectral_set_zero(); // = 0.0; //-ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	}
	else
	{
		// nonlinear div
		o_h_t = -i_h*(ops->diff_c_x(i_u) + ops->diff_c_y(i_v));
		// Smooth spectrum to avoid instability
		if (shackPDESWEPlane->use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
			SWEETError("Implicit diffusion only supported with spectral space activated");
#else
			o_h_t= ops->implicit_diffusion(o_h_t, shackTimestepControl->current_timestep_size*shackPDESWEPlane->viscosity, shackPDESWEPlane->viscosity_order);
#endif
		}

	}
}



void SWE_Plane_TS_l_rexi_na_sl_nr_etdrk::run_timestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_rexi_na_sl_nd_etdrk: Only constant time step size allowed (Please set --dt)");


	const sweet::PlaneDataConfig *planeDataConfig = io_h.planeDataConfig;

	// Tmp vars
	//h, u, v tmp
	sweet::PlaneData_Spectral h(planeDataConfig);
	sweet::PlaneData_Spectral u(planeDataConfig);
	sweet::PlaneData_Spectral v(planeDataConfig);
	//Nonlinear calculation of u,v,h from input
	sweet::PlaneData_Spectral FUn_h(planeDataConfig);
	sweet::PlaneData_Spectral FUn_u(planeDataConfig);
	sweet::PlaneData_Spectral FUn_v(planeDataConfig);


	//Departure points and arrival points
	sweet::ScalarDataArray posx_d(io_h.planeDataConfig->physical_array_data_number_of_elements);
	sweet::ScalarDataArray posy_d(io_h.planeDataConfig->physical_array_data_number_of_elements);


	Staggering staggering;
	assert(staggering.staggering_type == 'a');

	if (i_simulation_timestamp == 0)
	{
#if (!SWEET_PARAREAL) && (!SWEET_XBRAID)
		/*
		 * First time step
		 */
		h_prev = io_h;
		u_prev = io_u;
		v_prev = io_v;
#endif
	}


	//Preserve io unmodified
	u = io_u;
	v = io_v;
	h = io_h;

	// Calculate departure points - always force to be second order accurate!
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev.toPhys(),	v_prev.toPhys(),
			u.toPhys(),		v.toPhys(),
			posx_a,		posy_a,
			i_dt,
			posx_d,	posy_d,			// output
			shackPlaneDataOps->plane_domain_size,
			&staggering,
			2, //simVars.disc.timestepping_order,

			shackPDESWETimeDisc->semi_lagrangian_max_iterations,
			shackPDESWETimeDisc->semi_lagrangian_convergence_threshold
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
			sweet::PlaneData_Spectral psi1_FUn_h(planeDataConfig);
			sweet::PlaneData_Spectral psi1_FUn_u(planeDataConfig);
			sweet::PlaneData_Spectral psi1_FUn_v(planeDataConfig);

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
		sweet::PlaneData_Spectral phi0_Un_h(planeDataConfig);
		sweet::PlaneData_Spectral phi0_Un_u(planeDataConfig);
		sweet::PlaneData_Spectral phi0_Un_v(planeDataConfig);
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
		SWEETError("TODO: This order is not implemented, yet!");
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
		sweet::PlaneData_Spectral h1(planeDataConfig);
		sweet::PlaneData_Spectral u1(planeDataConfig);
		sweet::PlaneData_Spectral v1(planeDataConfig);
		h1 = h;
		u1 = u;
		v1 = v;

		//Calculate psi2NU_1
		//-----------------------

		//NU_1
		sweet::PlaneData_Spectral FU1_h(planeDataConfig);
		sweet::PlaneData_Spectral FU1_u(planeDataConfig);
		sweet::PlaneData_Spectral FU1_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				h1, u1, v1,
				FU1_h, FU1_u, FU1_v,
				i_simulation_timestamp
		);


		//Apply psi2
		sweet::PlaneData_Spectral psi2_FU1_h(planeDataConfig);
		sweet::PlaneData_Spectral psi2_FU1_u(planeDataConfig);
		sweet::PlaneData_Spectral psi2_FU1_v(planeDataConfig);
		ts_psi2_rexi.run_timestep(
				FU1_h, FU1_u, FU1_v,
				psi2_FU1_h, psi2_FU1_u, psi2_FU1_v,
				i_dt,
				i_simulation_timestamp
		);

		//Calculate psi2NUn
		sweet::PlaneData_Spectral psi2_FUn_h(planeDataConfig);
		sweet::PlaneData_Spectral psi2_FUn_u(planeDataConfig);
		sweet::PlaneData_Spectral psi2_FUn_v(planeDataConfig);

		ts_psi2_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				psi2_FUn_h, psi2_FUn_u, psi2_FUn_v,
				i_dt,
				i_simulation_timestamp
		);


		//Interpolate psi2NUn to departure points ( )_*

		sweet::PlaneData_Spectral psi2FUn_h_dep(planeDataConfig);
		sweet::PlaneData_Spectral psi2FUn_u_dep(planeDataConfig);
		sweet::PlaneData_Spectral psi2FUn_v_dep(planeDataConfig);
		psi2FUn_h_dep = sampler2D.bicubic_scalar(psi2_FUn_h, posx_d, posy_d, -0.5, -0.5);
		psi2FUn_u_dep = sampler2D.bicubic_scalar(psi2_FUn_u, posx_d, posy_d, -0.5, -0.5);
		psi2FUn_v_dep = sampler2D.bicubic_scalar(psi2_FUn_v, posx_d, posy_d, -0.5, -0.5);


		//psi2NU_1-psi2NUn_dep
		sweet::PlaneData_Spectral dif2_h(io_h.planeDataConfig);
		sweet::PlaneData_Spectral dif2_u(io_u.planeDataConfig);
		sweet::PlaneData_Spectral dif2_v(io_v.planeDataConfig);

		dif2_h = psi2_FU1_h-psi2FUn_h_dep;
		dif2_u = psi2_FU1_u-psi2FUn_u_dep;
		dif2_v = psi2_FU1_v-psi2FUn_v_dep;

		//Apply phi0
		//second order final forcing
		sweet::PlaneData_Spectral phi0_dif2_h(planeDataConfig);
		sweet::PlaneData_Spectral phi0_dif2_u(planeDataConfig);
		sweet::PlaneData_Spectral phi0_dif2_v(planeDataConfig);
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




bool SWE_Plane_TS_l_rexi_na_sl_nr_etdrk::registerShacks(
		sweet::ShackDictionary *io_shackDict
)
{
	PDESWEPlaneTS_BaseInterface::registerShacks(io_shackDict);

	ts_phi0_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_phi0_rexi);

	ts_phi1_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_phi1_rexi);

	ts_phi2_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_phi2_rexi);


	ts_ups0_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_ups0_rexi);

	ts_ups1_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_ups1_rexi);

	ts_ups2_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_ups2_rexi);

	ts_ups3_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_ups3_rexi);


	ts_psi1_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_psi1_rexi);

	ts_psi2_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_psi2_rexi);

	ts_psi3_rexi.registerShacks(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_psi3_rexi);

	return true;
}


/**
 * Setup
 */
bool SWE_Plane_TS_l_rexi_na_sl_nr_etdrk::setup(
		sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;


	ts_phi0_rexi.setup(io_ops, "phi0");

	if (timestepping_order == 1 && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(io_ops, "phi1");
		ts_psi1_rexi.setup(io_ops, "psi1");
	}
	else if (timestepping_order == 2 && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(io_ops, "phi1");
		ts_phi2_rexi.setup(io_ops, "phi2");

		ts_psi1_rexi.setup(io_ops, "psi1");
		ts_psi2_rexi.setup(io_ops, "psi2");
	}
	else if (timestepping_order == 4 && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(io_ops, "phi1");
		ts_phi2_rexi.setup(io_ops, "phi2");

		ts_psi1_rexi.setup(io_ops, "psi1");
		ts_psi2_rexi.setup(io_ops, "psi2");
		ts_psi3_rexi.setup(io_ops, "psi3");

		ts_ups0_rexi.setup(io_ops, "phi0");
		ts_ups1_rexi.setup(io_ops, "ups1");
		ts_ups2_rexi.setup(io_ops, "ups2");
		ts_ups3_rexi.setup(io_ops, "ups3");
	}

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("SWE_Plane_TS_l_rexi_na_sl_nd_etdrk: Staggering not supported for l_rexi_na_sl_nd_etdrk");

	// Setup sampler for future interpolations
	sampler2D.setup(shackPlaneDataOps->plane_domain_size, ops->planeDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(shackPlaneDataOps->plane_domain_size, ops->planeDataConfig);


	sweet::PlaneData_Physical tmp_x(ops->planeDataConfig);
	tmp_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)i)*shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0];
			},
			false
	);

	sweet::PlaneData_Physical tmp_y(ops->planeDataConfig);
	tmp_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)j)*shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1];
			},
			false
	);

	// Initialize arrival points with h position
	sweet::ScalarDataArray pos_x = sweet::Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_x);
	sweet::ScalarDataArray pos_y = sweet::Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_y);


	double cell_size_x = shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0];
	double cell_size_y = shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1];

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*cell_size_x;
	posy_a = pos_y+0.5*cell_size_y;


	h_prev.setup(ops->planeDataConfig);
	u_prev.setup(ops->planeDataConfig);
	v_prev.setup(ops->planeDataConfig);

	posx_a.setup(ops->planeDataConfig->physical_array_data_number_of_elements);
	posy_a.setup(ops->planeDataConfig->physical_array_data_number_of_elements);

	posx_d.setup(ops->planeDataConfig->physical_array_data_number_of_elements);
	posy_d.setup(ops->planeDataConfig->physical_array_data_number_of_elements);

	return true;
}

