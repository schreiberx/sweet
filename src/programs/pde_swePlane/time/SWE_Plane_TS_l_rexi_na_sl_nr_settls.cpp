/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_settls.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */


#include "SWE_Plane_TS_l_rexi_na_sl_nr_settls.hpp"


void SWE_Plane_TS_l_rexi_na_sl_nr_settls::runTimestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables - perturbation of height!
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables - zonal velocity
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables - meridional velocity

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_rexi_na_sl_nd_settls: Only constant time step size allowed (Please set --dt)");

	const sweet::PlaneData_Config *planeDataConfig = io_h.planeDataConfig;

	//Out vars
	sweet::PlaneData_Spectral h(io_h.planeDataConfig);
	sweet::PlaneData_Spectral u(io_h.planeDataConfig);
	sweet::PlaneData_Spectral v(io_h.planeDataConfig);

	//Temporary
	sweet::PlaneData_Spectral N_h(io_h.planeDataConfig);
	sweet::PlaneData_Spectral N_u(io_h.planeDataConfig);
	sweet::PlaneData_Spectral N_v(io_h.planeDataConfig);
	sweet::PlaneData_Spectral hdiv(io_h.planeDataConfig);
	sweet::PlaneData_Spectral N_h_prev(io_h.planeDataConfig);
	sweet::PlaneData_Spectral N_h_ext(io_h.planeDataConfig);

	// Departure points
	sweet::ScalarDataArray posx_d(io_h.planeDataConfig->physical_array_data_number_of_elements);
	sweet::ScalarDataArray posy_d(io_h.planeDataConfig->physical_array_data_number_of_elements);

	// Parameters
	double dt = i_dt;

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

	// Calculate departure points
	//Calculate departure points - always force to be second order accurate!
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev.toPhys(),	v_prev.toPhys(),
			u.toPhys(),		v.toPhys(),
			posx_a,		posy_a,
			dt,
			posx_d,	posy_d,			// output
			shackPlaneDataOps->plane_domain_size,
			&staggering,
			2, //shackDict.disc.timestepping_order

			shackPDESWETimeDisc->semi_lagrangian_max_iterations,
			shackPDESWETimeDisc->semi_lagrangian_convergence_threshold
	);

	N_u.spectral_set_zero();
	N_v.spectral_set_zero();
	N_h.spectral_set_zero();
	N_h_prev.spectral_set_zero();
	N_h_ext.spectral_set_zero();
	hdiv.spectral_set_zero();

	//Original more stable scheme

	//Calculate nonlinear terms - not done in case of only linear divergence (linear div is already in linear part)
	if (!use_only_linear_divergence) // Full nonlinear case
	{
		// Calculate nonlinear term for the previous time step
		N_h = -h_prev * (ops->diff_c_x(u_prev) + ops->diff_c_y(v_prev));

		if (shackPDESWEPlane->use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
			SWEETError("Implicit diffusion only supported with spectral space activated");
#else
			N_h = ops->implicit_diffusion(N_h, shackTimestepControl->current_timestep_size*shackPDESWEPlane->viscosity, shackPDESWEPlane->viscosity_order);
#endif
		}

		//Calculate exp(Ldt)N(n-1), relative to previous timestep
		ts_l_rexi.runTimestep(N_h, N_u, N_v, i_dt, i_simulation_timestamp);

		//Use N_h to store now the nonlinearity of the current time (prev will not be required anymore)
		//Update the nonlinear terms with the constants relative to dt
		// N=dtN^n-0.5dt exp(dtL)N^n-1 from paper
		// N=-h*div is calculate in cartesian space (pseudo-spectrally)
		hdiv =  - h * (ops->diff_c_x(io_u) + ops->diff_c_y(io_v));
		if (shackPDESWEPlane->use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
			SWEETError("Implicit diffusion only supported with spectral space activated");
#else
			hdiv = ops->implicit_diffusion(hdiv, shackTimestepControl->current_timestep_size*shackPDESWEPlane->viscosity, shackPDESWEPlane->viscosity_order);
#endif
		}

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
	if (!use_only_linear_divergence) // Full nonlinear case
	{
		h = h + 0.5 * dt * hdiv;
	}

	/*
	 * Calculate the exp(Ldt) of the resulting u,v,h
	 */
	//Calculate phi_0 of interpolated U
	sweet::PlaneData_Spectral phi0_Un_h(planeDataConfig);
	sweet::PlaneData_Spectral phi0_Un_u(planeDataConfig);
	sweet::PlaneData_Spectral phi0_Un_v(planeDataConfig);
	//ts_l_rexi.runTimestep(h, u, v, i_dt, i_simulation_timestamp);
	ts_l_rexi.runTimestep(
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
bool SWE_Plane_TS_l_rexi_na_sl_nr_settls::setup(
		sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	h_prev.setup(ops->planeDataConfig);
	u_prev.setup(ops->planeDataConfig);
	v_prev.setup(ops->planeDataConfig);

	posx_a.setup(ops->planeDataConfig->physical_array_data_number_of_elements);
	posy_a.setup(ops->planeDataConfig->physical_array_data_number_of_elements);

	posx_d.setup(ops->planeDataConfig->physical_array_data_number_of_elements);
	posy_d.setup(ops->planeDataConfig->physical_array_data_number_of_elements);


	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	ts_l_rexi.setup(io_ops, "phi0");

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("SWE_Plane_TS_l_rexi_na_sl_nd_settls: Staggering not supported for l_rexi_na_sl_nd_settls");

	//with_linear_div_only = i_use_linear_div;

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

	return true;
}



bool SWE_Plane_TS_l_rexi_na_sl_nr_settls::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	PDESWEPlaneTS_BaseInterface::shackRegistration(io_shackDict);

	ts_l_rexi.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_RETURN_BOOLEAN(ts_l_rexi);
	return true;
}
