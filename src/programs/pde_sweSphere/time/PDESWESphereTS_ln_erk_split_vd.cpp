/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */



#include "PDESWESphereTS_ln_erk_split_vd.hpp"


bool PDESWESphereTS_ln_erk_split_vd::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	timestepping_order = shackPDESWETimeDisc->timestepping_order;

	/*
	 * l_na
	 */
	if (timestepping_method == "l_na_erk_split_vd")
		return setup_main(io_ops, timestepping_order, true, true, true, false, false);

	if (timestepping_method == "l_na_erk_split_aa_vd")
		return setup_main(io_ops, timestepping_order, true, true, true, false, true);

	/*
	 * l
	 */
	if (timestepping_method == "l_erk_split_vd")
		return setup_main(io_ops, timestepping_order, true, true, false, false, false);

	if (timestepping_method == "l_erk_split_aa_vd")
		return setup_main(io_ops, timestepping_order, true, true, false, false, true);

	/*
	 * ln
	 */
	if (timestepping_method == "ln_erk_split_vd")
		return setup_main(io_ops, timestepping_order, true, true, true, true, false);

	if (timestepping_method == "ln_erk_split_aa_vd")
		return setup_main(io_ops, timestepping_order, true, true, true, true, true);

	SWEETError("Should never happen");
	return false;
}


bool PDESWESphereTS_ln_erk_split_vd::setup_main(
		sweet::SphereOperators *io_ops,
		int i_order,		///< order of RK time stepping method
		bool i_use_lg,
		bool i_use_lc,
		bool i_use_na,
		bool i_use_nr,

		bool i_antialiasing_for_each_term
)
{
	ops = io_ops;

	setupFG();

	timestepping_order = i_order;

	use_lg = i_use_lg;
	use_lc = i_use_lc;
	use_na = i_use_na;
	use_nr = i_use_nr;

	anti_aliasing_for_each_term = i_antialiasing_for_each_term;
	return true;
}



void PDESWESphereTS_ln_erk_split_vd::euler_timestep_update_lg(
		const sweet::SphereData_Spectral &i_U_phi,
		const sweet::SphereData_Spectral &i_U_vrt,
		const sweet::SphereData_Spectral &i_U_div,

		sweet::SphereData_Spectral &o_phi_t,
		sweet::SphereData_Spectral &o_vrt_t,
		sweet::SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	double gh0 = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

	o_phi_t -= gh0*i_U_div;
	o_div_t -= ops->laplace(i_U_phi);
}



void PDESWESphereTS_ln_erk_split_vd::euler_timestep_update_lc(
		const sweet::SphereData_Spectral &i_U_phi,
		const sweet::SphereData_Spectral &i_U_vrt,
		const sweet::SphereData_Spectral &i_U_div,

		sweet::SphereData_Spectral &io_phi_t,
		sweet::SphereData_Spectral &io_vrt_t,
		sweet::SphereData_Spectral &io_div_t,

		double i_simulation_timestamp
)
{
	sweet::SphereData_Physical U_u_phys, U_v_phys;
	ops->vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	// Apply f term
	sweet::SphereData_Physical fu_nl = fg*U_u_phys;
	sweet::SphereData_Physical fv_nl = fg*U_v_phys;

	sweet::SphereData_Spectral div, vrt;
	ops->uv_to_vrtdiv(fu_nl, fv_nl, vrt, div);

	io_vrt_t -= div;
	io_div_t += vrt;
}



/*
 * This is a version which only operates in spectral space.
 *
 * It doesn't necessarily reflect exactly 1:1 the application of the Coriolis effect in physical space.
 */
void PDESWESphereTS_ln_erk_split_vd::euler_timestep_update_lc_spectral_only(
		const sweet::SphereData_Spectral &i_U_phi,
		const sweet::SphereData_Spectral &i_U_vrt,
		const sweet::SphereData_Spectral &i_U_div,

		sweet::SphereData_Spectral &io_phi_t,
		sweet::SphereData_Spectral &io_vrt_t,
		sweet::SphereData_Spectral &io_div_t,

		double i_simulation_timestamp
)
{
	io_vrt_t -= ops->implicit_F(i_U_div, 2.0*shackPDESWESphere->sphere_rotating_coriolis_omega);
	io_div_t += ops->implicit_F(i_U_vrt, 2.0*shackPDESWESphere->sphere_rotating_coriolis_omega);
}



void PDESWESphereTS_ln_erk_split_vd::euler_timestep_update_na(
		const sweet::SphereData_Spectral &i_U_phi,
		const sweet::SphereData_Spectral &i_U_vrt,
		const sweet::SphereData_Spectral &i_U_div,

		sweet::SphereData_Spectral &o_phi_t,
		sweet::SphereData_Spectral &o_vrt_t,
		sweet::SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	sweet::SphereData_Physical U_u_phys, U_v_phys;
	ops->vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	sweet::SphereData_Physical U_div_phys = i_U_div.toPhys();
	o_phi_t -= ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_phi.toPhys());
	o_vrt_t -= ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_vrt.toPhys());
	o_div_t -= ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_div.toPhys());
}



void PDESWESphereTS_ln_erk_split_vd::euler_timestep_update_nr(
		const sweet::SphereData_Spectral &i_U_phi,
		const sweet::SphereData_Spectral &i_U_vrt,
		const sweet::SphereData_Spectral &i_U_div,

		sweet::SphereData_Spectral &o_phi_t,
		sweet::SphereData_Spectral &o_vrt_t,
		sweet::SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	sweet::SphereData_Physical U_u_phys, U_v_phys;
	ops->vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	// dt calculation starts here

	sweet::SphereData_Physical U_div_phys = i_U_div.toPhys();

	o_phi_t -= sweet::SphereData_Spectral(i_U_phi.toPhys()*i_U_div.toPhys());

	if (0)
	{
		o_vrt_t -= o_vrt_t.toPhys()*U_div_phys;
	}
	else
	{

		/*
		 * N from UV formulation
		 */
//		double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;


//		const sweet::SphereData_Spectral &U_phi = i_U_phi;
		//const sweet::SphereData_Spectral &U_vrt = i_U_vrt;
		const sweet::SphereData_Spectral &U_div = i_U_div;


		sweet::SphereData_Physical U_u_phys, U_v_phys;
		ops->vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

		sweet::SphereData_Physical U_div_phys = U_div.toPhys();

		/*
		 * Velocity
		 */
		sweet::SphereData_Physical vrtg = i_U_vrt.toPhys();

		sweet::SphereData_Physical u_nl = U_u_phys*vrtg;
		sweet::SphereData_Physical v_nl = U_v_phys*vrtg;

		sweet::SphereData_Spectral vrt, div;
		ops->uv_to_vrtdiv(u_nl, v_nl, vrt, div);
		//o_vrt_t -= div;


		/*
		 * NA part to be subtracted
		 */
		sweet::SphereData_Spectral phi_tmp(i_U_phi.sphereDataConfig);
		sweet::SphereData_Spectral vrt_tmp(i_U_vrt.sphereDataConfig);
		sweet::SphereData_Spectral div_tmp(i_U_div.sphereDataConfig);

		phi_tmp.spectral_set_zero();
		vrt_tmp.spectral_set_zero();
		div_tmp.spectral_set_zero();

		euler_timestep_update_na(
				i_U_phi, i_U_vrt, i_U_div,
				phi_tmp, vrt_tmp, div_tmp,
				i_simulation_timestamp
			);

		vrt_tmp = vrt_tmp.toPhys();

		o_vrt_t += -div - vrt_tmp;
	}

	const sweet::SphereData_Physical U_vrt_phys = i_U_vrt.toPhys();
	o_div_t += ops->uv_to_vort(U_vrt_phys*U_u_phys, U_vrt_phys*U_v_phys);
	o_div_t += ops->uv_to_div(U_div_phys*U_u_phys, U_div_phys*U_v_phys);
	o_div_t -= 0.5*ops->laplace(U_u_phys*U_u_phys + U_v_phys*U_v_phys);
	o_div_t -= U_div_phys*U_div_phys;
}



/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphereTS_ln_erk_split_vd::euler_timestep_set_tendencies(
		const sweet::SphereData_Spectral &i_U_phi,
		const sweet::SphereData_Spectral &i_U_vrt,
		const sweet::SphereData_Spectral &i_U_div,

		sweet::SphereData_Spectral &o_phi_t,
		sweet::SphereData_Spectral &o_vrt_t,
		sweet::SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	sweet::SphereData_Physical U_u_phys, U_v_phys;
	ops->vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	o_phi_t.spectral_set_zero();
	o_vrt_t.spectral_set_zero();
	o_div_t.spectral_set_zero();

	if (anti_aliasing_for_each_term)
	{
		sweet::SphereData_Spectral phi_tmp(i_U_phi.sphereDataConfig);
		sweet::SphereData_Spectral vrt_tmp(i_U_vrt.sphereDataConfig);
		sweet::SphereData_Spectral div_tmp(i_U_div.sphereDataConfig);


		/*
		 * See [SWEET]/doc/swe/swe_sphere_formulation/swe_on_sphere_formulation_in_sweet.pdf/lyx
		 */
		if (use_lg)
		{
			phi_tmp.spectral_set_zero();
			vrt_tmp.spectral_set_zero();
			div_tmp.spectral_set_zero();

			euler_timestep_update_lg(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.toPhys();
			o_vrt_t += vrt_tmp.toPhys();
			o_div_t += div_tmp.toPhys();
		}


		if (use_lc)
		{
			phi_tmp.spectral_set_zero();
			vrt_tmp.spectral_set_zero();
			div_tmp.spectral_set_zero();

			euler_timestep_update_lc(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.toPhys();
			o_vrt_t += vrt_tmp.toPhys();
			o_div_t += div_tmp.toPhys();
		}


		sweet::SphereData_Spectral vrt_backup = o_vrt_t;
		sweet::SphereData_Spectral div_backup = o_div_t;

		if (use_na)
		{
			phi_tmp.spectral_set_zero();
			vrt_tmp.spectral_set_zero();
			div_tmp.spectral_set_zero();

			euler_timestep_update_na(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp
				);

			o_phi_t += phi_tmp.toPhys();
			o_vrt_t += vrt_tmp.toPhys();
			o_div_t += div_tmp.toPhys();
		}


		if (use_nr)
		{
			phi_tmp.spectral_set_zero();
			vrt_tmp.spectral_set_zero();
			div_tmp.spectral_set_zero();

			euler_timestep_update_nr(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.toPhys();
			o_vrt_t += vrt_tmp.toPhys();
			o_div_t += div_tmp.toPhys();
		}

#if 0
		{
			sweet::SphereData_Spectral vrt_n(i_U_vrt.sphereDataConfig, 0);
			sweet::SphereData_Spectral div_n(i_U_div.sphereDataConfig, 0);

//			double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;


			sweet::SphereData_Physical U_u_phys, U_v_phys;
			ops->vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

			/*
			 * Velocity
			 */
			sweet::SphereData_Physical vrtg = i_U_vrt.toPhys();

			sweet::SphereData_Physical u_nl = U_u_phys*vrtg;
			sweet::SphereData_Physical v_nl = U_v_phys*vrtg;

			sweet::SphereData_Spectral vrt, div;
			ops->uv_to_vrtdiv(u_nl, v_nl, vrt, div);
			vrt_n -= div;

			div_n += vrt;
			div_n -= ops->laplace(0.5*(U_u_phys*U_u_phys+U_v_phys*U_v_phys));

			o_vrt_t = vrt_backup + vrt_n;
			//o_div_t = div_backup + div_n;
		}
#endif
	}
	else
	{

		/*
		 * See [SWEET]/doc/swe/swe_sphere_formulation/swe_on_sphere_formulation_in_sweet.pdf/lyx
		 */
		if (use_lg)
		{
			euler_timestep_update_lg(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_lc)
		{
			euler_timestep_update_lc(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_na)
		{
			euler_timestep_update_na(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_nr)
		{
			euler_timestep_update_nr(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}
	}
}



void PDESWESphereTS_ln_erk_split_vd::runTimestep(
		sweet::SphereData_Spectral &io_phi,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDESWESphereTS_ln_erk_split_vd::euler_timestep_set_tendencies,	///< pointer to function to compute euler time step updates
			io_phi, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



PDESWESphereTS_ln_erk_split_vd::PDESWESphereTS_ln_erk_split_vd()
{
}



PDESWESphereTS_ln_erk_split_vd::~PDESWESphereTS_ln_erk_split_vd()
{
}

