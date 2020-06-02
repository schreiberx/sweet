/*
 * SWE_Sphere_TS_split_lg_lc_na_nr_erk.cpp
 *
 *  Created on: 24 Apr 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_ln_erk_split_vd.hpp"



void SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_update_lg(
		const SphereData_Spectral &i_U_phi,
		const SphereData_Spectral &i_U_vrt,
		const SphereData_Spectral &i_U_div,

		SphereData_Spectral &o_phi_t,
		SphereData_Spectral &o_vrt_t,
		SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation*simVars.sim.h0;

	o_phi_t -= gh0*i_U_div;
	o_div_t -= op.laplace(i_U_phi);
}



void SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_update_lc(
		const SphereData_Spectral &i_U_phi,
		const SphereData_Spectral &i_U_vrt,
		const SphereData_Spectral &i_U_div,

		SphereData_Spectral &o_phi_t,
		SphereData_Spectral &o_vrt_t,
		SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	SphereData_Physical U_u_phys, U_v_phys;
	op.vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	// dt calculation starts here

	SphereData_Physical fu_nl = op.fg*U_u_phys;
	SphereData_Physical fv_nl = op.fg*U_v_phys;

	SphereData_Spectral div, vrt;
	op.uv_to_vrtdiv(fu_nl, fv_nl, vrt, div);

	o_vrt_t -= div;
	o_div_t += vrt;
}



void SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_update_na(
		const SphereData_Spectral &i_U_phi,
		const SphereData_Spectral &i_U_vrt,
		const SphereData_Spectral &i_U_div,

		SphereData_Spectral &o_phi_t,
		SphereData_Spectral &o_vrt_t,
		SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	SphereData_Physical U_u_phys, U_v_phys;
	op.vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	SphereData_Physical U_div_phys = i_U_div.toPhys();
	o_phi_t -= op.V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_phi.toPhys());
	o_vrt_t -= op.V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_vrt.toPhys());
	o_div_t -= op.V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U_div.toPhys());
}



void SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_update_nr(
		const SphereData_Spectral &i_U_phi,
		const SphereData_Spectral &i_U_vrt,
		const SphereData_Spectral &i_U_div,

		SphereData_Spectral &o_phi_t,
		SphereData_Spectral &o_vrt_t,
		SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	SphereData_Physical U_u_phys, U_v_phys;
	op.vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	// dt calculation starts here

	SphereData_Physical U_div_phys = i_U_div.toPhys();

	o_phi_t -= SphereData_Spectral(i_U_phi.toPhys()*i_U_div.toPhys());

	if (0)
	{
		o_vrt_t -= o_vrt_t.toPhys()*U_div_phys;
	}
	else
	{

		/*
		 * N from UV formulation
		 */
//		double gh0 = simVars.sim.gravitation * simVars.sim.h0;


//		const SphereData_Spectral &U_phi = i_U_phi;
		//const SphereData_Spectral &U_vrt = i_U_vrt;
		const SphereData_Spectral &U_div = i_U_div;


		SphereData_Physical U_u_phys, U_v_phys;
		op.vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

		SphereData_Physical U_div_phys = U_div.toPhys();

		/*
		 * Velocity
		 */
		SphereData_Physical vrtg = i_U_vrt.toPhys();

		SphereData_Physical u_nl = U_u_phys*vrtg;
		SphereData_Physical v_nl = U_v_phys*vrtg;

		SphereData_Spectral vrt, div;
		op.uv_to_vrtdiv(u_nl, v_nl, vrt, div);
		//o_vrt_t -= div;


		/*
		 * NA part to be subtracted
		 */
		SphereData_Spectral phi_tmp(i_U_phi.sphereDataConfig);
		SphereData_Spectral vrt_tmp(i_U_vrt.sphereDataConfig);
		SphereData_Spectral div_tmp(i_U_div.sphereDataConfig);

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

	const SphereData_Physical U_vrt_phys = i_U_vrt.toPhys();
	o_div_t += op.uv_to_vort(U_vrt_phys*U_u_phys, U_vrt_phys*U_v_phys);
	o_div_t += op.uv_to_div(U_div_phys*U_u_phys, U_div_phys*U_v_phys);
	o_div_t -= 0.5*op.laplace(U_u_phys*U_u_phys + U_v_phys*U_v_phys);
	o_div_t -= U_div_phys*U_div_phys;
}



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_set_tendencies(
		const SphereData_Spectral &i_U_phi,
		const SphereData_Spectral &i_U_vrt,
		const SphereData_Spectral &i_U_div,

		SphereData_Spectral &o_phi_t,
		SphereData_Spectral &o_vrt_t,
		SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	SphereData_Physical U_u_phys, U_v_phys;
	op.vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	o_phi_t.spectral_set_zero();
	o_vrt_t.spectral_set_zero();
	o_div_t.spectral_set_zero();

	if (anti_aliasing_for_each_term)
	{
		SphereData_Spectral phi_tmp(i_U_phi.sphereDataConfig);
		SphereData_Spectral vrt_tmp(i_U_vrt.sphereDataConfig);
		SphereData_Spectral div_tmp(i_U_div.sphereDataConfig);


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


		SphereData_Spectral vrt_backup = o_vrt_t;
		SphereData_Spectral div_backup = o_div_t;

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
			SphereData_Spectral vrt_n(i_U_vrt.sphereDataConfig, 0);
			SphereData_Spectral div_n(i_U_div.sphereDataConfig, 0);

//			double gh0 = simVars.sim.gravitation * simVars.sim.h0;


			SphereData_Physical U_u_phys, U_v_phys;
			op.vrtdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

			/*
			 * Velocity
			 */
			SphereData_Physical vrtg = i_U_vrt.toPhys();

			SphereData_Physical u_nl = U_u_phys*vrtg;
			SphereData_Physical v_nl = U_v_phys*vrtg;

			SphereData_Spectral vrt, div;
			op.uv_to_vrtdiv(u_nl, v_nl, vrt, div);
			vrt_n -= div;

			div_n += vrt;
			div_n -= op.laplace(0.5*(U_u_phys*U_u_phys+U_v_phys*U_v_phys));

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



void SWE_Sphere_TS_ln_erk_split_vd::run_timestep(
		SphereData_Spectral &io_phi,
		SphereData_Spectral &io_vrt,
		SphereData_Spectral &io_div,

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_set_tendencies,	///< pointer to function to compute euler time step updates
			io_phi, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}


/*
 * Setup
 */
void SWE_Sphere_TS_ln_erk_split_vd::setup(
		int i_order,		///< order of RK time stepping method

		bool i_use_lg,
		bool i_use_lc,
		bool i_use_na,
		bool i_use_nr,

		bool i_antialiasing_for_each_term
)
{
	timestepping_order = i_order;

	use_lg = i_use_lg;
	use_lc = i_use_lc;
	use_na = i_use_na;
	use_nr = i_use_nr;

	anti_aliasing_for_each_term = i_antialiasing_for_each_term;
}



SWE_Sphere_TS_ln_erk_split_vd::SWE_Sphere_TS_ln_erk_split_vd(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
}



SWE_Sphere_TS_ln_erk_split_vd::~SWE_Sphere_TS_ln_erk_split_vd()
{
}

