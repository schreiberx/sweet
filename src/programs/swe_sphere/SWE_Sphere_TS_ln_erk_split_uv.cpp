/*
 * SWE_Sphere_TS_split_lg_lc_na_nr_erk.cpp
 *
 *  Created on: 24 Apr 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "SWE_Sphere_TS_ln_erk_split_uv.hpp"
#include <sweet/sphere/SphereData_DebugContainer.hpp>



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk_split_uv::euler_timestep_update_pert_lg(
		const SphereData_Spectral &i_U_phi_pert,	///< prognostic variables
		const SphereData_Spectral &i_U_vrt,	///< prognostic variables
		const SphereData_Spectral &i_U_div,	///< prognostic variables

		SphereData_Spectral &o_phi_pert_t,	///< time updates
		SphereData_Spectral &o_vrt_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation * simVars.sim.h0;


	const SphereData_Spectral &U_phi_pert = i_U_phi_pert;
	//const SphereData_Spectral &U_vrt = i_U_vrt;
	const SphereData_Spectral &U_div = i_U_div;


	SphereData_Physical U_u_phys, U_v_phys;
	op.vortdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	o_phi_pert_t -= gh0*U_div;
	o_div_t -= op.laplace(U_phi_pert);
}



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk_split_uv::euler_timestep_update_pert_lc(
		const SphereData_Spectral &i_U_phi,	///< prognostic variables
		const SphereData_Spectral &i_U_vrt,	///< prognostic variables
		const SphereData_Spectral &i_U_div,	///< prognostic variables

		SphereData_Spectral &o_phi_pert_t,	///< time updates
		SphereData_Spectral &o_vrt_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
//	double gh0 = simVars.sim.gravitation * simVars.sim.h0;


//	const SphereData_Spectral &U_phi_pert = i_U_phi;
	//const SphereData_Spectral &U_vrt = i_U_vrt;
//	const SphereData_Spectral &U_div = i_U_div;


	SphereData_Physical U_u_phys, U_v_phys;
	op.vortdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);


	SphereData_Physical fu_nl = op.fg*U_u_phys;
	SphereData_Physical fv_nl = op.fg*U_v_phys;

	SphereData_Spectral div, vrt;
	op.uv_to_vortdiv(fu_nl, fv_nl, vrt, div);

	o_vrt_t -= div;
	o_div_t += vrt;
}



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk_split_uv::euler_timestep_update_pert_na(
		const SphereData_Spectral &i_U_phi,	///< prognostic variables
		const SphereData_Spectral &i_U_vrt,	///< prognostic variables
		const SphereData_Spectral &i_U_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vrt_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
//	double gh0 = simVars.sim.gravitation * simVars.sim.h0;


	const SphereData_Spectral &U_phi_pert = i_U_phi;
	//const SphereData_Spectral &U_vrt = i_U_vrt;
	const SphereData_Spectral &U_div = i_U_div;


	SphereData_Physical U_u_phys, U_v_phys;
	op.vortdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);



	SphereData_Physical U_div_phys = U_div.toPhys();
	o_phi_t -= op.V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, U_phi_pert.toPhys());

	/*
	 * Velocity
	 */
	SphereData_Physical vrtg = i_U_vrt.toPhys();

	SphereData_Physical u_nl = U_u_phys*vrtg;
	SphereData_Physical v_nl = U_v_phys*vrtg;

	SphereData_Spectral vrt, div;
	op.uv_to_vortdiv(u_nl, v_nl, vrt, div);
	o_vrt_t -= div;
	o_div_t += vrt;

	o_div_t -= op.laplace(0.5*(U_u_phys*U_u_phys+U_v_phys*U_v_phys));
}



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk_split_uv::euler_timestep_update_pert_nr(
		const SphereData_Spectral &i_U_phi,	///< prognostic variables
		const SphereData_Spectral &i_U_vrt,	///< prognostic variables
		const SphereData_Spectral &i_U_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vrt_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	o_phi_t -= SphereData_Spectral(i_U_phi.toPhys()*i_U_div.toPhys());
}



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk_split_uv::euler_timestep_update_pert(
		const SphereData_Spectral &i_U_phi,	///< prognostic variables
		const SphereData_Spectral &i_U_vrt,	///< prognostic variables
		const SphereData_Spectral &i_U_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vrt_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
//	double gh0 = simVars.sim.gravitation * simVars.sim.h0;


//	const SphereData_Spectral &U_phi = i_U_phi;
	//const SphereData_Spectral &U_vrt = i_U_vrt;
//	const SphereData_Spectral &U_div = i_U_div;


	SphereData_Physical U_u_phys, U_v_phys;
	op.vortdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

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

			euler_timestep_update_pert_lg(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.getSphereDataPhysical();
			o_vrt_t += vrt_tmp.getSphereDataPhysical();
			o_div_t += div_tmp.getSphereDataPhysical();
		}


		if (use_lc)
		{
			phi_tmp.spectral_set_zero();
			vrt_tmp.spectral_set_zero();
			div_tmp.spectral_set_zero();

			euler_timestep_update_pert_lc(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.getSphereDataPhysical();
			o_vrt_t += vrt_tmp.getSphereDataPhysical();
			o_div_t += div_tmp.getSphereDataPhysical();
		}


		if (use_na)
		{
			phi_tmp.spectral_set_zero();
			vrt_tmp.spectral_set_zero();
			div_tmp.spectral_set_zero();

			euler_timestep_update_pert_na(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.getSphereDataPhysical();
			o_vrt_t += vrt_tmp.getSphereDataPhysical();
			o_div_t += div_tmp.getSphereDataPhysical();
		}


		if (use_nr)
		{
			phi_tmp.spectral_set_zero();
			vrt_tmp.spectral_set_zero();
			div_tmp.spectral_set_zero();

			euler_timestep_update_pert_nr(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.getSphereDataPhysical();
			o_vrt_t += vrt_tmp.getSphereDataPhysical();
			o_div_t += div_tmp.getSphereDataPhysical();
		}
	}
	else
	{

		/*
		 * See [SWEET]/doc/swe/swe_sphere_formulation/swe_on_sphere_formulation_in_sweet.pdf/lyx
		 */
		if (use_lg)
		{
			euler_timestep_update_pert_lg(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_lc)
		{
			euler_timestep_update_pert_lc(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_na)
		{
			euler_timestep_update_pert_na(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_nr)
		{
			euler_timestep_update_pert_nr(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}
	}

#if 0
	{
		/*
		 * Regular way how to compute it (the vort/div formulation of nonlinear advection part)
		 *
		 * However, there seems to go something wrong
		 */

		SphereData_Spectral reg_vrt, reg_div;
		/*
		 * Nonlinear parts of U/V
		 */
		SphereData_Physical vrtg = i_U_vrt.toPhys();

		SphereData_Spectral vrt, div;
		op.uv_to_vortdiv(U_u_phys*vrtg, U_v_phys*vrtg, vrt, div);
		reg_vrt = -div;
		reg_div = vrt;
		reg_div -= op.laplace(0.5*(U_u_phys*U_u_phys+U_v_phys*U_v_phys));


		//////////////////////////////////////////////////////////////////////////////////////////

		/*
		 * Computation based on first
		 * V \cdot grad(u)
		 * and
		 * V \cdot grad(v)
		 *
		 * then compute vort/div
		 */
		SphereData_Spectral new_vrt, new_div;

		SphereData_Physical U_div_phys = U_div.toPhys();

		SphereData_Spectral u_t = -op.V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, U_u_phys);
		SphereData_Spectral v_t = -op.V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, U_v_phys);

		op.uv_to_vortdiv(u_t.toPhys(), v_t.toPhys(), vrt, div);

		new_vrt = vrt;
		new_div = div;

		SphereData_DebugContainer::clear();
		SphereData_DebugContainer::append(reg_vrt - new_vrt, "tmp_vrt");
		SphereData_DebugContainer::append(reg_div - new_div, "tmp_div");
		SphereData_DebugContainer::append(reg_vrt, "reg_vrt");
		SphereData_DebugContainer::append(new_vrt, "new_vrt");
		SphereData_DebugContainer::append(reg_div, "reg_div");
		SphereData_DebugContainer::append(new_div, "new_div");
	}
#endif
}



void SWE_Sphere_TS_ln_erk_split_uv::run_timestep_pert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_ln_erk_split_uv::euler_timestep_update_pert,	///< pointer to function to compute euler time step updates
			io_phi, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Sphere_TS_ln_erk_split_uv::setup(
		int i_order,	///< order of RK time stepping method
		bool i_use_lg,
		bool i_use_lc,
		bool i_use_na,
		bool i_use_nr,

		bool i_antialiasing_for_each_term
)
{
	use_lg = i_use_lg;
	use_lc = i_use_lc;
	use_na = i_use_na;
	use_nr = i_use_nr;

	anti_aliasing_for_each_term = i_antialiasing_for_each_term;

	timestepping_order = i_order;
}



SWE_Sphere_TS_ln_erk_split_uv::SWE_Sphere_TS_ln_erk_split_uv(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
}



SWE_Sphere_TS_ln_erk_split_uv::~SWE_Sphere_TS_ln_erk_split_uv()
{
}

