/*
 * SWE_Sphere_TS_split_lg_lc_na_nr_erk.cpp
 *
 *  Created on: 24 Apr 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "SWE_Sphere_TS_ln_erk_split_vd.hpp"
#include <sweet/sphere/SphereData_DebugContainer.hpp>



/*
 * Compute Nonlinear advection terms
 *
 * U \cdot \grad phi = \div \cdot (V*phi) - \nabla
 */
SphereData_Spectral SWE_Sphere_TS_ln_erk_split_vd::V_dot_grad_scalar(
		const SphereData_Physical &i_u_phys,		///< u velocity
		const SphereData_Physical &i_v_phys,		///< v velocity
		const SphereData_Physical &i_div_phys,		///< divergence in physical space to avoid transformation
		const SphereData_Physical &i_scalar_phys	///< scalar
)
{
	return op.uv_to_div(
			i_u_phys*i_scalar_phys,
			i_v_phys*i_scalar_phys
		)
			- i_div_phys*i_scalar_phys;
}



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_update_pert(
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
	const SphereData_Config *sphereDataConfig = i_U_phi_pert.sphereDataConfig;


	const SphereData_Spectral &U_phi_pert = i_U_phi_pert;
	const SphereData_Spectral &U_vrt = i_U_vrt;
	const SphereData_Spectral &U_div = i_U_div;


	SphereData_Physical U_u_phys, U_v_phys;
	op.vortdiv_to_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	o_phi_pert_t.spectral_set_zero();
	o_vrt_t.spectral_set_zero();
	o_div_t.spectral_set_zero();


	/*
	 * See [SWEET]/doc/swe/swe_sphere_formulation/swe_on_sphere_formulation_in_sweet.pdf/lyx
	 */

	if (use_lg)
	{
		o_phi_pert_t -= gh0*U_div;
		o_div_t -= op.laplace(U_phi_pert);
	}


	if (use_lc)
	{
		SphereData_Physical fu_nl = op.fg*U_u_phys;
		SphereData_Physical fv_nl = op.fg*U_v_phys;

		SphereData_Spectral div, vrt;
		op.uv_to_vortdiv(fu_nl, fv_nl, vrt, div);

		o_vrt_t -= div;
		o_div_t += vrt;
	}


	if (use_na)
	{
		SphereData_Physical U_div_phys = U_div.toPhys();
		o_phi_pert_t -= V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, U_phi_pert.toPhys());
		o_vrt_t -= V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, U_vrt.toPhys());
		o_div_t -= V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, U_div.toPhys());
	}


	if (use_nr)
	{
		SphereData_Physical U_div_phys = U_div.toPhys();

		o_phi_pert_t -= SphereData_Spectral(U_phi_pert.toPhys()*U_div.toPhys());

		o_vrt_t -= o_vrt_t.toPhys()*U_div_phys;

		const SphereData_Physical U_vrt_phys = U_vrt.toPhys();;
		o_div_t += op.uv_to_vort(U_vrt_phys*U_u_phys, U_vrt_phys*U_v_phys);
		o_div_t += op.uv_to_div(U_div_phys*U_u_phys, U_div_phys*U_v_phys);
		o_div_t -= 0.5*op.laplace(U_u_phys*U_u_phys + U_v_phys*U_v_phys);
		o_div_t -= U_div_phys*U_div_phys;
	}
}



void SWE_Sphere_TS_ln_erk_split_vd::run_timestep_pert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,		///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_ln_erk_split_vd::euler_timestep_update_pert,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
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
		bool i_use_nr
)
{
	use_lg = i_use_lg;
	use_lc = i_use_lc;
	use_na = i_use_na;
	use_nr = i_use_nr;

	timestepping_order = i_order;
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

