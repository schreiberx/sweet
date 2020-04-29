/*
 * SWE_Sphere_TS_split_lg_lc_na_nr_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_ERK_SPLIT_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_ERK_SPLIT_UV_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "../../programs/swe_sphere/SWE_Sphere_TS_interface.hpp"


class SWE_Sphere_TS_ln_erk_split_uv	: public SWE_Sphere_TS_interface
{
public:
	static bool implements_timestepping_method(const std::string &i_timestepping_method)
	{
		if (i_timestepping_method == "ln_erk_split_uv")
			return true;

		return false;
	}

	std::string string_id()
	{
		return "ln_erk_split_uv";
	}

	void setup_auto()
	{
		setup(simVars.disc.timestepping_order, true, true, true, true);
	}

private:
	SimulationVariables &simVars;

	SphereOperators_SphereData &op;

	int timestepping_order;

	bool use_lg = false;
	bool use_lc = false;
	bool use_na = false;
	bool use_nr = false;


	// Sampler
	SphereTimestepping_ExplicitRK timestepping_rk;

private:
	SphereData_Spectral V_dot_grad_scalar(
			const SphereData_Physical &i_u_phys,		///< u velocity
			const SphereData_Physical &i_v_phys,		///< v velocity
			const SphereData_Physical &i_div_phys,		///< divergence in physical space to avoid transformation
			const SphereData_Physical &i_scalar_phys	///< scalar
	);

public:
	void euler_timestep_update_pert_lg(
			const SphereData_Spectral &i_U_phi,	///< prognostic variables
			const SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const SphereData_Spectral &i_U_div,	///< prognostic variables

			SphereData_Spectral &o_U_phi_t,	///< time updates
			SphereData_Spectral &o_U_vrt_t,	///< time updates
			SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_update_pert_lc(
			const SphereData_Spectral &i_U_phi,	///< prognostic variables
			const SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const SphereData_Spectral &i_U_div,	///< prognostic variables

			SphereData_Spectral &o_U_phi_t,	///< time updates
			SphereData_Spectral &o_U_vrt_t,	///< time updates
			SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_update_pert_na(
			const SphereData_Spectral &i_U_phi,	///< prognostic variables
			const SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const SphereData_Spectral &i_U_div,	///< prognostic variables

			SphereData_Spectral &o_U_phi_t,	///< time updates
			SphereData_Spectral &o_U_vrt_t,	///< time updates
			SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_update_pert_nr(
			const SphereData_Spectral &i_U_phi,	///< prognostic variables
			const SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const SphereData_Spectral &i_U_div,	///< prognostic variables

			SphereData_Spectral &o_U_phi_t,	///< time updates
			SphereData_Spectral &o_U_vrt_t,	///< time updates
			SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);


private:
	void euler_timestep_update_pert(
			const SphereData_Spectral &i_U_phi,	///< prognostic variables
			const SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const SphereData_Spectral &i_U_div,	///< prognostic variables

			SphereData_Spectral &o_U_phi_t,	///< time updates
			SphereData_Spectral &o_U_vrt_t,	///< time updates
			SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	SWE_Sphere_TS_ln_erk_split_uv(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			bool i_lg,
			bool i_lc,
			bool i_na,
			bool i_nr
	);

	void run_timestep_pert(
			SphereData_Spectral &io_U_phi,	///< prognostic variables
			SphereData_Spectral &io_U_vrt,	///< prognostic variables
			SphereData_Spectral &io_U_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);




	virtual ~SWE_Sphere_TS_ln_erk_split_uv();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
