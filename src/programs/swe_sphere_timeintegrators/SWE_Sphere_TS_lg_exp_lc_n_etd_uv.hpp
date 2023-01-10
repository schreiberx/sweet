/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_LC_N_ETD_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_LC_N_ETD_UV_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"
#include "SWE_Sphere_TS_ln_erk_split_uv.hpp"


class SWE_Sphere_TS_lg_exp_lc_n_etd_uv	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();
	void setup_auto();
	void print_help();

private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &ops;

	SWE_Sphere_TS_ln_erk_split_uv ts_ln_erk_split_uv;

	SphereData_Spectral NU_phi_prev, NU_vrt_prev, NU_div_prev;
	SphereData_Spectral NU_phi_prev_2, NU_vrt_prev_2, NU_div_prev_2;

	bool with_na;
	bool with_nr;

	SWE_Sphere_TS_l_exp ts_phi0_exp;
	SWE_Sphere_TS_l_exp ts_phi1_exp;
	SWE_Sphere_TS_l_exp ts_phi2_exp;
	SWE_Sphere_TS_l_exp ts_phi3_exp;

	int timestepping_order;

private:
	void euler_timestep_update_linear(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			SphereData_Spectral &o_h_t,	///< time updates
			SphereData_Spectral &o_u_t,	///< time updates
			SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);



private:
	void euler_timestep_update_nonlinear(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			SphereData_Spectral &o_h_t,	///< time updates
			SphereData_Spectral &o_u_t,	///< time updates
			SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);


public:
	SWE_Sphere_TS_lg_exp_lc_n_etd_uv(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			EXP_SimulationVariables &i_rexi,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vrt,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_lg_exp_lc_n_etd_uv();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
