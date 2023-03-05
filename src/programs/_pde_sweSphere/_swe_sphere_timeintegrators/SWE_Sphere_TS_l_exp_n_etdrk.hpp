/*
 * SWE_Sphere_TS_l_phi0_n_edt.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_Sphere_TS_l_rexi_n_etdrk_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_Sphere_TS_l_rexi_n_etdrk_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"


class SWE_Sphere_TS_l_exp_n_etdrk	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();
	void setup_auto();

private:
	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &op;

	SWE_Sphere_TS_l_erk_n_erk ts_l_erk_n_erk;

	SWE_Sphere_TS_l_exp ts_phi0_exp;
	SWE_Sphere_TS_l_exp ts_phi1_exp;
	SWE_Sphere_TS_l_exp ts_phi2_exp;

	SWE_Sphere_TS_l_exp ts_ups0_exp;
	SWE_Sphere_TS_l_exp ts_ups1_exp;
	SWE_Sphere_TS_l_exp ts_ups2_exp;
	SWE_Sphere_TS_l_exp ts_ups3_exp;

	int timestepping_order;
	int timestepping_order2;

private:
	void euler_timestep_update_linear(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);



private:
	void euler_timestep_update_nonlinear(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);


public:
	SWE_Sphere_TS_l_exp_n_etdrk(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			EXP_sweet::ShackDictionary &i_rexi,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size
	);

	void run_timestep(
			sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_l_exp_n_etdrk();
};

#endif
