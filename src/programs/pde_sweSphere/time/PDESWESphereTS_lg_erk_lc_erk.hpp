/*
 * PDESWESphereTS_lg_erk_lc_erk.hpp
 *
 *  Created on: 11 November 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_LC_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_LC_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"



class PDESWESphereTS_lg_erk_lc_erk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		);

	bool setup_main(
			sweet::SphereOperators *io_ops,
			int i_order	///< order of RK time stepping method
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method)
	{
		timestepping_method = i_timestepping_method;
		return i_timestepping_method == "lg_erk_lc_erk";
	}

public:
	std::string getIDString()
	{
		return "lg_erk_lc_erk";
	}

	sweet::TimesteppingExplicitRKSphereData timestepping_rk_linear;
	sweet::TimesteppingExplicitRKSphereData timestepping_rk_nonlinear;


public:
	void euler_timestep_update_lg(
			const sweet::SphereData_Spectral &i_h,	///< prognostic variables
			const sweet::SphereData_Spectral &i_u,	///< prognostic variables
			const sweet::SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_lc(
			const sweet::SphereData_Spectral &i_h,	///< prognostic variables
			const sweet::SphereData_Spectral &i_u,	///< prognostic variables
			const sweet::SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_lc(
			sweet::SphereData_Spectral &io_phi,		///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,		///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


private:
	void euler_timestep_update(
			const sweet::SphereData_Spectral &i_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_vort,	///< prognostic variables
			const sweet::SphereData_Spectral &i_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_vort_t,	///< time updates
			sweet::SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	PDESWESphereTS_lg_erk_lc_erk();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~PDESWESphereTS_lg_erk_lc_erk();
};

#endif
