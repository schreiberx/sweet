/*
 * PDESWESphereTS_l_rexi_n_erk.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SWE_SPHERE_TS_L_REXI_N_ERK_HPP_
#define SWE_SPHERE_TS_L_REXI_N_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_erk_n_erk.hpp"
#include "PDESWESphereTS_l_exp.hpp"



class PDESWESphereTS_l_exp_n_erk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		);

	bool setup_main(
			sweet::SphereOperators *io_ops,
			sweet::ShackExpIntegration *i_shackExpIntegration,
			const std::string &i_exp_method,
			int i_order,	///< order of RK time stepping method
			int i_order2,	///< order of RK time stepping method of non-linear parts
			double i_timestep_size,
			bool i_use_f_sphere,
			int i_version_id,
			bool i_use_rexi_sphere_solver_preallocation
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method)
	{
		timestepping_method = i_timestepping_method;

		if (
			i_timestepping_method == "l_exp_n_erk" || i_timestepping_method == "l_exp_n_erk_ver0" ||
			i_timestepping_method == "l_exp_n_erk_ver1"
		)
			return true;

		return false;
	}

public:
	std::string getIDString()
	{
		std::string s = "l_exp_n_erk_ver";

		if (version_id == 0)
			s += "0";
		else if (version_id == 1)
			s += "1";
		else
			SWEETError("Version ID");

		return s;
	}


	double timestep_size;

	/*
	 * Linear time steppers
	 */
	PDESWESphereTS_l_exp timestepping_l_rexi;

	/*
	 * Non-linear time steppers
	 */
	PDESWESphereTS_l_erk_n_erk timestepping_l_erk_n_erk;

	sweet::TimesteppingExplicitRKSphereData timestepping_rk_nonlinear;

	int version_id;


public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_l_rexi.shackRegistration(io_shackDict);
		timestepping_l_erk_n_erk.shackRegistration(io_shackDict);
		return true;
	}


public:
	PDESWESphereTS_l_exp_n_erk();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~PDESWESphereTS_l_exp_n_erk();
};

#endif
