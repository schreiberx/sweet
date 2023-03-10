/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_REXI_LC_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_REXI_LC_N_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_exp.hpp"
#include "PDESWESphereTS_lg_erk_lc_n_erk.hpp"



class PDESWESphereTS_lg_exp_lc_n_erk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		);

	bool setup_main(
			sweet::SphereOperators *io_ops,
			sweet::ShackExpIntegration *i_shackExpIntegration,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size,
			int i_version_id
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		if (
			i_timestepping_method == "lg_exp_lc_n_erk" || i_timestepping_method == "lg_exp_lc_n_erk_ver0" ||
			i_timestepping_method == "lg_exp_lc_n_erk_ver1"
		)
			return true;

		return false;
	}

	std::string getIDString()
	{
		std::string s = "lg_exp_lc_n_erk_ver";

		if (version_id == 0)
			s += "0";
		else if (version_id == 1)
			s += "1";
		else
			SWEETError("Version ID");

		return s;
	}

private:
	int version_id;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	PDESWESphereTS_l_exp timestepping_lg_rexi;

	/*
	 * Non-linear time steppers
	 */
	PDESWESphereTS_lg_erk_lc_n_erk timestepping_lg_erk_lc_n_erk;

	sweet::TimesteppingExplicitRKSphereData timestepping_rk_nonlinear;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_lg_rexi.shackRegistration(io_shackDict);
		timestepping_lg_erk_lc_n_erk.shackRegistration(io_shackDict);
		return true;
	}

public:
	PDESWESphereTS_lg_exp_lc_n_erk();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,		///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,		///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~PDESWESphereTS_lg_exp_lc_n_erk();
};

#endif
