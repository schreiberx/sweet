/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_IRK_LC_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_IRK_LC_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_lg_erk_lc_erk.hpp"
#include "PDESWESphereTS_lg_irk.hpp"



class PDESWESphereTS_lg_irk_lc_erk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup_main(
			sweet::SphereOperators *io_ops,
			int i_order,	///< order of RK time stepping method
			int i_version_id
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		return (
				i_timestepping_method == "lg_irk_lc_erk" ||
				i_timestepping_method == "lg_irk_lc_erk_ver0"	||
				i_timestepping_method == "lg_irk_lc_erk_ver1"
		);
	}

public:
	std::string getIDString() override
	{
		std::string s = "lg_irk_lc_erk_ver";

		if (version_id == 0)
			s += "0";
		else if (version_id == 1)
			s += "1";
		else
			SWEETError("Version ID");

		return s;
	}

	int version_id;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	PDESWESphereTS_lg_irk timestepping_lg_irk;

	/*
	 * Non-linear time steppers
	 */
	PDESWESphereTS_lg_erk_lc_erk timestepping_lg_erk_lc_erk;

	sweet::TimesteppingExplicitRKSphereData timestepping_rk_nonlinear;


public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_lg_irk.shackRegistration(io_shackDict);
		timestepping_lg_erk_lc_erk.shackRegistration(io_shackDict);
		return true;
	}



public:
	PDESWESphereTS_lg_irk_lc_erk();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vort,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	virtual ~PDESWESphereTS_lg_irk_lc_erk();
};

#endif
