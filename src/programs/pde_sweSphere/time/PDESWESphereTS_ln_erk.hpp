/*
 * PDESWESphereTS_ln_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"


class PDESWESphereTS_ln_erk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		);

	bool setup(
			sweet::SphereOperators *io_ops,
			int i_order	///< order of RK time stepping method
	);


public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		if (i_timestepping_method == "ln_erk")
			return true;

		return false;
	}

	std::string getIDString()
	{
		return "ln_erk";
	}


private:
	// Sampler
	sweet::TimesteppingExplicitRKSphereData timestepping_rk;


public:
	void euler_timestep_update_pert(
			const sweet::SphereData_Spectral &i_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_vort,	///< prognostic variables
			const sweet::SphereData_Spectral &i_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_vort_t,	///< time updates
			sweet::SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	PDESWESphereTS_ln_erk();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~PDESWESphereTS_ln_erk();
};

#endif
