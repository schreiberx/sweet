/*
 * PDESWESphereTS_lg_erk.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"



class PDESWESphereTS_lg_erk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup_main(
			sweet::SphereOperators *io_ops,
			int i_timestepping_order
		);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		return i_timestepping_method == "lg_erk";
	}

	std::string getIDString() override
	{
		return "lg_erk";
	}


private:
	// Sampler
	sweet::TimesteppingExplicitRKSphereData timestepping_rk;

public:
	void euler_timestep_update(
			const sweet::SphereData_Spectral &i_phi,
			const sweet::SphereData_Spectral &i_vort,
			const sweet::SphereData_Spectral &i_div,

			sweet::SphereData_Spectral &o_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_vort_t,	///< time updates
			sweet::SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	PDESWESphereTS_lg_erk();

	bool setup(
			int i_order	///< order of RK time stepping method
	);

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vort,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	virtual ~PDESWESphereTS_lg_erk();
};

#endif
