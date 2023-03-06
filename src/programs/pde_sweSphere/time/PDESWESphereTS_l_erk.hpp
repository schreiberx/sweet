/*
 * PDESWESphereTS_l_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"



class PDESWESphereTS_l_erk final	: public PDESWESphereTS_BaseInterface
{

public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackDict.disc.timestepping_order;
		//timestepping_order2 = shackDict.disc.timestepping_order2;
		return i_timestepping_method == "l_erk";
	}

public:
	std::string string_id()
	{
		return "l_erk";
	}

private:
	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &op;

	int timestepping_order;

	// Sampler
	SphereTimestepping_ExplicitRK timestepping_rk;

public:
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
	PDESWESphereTS_l_erk(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup_auto();

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~PDESWESphereTS_l_erk();
};

#endif
