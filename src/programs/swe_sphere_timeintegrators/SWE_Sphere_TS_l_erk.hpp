/*
 * SWE_Sphere_TS_l_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



class SWE_Sphere_TS_l_erk final	: public SWE_Sphere_TS_interface
{

public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = simVars.disc.timestepping_order;
		//timestepping_order2 = simVars.disc.timestepping_order2;
		return i_timestepping_method == "l_erk";
	}

public:
	std::string string_id()
	{
		return "l_erk";
	}

private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int timestepping_order;

	// Sampler
	SphereTimestepping_ExplicitRK timestepping_rk;

public:
	void euler_timestep_update(
			const SphereData_Spectral &i_phi,	///< prognostic variables
			const SphereData_Spectral &i_vort,	///< prognostic variables
			const SphereData_Spectral &i_div,	///< prognostic variables

			SphereData_Spectral &o_phi_t,	///< time updates
			SphereData_Spectral &o_vort_t,	///< time updates
			SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);


public:
	SWE_Sphere_TS_l_erk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup_auto();

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereData_Spectral &io_phi_pert,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_l_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
