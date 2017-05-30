/*
 * SWE_Sphere_TS_l_lf.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_LF_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_LF_HPP_

#include <limits>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataTimesteppingExplicitLeapfrog.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



class SWE_Sphere_TS_l_lf	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators &op;

	int timestepping_order;

	// Sampler
	SphereDataTimesteppingExplicitLeapfrog timestepping_lf;

	// Coriolis effect
	SphereDataPhysical fg;

private:
	void euler_timestep_update(
			const SphereData &i_phi,	///< prognostic variables
			const SphereData &i_vort,	///< prognostic variables
			const SphereData &i_div,	///< prognostic variables

			SphereData &o_phi_t,	///< time updates
			SphereData &o_vort_t,	///< time updates
			SphereData &o_div_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);

public:
	SWE_Sphere_TS_l_lf(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereData &io_phi,	///< prognostic variables
			SphereData &io_vort,	///< prognostic variables
			SphereData &io_div,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	);



	virtual ~SWE_Sphere_TS_l_lf();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
