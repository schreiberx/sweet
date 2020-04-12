/*
 * SWE_Sphere_TS_l_lf.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_LF_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_LF_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitLeapfrog.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



class SWE_Sphere_TS_l_lf	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	// Sampler
	SphereTimestepping_ExplicitLeapfrog timestepping_lf;

	int timestepping_order;
	double robert_asselin_filter;

private:
	void euler_timestep_update(
			const SphereData_Spectral &i_phi,	///< prognostic variables
			const SphereData_Spectral &i_vort,	///< prognostic variables
			const SphereData_Spectral &i_div,	///< prognostic variables

			SphereData_Spectral &o_phi_t,	///< time updates
			SphereData_Spectral &o_vort_t,	///< time updates
			SphereData_Spectral &o_div_t,	///< time updates

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);


public:
	SWE_Sphere_TS_l_lf(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);


	void setup(
			int i_order,	///< order of RK time stepping method
			double i_robert_asselin_filter
	);


	void run_timestep_pert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);


	void run_timestep_nonpert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_l_lf();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
