/*
 * SWE_Sphere_TS_PFASST_l_irk_n_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_PFASST_L_IRK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_PFASST_L_IRK_N_ERK_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_PFASST_interface.hpp"
#include "SWE_Sphere_TS_PFASST_l_irk.hpp"
#include "SWE_Sphere_TS_PFASST_l_cn.hpp"
#include "SWE_Sphere_TS_PFASST_l_erk_n_erk.hpp"



class SWE_Sphere_TS_PFASST_l_irk_n_erk	: public SWE_Sphere_TS_PFASST_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_PFASST_l_irk timestepping_l_irk;
	SWE_Sphere_TS_PFASST_l_cn timestepping_l_cn;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_PFASST_l_erk_n_erk timestepping_l_erk_n_erk;

	SphereTimestepping_ExplicitRK timestepping_rk_nonlinear;

	int version_id;

	int timestepping_order;
	int timestepping_order2;


public:
	SWE_Sphere_TS_PFASST_l_irk_n_erk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method for linear parts
			int i_order2,	///< order of RK time stepping method for non-linear parts
			int i_version_id
	);

	void run_timestep_nonpert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_PFASST_l_irk_n_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_PFASST_LN_ERK_HPP_ */
