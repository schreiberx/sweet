/*
 * SWE_Sphere_TS_l_irk_n_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

# pragma once

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"



class SWE_Sphere_TS_ln_imex_sdc	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();

	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_l_irk timestepping_l_irk;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_l_erk_n_erk timestepping_l_erk_n_erk;

	/*
	 * To be specified ...
 	 */
	int version_id;

	int timestepping_order;
	int timestepping_order2;

private:
	/*
	 * SDC specific attributes
 	 */
	
	/*
	 * SDC specific methods
 	 */

	// Wrapper evaluating linear terms and storing them in separate variables
	void evalLinearTerms(
			SphereData_Spectral &phi_pert,	///< prognostic variables
			SphereData_Spectral &vort,	    ///< prognostic variables
			SphereData_Spectral &div,	    ///< prognostic variables
			SphereData_Spectral &phi_pert_L,	///< evaluation
			SphereData_Spectral &vort_L,	    ///< evaluation
			SphereData_Spectral &div_L,	        ///< evaluation
			double simulation_timestamp = -1
	);

	// Wrapper evaluating non-linear terms and storing them in separate variables
	void evalNonLinearTerms(
			SphereData_Spectral &phi_pert,	///< prognostic variables
			SphereData_Spectral &vort,	    ///< prognostic variables
			SphereData_Spectral &div,	    ///< prognostic variables
			SphereData_Spectral &phi_pert_NL,	///< evaluation
			SphereData_Spectral &vort_NL,	    ///< evaluation
			SphereData_Spectral &div_NL,	    ///< evaluation
			double simulation_timestamp = -1
	);

	/* Wrapper solving the implicit system built from the linear term :
	 u - dt*L(u) = rhs
	 LHS is always updated, independing on the dt value
	 WARNING : rhs variables are overwritten with u 
	*/ 
	void solveImplicit(
		SphereData_Spectral &rhs_phi,	///< rhs variables
		SphereData_Spectral &rhs_vrt,	///< rhs variables
		SphereData_Spectral &rhs_div,	///< rhs variables

		double dt
	);

public:
	SWE_Sphere_TS_ln_imex_sdc(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method for linear parts
			int i_order2,	///< order of RK time stepping method for non-linear parts
			int i_version_id
	);

	void setup_auto();

	void run_timestep(
			SphereData_Spectral &io_phi_pert,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_ln_imex_sdc();
};
