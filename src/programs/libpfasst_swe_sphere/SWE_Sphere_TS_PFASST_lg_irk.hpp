/*
 * SWEImplicit_SPHRobert.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_PFASST_LG_IRK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_PFASST_LG_IRK_HPP_



#include <complex>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_PFASST_interface.hpp"



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWE_Sphere_TS_PFASST_lg_irk	: public SWE_Sphere_TS_PFASST_interface
{
	/// Simulation variables
	SimulationVariables &simVars;

	/// Operators for sphere
	SphereOperators_SphereData &op;

	/// SPH configuration
	const SphereData_Config *sphereDataConfig;

	/// alpha/beta (time step related component for implicit solver)
	double alpha;
	double beta;

	// Order of time stepping.
	int timestepping_order;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// Average geopotential
	double gh;

public:
	SWE_Sphere_TS_PFASST_lg_irk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
	);


	/**
	 * Setup the SWE REXI solver with SPH
	 */
public:
	void setup(
		int i_timestep_order,
		double i_timestep_size
	);


	/**
	 * Solve a REXI time step for the given initial conditions
	 */
public:
	void run_timestep_nonpert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp = -1
	);

	virtual ~SWE_Sphere_TS_PFASST_lg_irk();
};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
