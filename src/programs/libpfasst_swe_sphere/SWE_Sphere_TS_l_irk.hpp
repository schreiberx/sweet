/*
 * SWEImplicit_SPHRobert.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_IRK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_IRK_HPP_


#include <complex>
#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWE_Sphere_TS_l_irk	: public SWE_Sphere_TS_interface
{
	/// Simulation variables
	SimulationVariables &simVars;

	/// Operators for sphere
	SphereOperators &op;

	/// SPH configuration
	const SphereDataConfig *sphereDataConfig;

	/// SPH configuration used for solver (maybe extended modes)
	const SphereDataConfig *sphereDataConfigSolver;
	SphereDataConfig sphereDataConfigSolverAddedModes;

	/// Solvers for alpha=Identity
	/// Template parameter is still complex-valued!!!
	/// This is because the spectral space is complex valued
	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolverPhi;
	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolverVel;

	bool use_extended_modes;

	/// alpha/beta (time step related component for implicit solver)
	double alpha;
	double beta;

	bool use_f_sphere;

	// Order of time stepping.
	int timestepping_order;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// f0
	double f0;

	/// Coriolis effect
	double two_coriolis;

	/// Average geopotential
	double gh;

	SphereDataPhysical mug;

public:
	SWE_Sphere_TS_l_irk(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
	);

private:
	void update_coefficients();


	/**
	 * Setup the SWE REXI solver with SPH
	 */
public:
	void setup(
			int i_timestep_order,
			double i_timestep_size,
			int i_use_extended_modes
	);


	/**
	 * Solve a REXI time step for the given initial conditions
	 */
public:
	void run_timestep(
			SphereDataSpectral &io_phi,		///< prognostic variables
			SphereDataSpectral &io_vort,	///< prognostic variables
			SphereDataSpectral &io_div,		///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Sphere_TS_l_irk();
};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
