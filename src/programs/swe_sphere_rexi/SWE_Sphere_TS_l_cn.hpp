/*
 * SWE_Sphere_TS_l_cn.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_HPP_


#include <complex>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWE_Sphere_TS_l_cn	: public SWE_Sphere_TS_interface
{
	/// Simulation variables
	SimulationVariables &simVars;

	/// Operators for sphere
	SphereOperators &op;

	/// SPH configuration
	SphereDataConfig *sphereDataConfig;

	/// SPH configuration used for solver (maybe extended modes)
	SphereDataConfig *sphereDataConfigSolver;
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

	/// Crank-Nicolson damping factor
	double crank_nicolson_damping_factor = 0.5;

	bool use_f_sphere;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// set to f0 for f-sphere and 2*coriolis_omega for non-f-sphere
	double coriolis;

	/// Average geopotential
	double gh;

	SphereDataPhysical fg;
	SphereDataPhysical mug;

public:
	SWE_Sphere_TS_l_cn(
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
			double i_crank_nicolson_damping_factor,
			double i_timestep_size,
			int i_use_extended_modes
	);


	/**
	 * Solve a REXI time step for the given initial conditions
	 */
public:
	void run_timestep(
			SphereData &io_phi,		///< prognostic variables
			SphereData &io_vort,	///< prognostic variables
			SphereData &io_div,		///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	);

	virtual ~SWE_Sphere_TS_l_cn();
};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
