/*
 * SWE_Sphere_TS_l_cn.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_HPP_


#include <complex>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWE_Sphere_TS_l_cn_DEPRECATED	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method)
	{
		if (i_timestepping_method == "l_cn_DEPRECATED" || i_timestepping_method == "l_irk_DEPRECATED")
			return true;

		return false;
	}

	std::string string_id()
	{
		return "l_irk_DEPRECATED";
	}

	void setup_auto()
	{
		setup(
				simVars.disc.timestepping_crank_nicolson_filter,
				simVars.timecontrol.current_timestep_size,
				simVars.rexi.use_sphere_extended_modes
			);
	}


private:
	/// Simulation variables
	SimulationVariables &simVars;

	/// Operators for sphere
	SphereOperators_SphereData &op;

	/// SPH configuration
	const SphereData_Config *sphereDataConfig;

	/// SPH configuration used for solver (maybe extended modes)
	const SphereData_Config *sphereDataConfigSolver;
	SphereData_Config sphereDataConfigSolverAddedModes;

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

	/// coriolis
	double two_coriolis;

	/// f0
	double f0;

	/// Average geopotential
	double gh;

	SphereData_Physical mug;

public:
	SWE_Sphere_TS_l_cn_DEPRECATED(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
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
	void run_timestep_pert(
			SphereData_Spectral &io_phi,		///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,		///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
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

	virtual ~SWE_Sphere_TS_l_cn_DEPRECATED();
};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
