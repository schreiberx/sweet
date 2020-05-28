/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LG_IRK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LG_IRK_HPP_



#include <complex>
#include "helpers/SWESphBandedMatrixPhysicalReal.hpp"
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_lg_erk.hpp"


class SWE_Sphere_TS_lg_irk	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method);
	std::string string_id();
	void setup_auto();


private:
	/// Simulation variables
	SimulationVariables &simVars;

	/// Operators for sphere
	SphereOperators_SphereData &op;

	/// SPH configuration
	const SphereData_Config *sphereDataConfig;

	SWE_Sphere_TS_lg_erk *lg_erk = nullptr;

	/// alpha/beta (time step related component for implicit solver)
	double alpha;
	double beta;

	/// Crank-Nicolson damping factor
	double crank_nicolson_damping_factor = 0.5;

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
	SWE_Sphere_TS_lg_irk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
	);


public:
	void update_coefficients();


public:
	void setup(
		int i_timestep_order,
		double i_timestep_size
	);


public:
	void setup(
		int i_timestep_order,
		double i_timestep_size,
		double i_crank_nicolson_damping_factor
	);


public:
	void run_timestep(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp = -1
	);

	virtual ~SWE_Sphere_TS_lg_irk();
};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
