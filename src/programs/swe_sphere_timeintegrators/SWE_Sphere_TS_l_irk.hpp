/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_IRK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_IRK_HPP_


#include <complex>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>

#include "helpers/SWESphBandedMatrixPhysicalReal.hpp"
#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_erk.hpp"
#include "SWE_Sphere_TS_lg_erk.hpp"



/**
 * Implicit solver
 */
class SWE_Sphere_TS_l_irk	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();
	void setup_auto();

	std::string timestepping_method;

private:
	/// Simulation variables
	SimulationVariables &simVars;

	/// Operators for sphere
	SphereOperators_SphereData &ops;

	/// SPH configuration
	const SphereData_Config *sphereDataConfig;

	SWE_Sphere_TS_lg_erk *lg_erk = nullptr;
	SWE_Sphere_TS_l_erk *l_erk = nullptr;

	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolverDiv;

	// Order of time stepping.
	int timestepping_order;

	double crank_nicolson_damping_factor;

	/// timestep size
	double timestep_size;

	/// individual time step size
	double dt_explicit = -1;
	double dt_implicit = -1;

	/// earth radius
	double sphere_radius;

	bool use_f_sphere;

	bool no_coriolis;

	/// f0
	double f0;

	/// Coriolis effect
	double two_coriolis;

	SphereData_Physical mug;

public:
	SWE_Sphere_TS_l_irk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
	);

private:
	void update_coefficients(double i_timestep_size);


public:
	void setup(
			int i_timestep_order,
			double i_timestep_size,
			double i_crank_nicolson_damping_factor,
			bool i_no_coriolis
	);

public:
	void setup(
			int i_timestep_order,
			double i_timestep_size
	);

public:
	void free();

public:
	void run_timestep(
			SphereData_Spectral &io_phi,
			SphereData_Spectral &io_vrt,
			SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void solveImplicit(
		SphereData_Spectral &rhs_phi,
		SphereData_Spectral &rhs_vrt,
		SphereData_Spectral &rhs_div,

		double dt
	);


	virtual ~SWE_Sphere_TS_l_irk();
};


#endif
