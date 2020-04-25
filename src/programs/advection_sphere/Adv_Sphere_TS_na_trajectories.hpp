/*
 * Adv_Sphere_TS_na_trajectories.hpp
 *
 *  Created on: 17 April 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_TRAJECTORIES_HPP_
#define SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_TRAJECTORIES_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include <sweet/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>
#include <benchmarks_sphere/SWESphereBenchmarksCombined.hpp>

#include "Adv_Sphere_TS_interface.hpp"


/*
 * This time integration method is based on given benchmark solutions
 *
 * See e.g. R. Nair, P. Lauritzen
 * "A class of deformational flow test cases for linear transport problems on the sphere"
 */
class Adv_Sphere_TS_na_trajectories	: public Adv_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int timestepping_order;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical &sampler2D;

	SphereData_Spectral U_phi_prev, U_vrt_prev, U_div_prev;


public:
	Adv_Sphere_TS_na_trajectories(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt,				///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp,

			// for varying velocity fields
			const SWESphereBenchmarksCombined *i_sphereBenchmarks,
			SphereData_Physical &io_U_phi_phys
	);



	virtual ~Adv_Sphere_TS_na_trajectories();
};

#endif
