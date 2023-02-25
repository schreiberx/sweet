/*
 * Adv_Plane_TS_na_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_ERK_HPP_
#define SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/SimulationVariables.hpp>

#include "Adv_Plane_TS_interface.hpp"



class Adv_Plane_TS_na_erk	: public Adv_Plane_TS_interface
{
	SimulationVariables &simVars;
	sweet::PlaneOperators &op;

	int timestepping_order;

	// Sampler
	PlaneDataTimesteppingExplicitRK timestepping_rk;


private:
	void euler_timestep_update(
			const sweet::PlaneData_Spectral &i_phi,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_vort,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_div,	///< prognostic variables

			sweet::PlaneData_Spectral &o_phi_t,	///< time updates
			sweet::PlaneData_Spectral &o_vort_t,	///< time updates
			sweet::PlaneData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	Adv_Plane_TS_na_erk(
			SimulationVariables &i_simVars,
			sweet::PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			sweet::PlaneData_Spectral &io_phi,	///< prognostic variables
			sweet::PlaneData_Spectral &io_vort,	///< prognostic variables
			sweet::PlaneData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~Adv_Plane_TS_na_erk();
};

#endif /* SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_LN_ERK_HPP_ */
