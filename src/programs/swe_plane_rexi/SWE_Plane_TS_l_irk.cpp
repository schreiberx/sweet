/*
 * SWE_Plane_TS_l_irk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Plane_TS_l_irk.hpp"

#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>




/**
 * Solve SWE with implicit time stepping
 *
 * U_t = L U(0)
 *
 * (U(tau) - U(0)) / tau = L U(tau)
 *
 * <=> U(tau) - U(0) = L U(tau) tau
 *
 * <=> U(tau) - L tau U(tau) = U(0)
 *
 * <=> (1 - L tau) U(tau) = U(0)
 *
 * <=> (1/tau - L) U(tau) = U(0)/tau
 */
void SWE_Plane_TS_l_irk::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt <= 0)
		FatalError("OSWE_Plane_TS_l_irk: nly constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time-i_simulation_timestamp;

	PlaneData &eta0 = io_h;
	PlaneData &u0 = io_u;
	PlaneData &v0 = io_v;

	double alpha = 1.0/i_fixed_dt;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + simVars.sim.f0*simVars.sim.f0;

	double eta_bar = simVars.sim.h0;
	double g = simVars.sim.gravitation;

	PlaneData rhs =
			(kappa/alpha) * eta0
			- eta_bar*(op.diff_c_x(u0) + op.diff_c_y(v0))
			- (simVars.sim.f0*eta_bar/alpha) * (op.diff_c_x(v0) - op.diff_c_y(u0))
		;

	PlaneData lhs = (-g*eta_bar*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(kappa);
	io_h = rhs.spectral_div_element_wise(lhs);

	PlaneData uh = u0 - g*op.diff_c_x(io_h);
	PlaneData vh = v0 - g*op.diff_c_y(io_h);

	io_u = alpha/kappa * uh     + simVars.sim.f0/kappa * vh;
	io_v =    -simVars.sim.f0/kappa * uh + alpha/kappa * vh;

	o_dt = i_fixed_dt;
}



/*
 * Setup
 */
void SWE_Plane_TS_l_irk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	if (timestepping_order != 1)
		FatalError("SWE_Plane_TS_l_irk: Only 1st order IRK is supported");
}



SWE_Plane_TS_l_irk::SWE_Plane_TS_l_irk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(1);
}



SWE_Plane_TS_l_irk::~SWE_Plane_TS_l_irk()
{
}

