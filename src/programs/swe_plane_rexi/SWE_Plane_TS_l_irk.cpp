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


	PlaneDataComplex eta(io_h.planeDataConfig);

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
#warning "WARNING: Not doing this in spectral space leads to a loss of the highest mode"
	/*
	 * WARNING: This leads to a loss of precision due to the highest mode which cannot be
	 * tracked due to the Nyquist theorem
	 */
	PlaneDataComplex eta0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_h);
	PlaneDataComplex u0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_u);
	PlaneDataComplex v0 = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_v);
#else
	PlaneDataComplex eta0 = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_h);
	PlaneDataComplex u0 = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_u);
	PlaneDataComplex v0 = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_v);
#endif

	double alpha = 1.0/i_fixed_dt;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + simVars.sim.f0*simVars.sim.f0;

	double eta_bar = simVars.sim.h0;
	double g = simVars.sim.gravitation;

	PlaneDataComplex rhs =
			(kappa/alpha) * eta0
			- eta_bar*(opComplex.diff_c_x(u0) + opComplex.diff_c_y(v0))
			- (simVars.sim.f0*eta_bar/alpha) * (opComplex.diff_c_x(v0) - opComplex.diff_c_y(u0))
		;

	helmholtz_spectral_solver_spec(kappa, g*eta_bar, rhs, eta, 0);

	PlaneDataComplex uh = u0 - g*opComplex.diff_c_x(eta);
	PlaneDataComplex vh = v0 - g*opComplex.diff_c_y(eta);

	PlaneDataComplex u1 = alpha/kappa * uh     + simVars.sim.f0/kappa * vh;
	PlaneDataComplex v1 =    -simVars.sim.f0/kappa * uh + alpha/kappa * vh;

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
#warning "WARNING: Not doing this in spectral space leads to a loss of the highest mode"
	/*
	 * WARNING: This leads to a loss of precision due to the highest mode which cannot be
	 * tracked due to the Nyquist theorem
	 */
	io_h = Convert_PlaneDataComplex_To_PlaneData::physical_convert(eta);
	io_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(u1);
	io_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(v1);
#else
	io_h = Convert_PlaneDataComplex_To_PlaneData::spectral_convert(eta);
	io_u = Convert_PlaneDataComplex_To_PlaneData::spectral_convert(u1);
	io_v = Convert_PlaneDataComplex_To_PlaneData::spectral_convert(v1);
#endif

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
		op(i_op),
		opComplex(op.planeDataConfig, simVars.sim.domain_size)
{
	setup(1);
}



SWE_Plane_TS_l_irk::~SWE_Plane_TS_l_irk()
{
}

