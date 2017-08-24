/*
 * SWE_Sphere_TS_l_irk.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */


#include "SWE_Sphere_TS_l_irk.hpp"
#include <complex>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>


SWE_Sphere_TS_l_irk::SWE_Sphere_TS_l_irk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
)	:
	simVars(i_simVars),
	op(i_op),
	sphereDataConfig(op.sphereDataConfig)
{
}


/**
 * Setup the SWE REXI solver with SPH
 */
void SWE_Sphere_TS_l_irk::setup(
		int i_timestep_order,
		double i_timestep_size,
		int i_use_extended_modes
)
{
	if (i_timestep_order != 1)
		FatalError("Only 1st order IRK supported so far!");

	use_extended_modes = i_use_extended_modes;

	if (use_extended_modes == 0)
	{
		sphereDataConfigSolver = sphereDataConfig;
	}
	else
	{
		// Add modes only along latitude since these are the "problematic" modes
		sphereDataConfigSolverAddedModes.setupAdditionalModes(
				sphereDataConfig,
				use_extended_modes,	// TODO: Extend SPH wrapper to also support m != n to set this guy to 0
				use_extended_modes
		);

		sphereDataConfigSolver = &sphereDataConfigSolverAddedModes;
	}

	timestep_size = i_timestep_size;
	use_f_sphere = simVars.sim.f_sphere;

	if (use_f_sphere)
		f0 = simVars.sim.f0;
	else
		two_coriolis = 2.0*simVars.sim.coriolis_omega;

	alpha = -1.0/timestep_size;
	beta = -1.0/timestep_size;

	r = simVars.sim.earth_radius;
	inv_r = 1.0/r;

	gh = simVars.sim.gravitation*simVars.sim.h0;


	update_coefficients();
}


void SWE_Sphere_TS_l_irk::update_coefficients()
{

	if (!use_f_sphere)
	{
		sphSolverPhi.setup(sphereDataConfigSolver, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
		sphSolverPhi.solver_component_rexi_z2(	2.0*two_coriolis*two_coriolis*alpha*alpha, r);
		sphSolverPhi.solver_component_rexi_z3(	(two_coriolis*two_coriolis)*(two_coriolis*two_coriolis), r);
		sphSolverPhi.solver_component_rexi_z4robert(	-gh*alpha*two_coriolis, r);
		sphSolverPhi.solver_component_rexi_z5robert(	gh/alpha*two_coriolis*two_coriolis*two_coriolis, r);
		sphSolverPhi.solver_component_rexi_z6robert(	gh*2.0*two_coriolis*two_coriolis, r);
		sphSolverPhi.solver_component_rexi_z7(	-gh*alpha*alpha, r);
		sphSolverPhi.solver_component_rexi_z8(	-gh*two_coriolis*two_coriolis, r);

		mug.setup(sphereDataConfig);
		mug.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu;
			}
		);
	}
}


/**
 * Solve a REXI time step for the given initial conditions
 */

void SWE_Sphere_TS_l_irk::run_timestep(
		SphereData &io_phi,		///< prognostic variables
		SphereData &io_vort,	///< prognostic variables
		SphereData &io_div,		///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	SphereData phi0 = io_phi;
	SphereData vort0 = io_vort;
	SphereData div0 = io_div;

	SphereData phi(sphereDataConfig);
	SphereData vort(sphereDataConfig);
	SphereData div(sphereDataConfig);

	if (use_f_sphere)
	{
		SphereData rhs = gh*(div0 - f0/alpha*vort0) + (alpha+f0*f0/alpha)*phi0;
		phi = rhs.spectral_solve_helmholtz(alpha*alpha + f0*f0, -gh, r);

		div = -1.0/gh*(phi0 - alpha*phi);
		vort = (1.0/alpha)*(vort0 + f0*(div));
	}
	else
	{
		SphereData rhs(sphereDataConfig);

		SphereDataPhysical u0g(sphereDataConfig);
		SphereDataPhysical v0g(sphereDataConfig);
		op.robert_vortdiv_to_uv(vort0, div0, u0g, v0g);

		SphereDataPhysical phi0g = phi0.getSphereDataPhysical();

		SphereDataPhysical Fc_k =
				two_coriolis*inv_r*(
						-(-two_coriolis*two_coriolis*mug*mug + alpha*alpha)*u0g
						+ 2.0*alpha*two_coriolis*mug*v0g
				);

		SphereDataPhysical foo =
				(gh*(div0.getSphereDataPhysical() - (1.0/alpha)*two_coriolis*mug*vort0.getSphereDataPhysical())) +
				(alpha*phi0g + (1.0/alpha)*two_coriolis*two_coriolis*mug*mug*phi0g);

		SphereDataPhysical rhsg =
				alpha*alpha*foo +
				two_coriolis*two_coriolis*mug*mug*foo
				- (gh/alpha)*Fc_k;

		rhs = rhsg;

		phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);


		SphereDataPhysical u0(sphereDataConfig);
		SphereDataPhysical v0(sphereDataConfig);

		op.robert_vortdiv_to_uv(vort0, div0, u0, v0);

		SphereDataPhysical gradu(sphereDataConfig);
		SphereDataPhysical gradv(sphereDataConfig);
		op.robert_grad_to_vec(phi, gradu, gradv, r);

		SphereDataPhysical a = u0 + gradu;
		SphereDataPhysical b = v0 + gradv;

		SphereDataPhysical k = (two_coriolis*two_coriolis*mug*mug+alpha*alpha);
		SphereDataPhysical u = (alpha*a - two_coriolis*mug*(b))/k;
		SphereDataPhysical v = (two_coriolis*mug*(a) + alpha*b)/k;

		op.robert_uv_to_vortdiv(u, v, vort, div);
	}

	io_phi = phi * beta;
	io_vort = vort * beta;
	io_div = div * beta;
}


SWE_Sphere_TS_l_irk::~SWE_Sphere_TS_l_irk()
{
}
