/*
 * SWE_Sphere_TS_l_irk.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "SWE_Sphere_TS_l_irk.hpp"
#include <complex>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>



void SWE_Sphere_TS_l_irk::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 2)
	{
		// First run a half explicit time step for 2nd order accurate IRK (Crank-Nicolson)
		l_erk->run_timestep_pert(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
	}


	double gh0 = simVars.sim.gravitation*simVars.sim.h0;
	io_phi_pert += gh0;
	run_timestep_nonpert_backward_euler_private(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
	io_phi_pert -= gh0;
}



/**
 * Solve an implicit time step for the given initial conditions
 */
void SWE_Sphere_TS_l_irk::run_timestep_nonpert_backward_euler_private(
		SphereData_Spectral &io_phi,
		SphereData_Spectral &io_vrt,
		SphereData_Spectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (std::abs(timestep_size - i_fixed_dt)/std::max(timestep_size, i_fixed_dt) > 1e-10)
	{
	        std::cout << "Warning: Reducing time step size from " << i_fixed_dt << " to " << timestep_size << std::endl;

	        timestep_size = i_fixed_dt;
	        update_coefficients();
	}

	SphereData_Spectral phi0 = io_phi;
	SphereData_Spectral vort0 = io_vrt;
	SphereData_Spectral div0 = io_div;

	SphereData_Spectral phi(sphereDataConfig);
	SphereData_Spectral vort(sphereDataConfig);
	SphereData_Spectral div(sphereDataConfig);


	double gh0 = simVars.sim.gravitation*simVars.sim.h0;
	if (use_f_sphere)
	{
		SphereData_Spectral rhs = gh0*(div0 - f0/alpha*vort0) + (alpha+f0*f0/alpha)*phi0;
		phi = rhs.spectral_solve_helmholtz(alpha*alpha + f0*f0, -gh0, r);

		div = -1.0/gh0*(phi0 - alpha*phi);
		vort = (1.0/alpha)*(vort0 + f0*(div));
	}
	else
	{
		if (!simVars.misc.sphere_use_robert_functions)
			SWEETError("Using no Robert formulation is not yet supported in this time integrator");

		SphereData_Spectral rhs(sphereDataConfig);

		SphereData_Physical u0g(sphereDataConfig);
		SphereData_Physical v0g(sphereDataConfig);
		ops.robert_vortdiv_to_uv(vort0, div0, u0g, v0g);

		SphereData_Physical phi0g = phi0.getSphereDataPhysical();

		SphereData_Physical Fc_k =
				two_coriolis*inv_r*(
						-(-two_coriolis*two_coriolis*mug*mug + alpha*alpha)*u0g
						+ 2.0*alpha*two_coriolis*mug*v0g
				);

		SphereData_Physical foo =
				(gh0*(div0.getSphereDataPhysical() - (1.0/alpha)*two_coriolis*mug*vort0.getSphereDataPhysical())) +
				(alpha*phi0g + (1.0/alpha)*two_coriolis*two_coriolis*mug*mug*phi0g);

		SphereData_Physical rhsg =
				alpha*alpha*foo +
				two_coriolis*two_coriolis*mug*mug*foo
				- (gh0/alpha)*Fc_k;

		rhs = rhsg;

		phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

		SphereData_Physical u0(sphereDataConfig);
		SphereData_Physical v0(sphereDataConfig);

		ops.robert_vortdiv_to_uv(vort0, div0, u0, v0);

		SphereData_Physical gradu(sphereDataConfig);
		SphereData_Physical gradv(sphereDataConfig);
		ops.robert_grad_to_vec(phi, gradu, gradv, r);

		SphereData_Physical a = u0 + gradu;
		SphereData_Physical b = v0 + gradv;

		SphereData_Physical k = (two_coriolis*two_coriolis*mug*mug+alpha*alpha);
		SphereData_Physical u = (alpha*a - two_coriolis*mug*(b))/k;
		SphereData_Physical v = (two_coriolis*mug*(a) + alpha*b)/k;

		ops.robert_uv_to_vortdiv(u, v, vort, div);
	}

	io_phi = phi * beta;
	io_vrt = vort * beta;
	io_div = div * beta;
}



void SWE_Sphere_TS_l_irk::update_coefficients()
{
    alpha = -1.0/timestep_size;
    beta = -1.0/timestep_size;

	alpha /= crank_nicolson_damping_factor;
	beta /= crank_nicolson_damping_factor;

	if (!use_f_sphere)
	{
		double gh0 = simVars.sim.gravitation*simVars.sim.h0;
		sphSolverPhi.setup(sphereDataConfigSolver, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
		sphSolverPhi.solver_component_rexi_z2(	2.0*two_coriolis*two_coriolis*alpha*alpha, r);
		sphSolverPhi.solver_component_rexi_z3(	(two_coriolis*two_coriolis)*(two_coriolis*two_coriolis), r);
		sphSolverPhi.solver_component_rexi_z4robert(	-gh0*alpha*two_coriolis, r);
		sphSolverPhi.solver_component_rexi_z5robert(	gh0/alpha*two_coriolis*two_coriolis*two_coriolis, r);
		sphSolverPhi.solver_component_rexi_z6robert(	gh0*2.0*two_coriolis*two_coriolis, r);
		sphSolverPhi.solver_component_rexi_z7(	-gh0*alpha*alpha, r);
		sphSolverPhi.solver_component_rexi_z8(	-gh0*two_coriolis*two_coriolis, r);

		mug.setup_if_required(sphereDataConfig);
		mug.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu;
			}
		);
	}
}



SWE_Sphere_TS_l_irk::SWE_Sphere_TS_l_irk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
	simVars(i_simVars),
	ops(i_op),
	sphereDataConfig(ops.sphereDataConfig)
{
}


/**
 * Setup the SWE implicit solver with SPH
 */
void SWE_Sphere_TS_l_irk::setup(
		int i_timestep_order,
		double i_timestep_size,
		int i_use_extended_modes,
		double i_crank_nicolson_damping_factor
)
{
	timestepping_order = i_timestep_order;
	timestep_size = i_timestep_size;
	use_f_sphere = simVars.sim.sphere_use_fsphere;

	if (timestepping_order == 1)
	{
		// set this to 1 to ignore it
		crank_nicolson_damping_factor = 1.0;
	}
	else if (timestepping_order == 2)
	{
		// set this to 0.5 to take 1/2 TS with forward Euler and another 1/2 one with backward Euler
		crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;

		l_erk = new SWE_Sphere_TS_l_erk(simVars, ops);
		l_erk->setup(1);
	}
	else
	{
		SWEETError("Only 1st and 2nd order IRK supported so far with this implementation! Use l_cn if you want to have 2nd order Crank-Nicolson!");
	}



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
				use_extended_modes,
				simVars.misc.reuse_spectral_transformation_plans
		);

		sphereDataConfigSolver = &sphereDataConfigSolverAddedModes;
	}

	if (use_f_sphere)
	{
		f0 = simVars.sim.sphere_fsphere_f0;
		two_coriolis = 0;
	}
	else
	{
		two_coriolis = 2.0*simVars.sim.sphere_rotating_coriolis_omega;
		f0 = 0;
	}

	r = simVars.sim.sphere_radius;
	inv_r = 1.0/r;

	update_coefficients();
}


SWE_Sphere_TS_l_irk::~SWE_Sphere_TS_l_irk()
{
	delete l_erk;
	l_erk = nullptr;
}
