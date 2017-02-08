/*
 * SWEImplicit_SPHRobert.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SWEIMPLICIT_SPHROBERT_HPP_
#define SRC_SWEIMPLICIT_SPHROBERT_HPP_

#include <complex>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>

#include "SWESphBandedMatrixPhysicalReal.hpp"



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWEImplicit_SPHRobert
{
	/// Operators for sphere
	SphereOperators &op;

	/// SPH configuration
	SphereDataConfig *sphereDataConfig;

	/// SPH configuration
	SphereDataConfig *sphereDataConfigSolver;

	/// Solvers for alpha=Identity
	/// Template parameter is still complex-valued!!!
	/// This is because the spectral space is complex valued
	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolverPhi;
	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolverVel;

	/// scalar infront of RHS
	std::complex<double> rhs_scalar;

	/// alpha/beta (time step related component for implicit solver)
	double alpha;
	double beta;

	bool use_formulation_with_coriolis_effect;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// Coriolis omega
	double coriolis_omega;

	/// 2*\Omega
	double two_omega;

	/// Average geopotential
	double avg_geopotential;

	std::complex<double> I = 1;


public:
	SWEImplicit_SPHRobert(SphereOperators &i_op)	:
		op(i_op),
		sphereDataConfig(nullptr),
		sphereDataConfigSolver(nullptr)
	{
	}


	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup(
			SphereDataConfig *i_sphereDataConfig,
			SphereDataConfig *i_sphereDataConfigSolver,

			double i_radius,
			double i_coriolis_omega,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_use_formulation_with_coriolis_effect = true
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereDataConfigSolver = i_sphereDataConfigSolver;

		use_formulation_with_coriolis_effect = i_use_formulation_with_coriolis_effect;
		timestep_size = i_timestep_size;

		alpha = -1.0/timestep_size;
		beta = -1.0/timestep_size;

		r = i_radius;
		inv_r = 1.0/r;

		coriolis_omega = i_coriolis_omega;

		two_omega = 2.0*coriolis_omega;
		avg_geopotential = i_avg_geopotential;

		sphSolverPhi.setup(sphereDataConfigSolver, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);

		if (use_formulation_with_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega*alpha*alpha, r);
			sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);

			sphSolverPhi.solver_component_rexi_z4robert(	-avg_geopotential*alpha*two_omega, r);
			sphSolverPhi.solver_component_rexi_z5robert(	avg_geopotential/alpha*two_omega*two_omega*two_omega, r);
			sphSolverPhi.solver_component_rexi_z6robert(	avg_geopotential*2.0*two_omega*two_omega, r);
		}

		sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential*alpha*alpha, r);
		if (use_formulation_with_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z8(	-avg_geopotential*two_omega*two_omega, r);
		}



		sphSolverVel.setup(sphereDataConfigSolver, 2);
		sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
		if (use_formulation_with_coriolis_effect)
		{
			sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);
		}
	}




	SphereData kappa(
			const SphereData &i_data
	)	const
	{
		return (alpha*alpha)*i_data + two_omega*two_omega*op.mu2(i_data);
	}


	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	void solve(
			const SphereData &i_phi0,
			const SphereData &i_u0,
			const SphereData &i_v0,

			SphereData &o_phi,
			SphereData &o_u,
			SphereData &o_v,

			double &o_dt
	)
	{
		o_dt = timestep_size;

		const SphereData &phi0 = i_phi0;
		const SphereData &u0 = i_u0;
		const SphereData &v0 = i_v0;

		SphereData div0 = inv_r*op.robert_div(u0, v0);

		SphereData eta0 = inv_r*op.robert_vort(u0, v0);

		SphereData phi(sphereDataConfig);
		SphereData u(sphereDataConfig);
		SphereData v(sphereDataConfig);

		if (use_formulation_with_coriolis_effect)
		{
#if 1
			/**
			 * Both versions (this and the version below) results of similar accuracy
			 */
			// only valid for Robert formulation!
			SphereData Fc_k =	two_omega*inv_r*(
										-(alpha*alpha*u0 - two_omega*two_omega*op.mu2(u0)) +
										2.0*alpha*two_omega*op.mu(v0)
									);

			SphereData foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*op.mu(eta0)) +
										(alpha*phi0 + two_omega*two_omega*(1.0/alpha)*op.mu2(phi0));

			SphereData rhs =	alpha*alpha*foo +
									two_omega*two_omega*op.mu2(foo)
									- (avg_geopotential/alpha)*Fc_k;

#else

			double fj = inv_r*two_omega;
			double phi_bar = avg_geopotential;

			SphereData f(sphereDataConfig);
			f.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, std::complex<double> &o_data)
					{
						o_data = mu*two_omega;
					}
				);

			SphereData Fp_i = fj*(-(alpha*alpha-f*f));
			SphereData Fp_j = fj*(2.0*alpha*f);

			SphereData Fck = Fp_i*u0 + Fp_j*v0;

			SphereData rhs =
					kappa(
							phi_bar*(div0 - f*(1.0/alpha)*eta0)
							+ (alpha + f*f*(1.0/alpha))*phi0
					)
					- phi_bar/alpha*Fck;
#endif

			phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

			SphereData a = u0 + inv_r*op.robert_grad_lon(phi);
			SphereData b = v0 + inv_r*op.robert_grad_lat(phi);

			SphereData rhsa = alpha*a - two_omega*op.mu(b);
			SphereData rhsb = two_omega*op.mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
		}
		else
		{
			SphereData rhs = avg_geopotential*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			u = (1.0/alpha) * (u0 + inv_r*op.robert_grad_lon(phi));
			v = (1.0/alpha) * (v0 + inv_r*op.robert_grad_lat(phi));
		}

		o_phi = phi * beta;
		o_u = u * beta;
		o_v = v * beta;
	}



	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	void solve_advection(
			const SphereData &i_phi0,
			const SphereData &i_u0,
			const SphereData &i_v0,

			SphereData &o_phi,

			double &o_dt
	)
	{
		o_dt = timestep_size;

		const SphereData &phi0 = i_phi0;
		const SphereData &u0 = i_u0;
		const SphereData &v0 = i_v0;

		SphereData div0 = inv_r*op.robert_div(u0, v0);
		SphereData eta0 = inv_r*op.robert_vort(u0, v0);

		SphereData phi(sphereDataConfig);
/*
		SphereData u(sphereDataConfig);
		SphereData v(sphereDataConfig);

		if (use_formulation_with_coriolis_effect)
		{
			*
			 * Both versions (this and the version below) results of similar accuracy

			// only valid for Robert formulation!
			SphereData Fc_k =	two_omega*inv_r*(
										-(alpha*alpha*u0 - two_omega*two_omega*op.mu2(u0)) +
										2.0*alpha*two_omega*op.mu(v0)
									);

			SphereData foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*op.mu(eta0)) +
										(alpha*phi0 + two_omega*two_omega*(1.0/alpha)*op.mu2(phi0));

			SphereData rhs =	alpha*alpha*foo +
									two_omega*two_omega*op.mu2(foo)
									- (avg_geopotential/alpha)*Fc_k;

			phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
		}
		else
*/
		{
			SphereData rhs = avg_geopotential*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);
		}

		o_phi = phi * beta;
	}
};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
