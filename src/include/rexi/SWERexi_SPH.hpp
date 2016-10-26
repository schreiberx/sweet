/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SWEREXI_SPH_HPP_
#define SRC_SWEREXI_SPH_HPP_

#include <sweet/sphere/SphBandedMatrixComplex.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>


/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWERexi_SPH
{
	/// SPH configuration
	SphereDataConfig *sphConfig;

	/// Solver for given alpha
	SphBandedMatrixComplex< std::complex<double> > sphSolverPhi;
	SphBandedMatrixComplex< std::complex<double> > sphSolverVel;

	/// scalar infront of RHS
	std::complex<double> rhs_scalar;

	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;

	bool include_coriolis_effect;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double ir;

	/// Coriolis omega
	double coriolis_omega;

	/// 2*\Omega
	double two_omega;

	/// Average geopotential
	double avg_geopotential;

public:
	SWERexi_SPH()	:
		sphConfig(nullptr)
	{
	}


	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup(
			SphereDataConfig *i_sphConfig,
			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,
			double i_radius,
			double i_coriolis_omega,
			double i_avg_geopotential,
			double i_timestep_size,
			bool i_include_coriolis_effect = true
	)
	{
		include_coriolis_effect = i_include_coriolis_effect;
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		ir = 1.0/r;

		coriolis_omega = i_coriolis_omega;

#if 0
		/*
		 * Setting this to zero is useful to check the
		 * (a*a+f*f)^{-1} oriented solver without any
		 * Coriolis effect it in
		 */
		// TODO: REMOVE ME!!!!
		// TODO: REMOVE ME!!!!
		// TODO: REMOVE ME!!!!
		coriolis_omega = 0;
#endif

		two_omega = 2.0*coriolis_omega;
		avg_geopotential = i_avg_geopotential;

		sphConfig = i_sphConfig;

		sphSolverPhi.setup(sphConfig, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);

		if (include_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega*alpha*alpha, r);
			sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);
			sphSolverPhi.solver_component_rexi_z4(	-avg_geopotential*alpha*two_omega, r);
			sphSolverPhi.solver_component_rexi_z5(	avg_geopotential/alpha*two_omega*two_omega*two_omega, r);
			sphSolverPhi.solver_component_rexi_z6(	avg_geopotential*2.0*two_omega*two_omega, r);
		}
		sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential*alpha*alpha, r);
		if (include_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z8(	-avg_geopotential*two_omega*two_omega, r);
		}

		sphSolverVel.setup(sphConfig, 2);
		sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
		if (include_coriolis_effect)
		{
			sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);
		}
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
			SphereData &o_v
	)
	{
		// TODO: replace with spectral operation
		SphereDataComplex mu(i_phi0.sphereDataConfig);
		mu.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, std::complex<double> &o_data)
				{
					o_data = mu;
				}
			);

		SphereDataComplex phi0(i_phi0);
		SphereDataComplex u0(i_u0);
		SphereDataComplex v0(i_v0);

		SphereDataComplex div0 = ir*SphereOperatorsComplex::div(u0, v0);
		SphereDataComplex eta0 = ir*SphereOperatorsComplex::vort(u0, v0);

		SphereDataComplex phi(sphConfig);
		SphereDataComplex u(sphConfig);
		SphereDataComplex v(sphConfig);

		if (include_coriolis_effect)
		{

#if 0
			// only works for Robert formulation!
			SphereDataComplex tmp = (
					-(alpha*alpha*u0 - two_omega*two_omega*SphereOperatorsComplex::mu2(u0)) +
					2.0*alpha*two_omega*SphereOperatorsComplex::mu(v0)
				);

			SphereDataComplex Fc_k =	two_omega*ir*(tmp-SphereOperatorsComplex::mu2(tmp));

#else

			SphereDataComplex Fc_k =	two_omega*ir*SphereOperatorsComplex::grad_lat(mu)*(
										-(alpha*alpha*u0 - two_omega*two_omega*SphereOperatorsComplex::mu2(u0)) +
										2.0*alpha*two_omega*SphereOperatorsComplex::mu(v0)
									);

#endif

			SphereDataComplex foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*SphereOperatorsComplex::mu(eta0)) +
									(alpha*i_phi0 + two_omega*two_omega*(1.0/alpha)*SphereOperatorsComplex::mu2(i_phi0));

			SphereDataComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*SphereOperatorsComplex::mu2(foo) +
									(avg_geopotential/alpha)*Fc_k;


			phi = sphSolverPhi.solve(rhs);

			SphereDataComplex a = u0 + ir*SphereOperatorsComplex::grad_lon(phi);
			SphereDataComplex b = v0 + ir*SphereOperatorsComplex::grad_lat(phi);

			SphereDataComplex rhsa = alpha*a - two_omega*SphereOperatorsComplex::mu(b);
			SphereDataComplex rhsb = two_omega*SphereOperatorsComplex::mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa);
			v = sphSolverVel.solve(rhsb);
		}
		else
		{
#if 0

			SphereDataComplex foo = avg_geopotential*div0 + alpha*i_phi0;
			SphereDataComplex rhs = (alpha*alpha)*foo;

			phi = sphSolverPhi.solve(rhs);


			SphereDataComplex rhsa = alpha*(u0 + ir*SphereOperatorsComplex::grad_lon(phi));
			SphereDataComplex rhsb = alpha*(v0 + ir*SphereOperatorsComplex::grad_lat(phi));

			// sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
			u = sphSolverVel.solve(rhsa);
			v = sphSolverVel.solve(rhsb);

#else

			SphereDataComplex rhs = avg_geopotential*div0 + alpha*i_phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			u = (1.0/alpha) * (u0 + ir*SphereOperatorsComplex::grad_lon(phi));
			v = (1.0/alpha) * (v0 + ir*SphereOperatorsComplex::grad_lat(phi));

#endif
		}

		phi *= beta;
		u *= beta;
		v *= beta;

		phi.physical_RealToSphereData(o_phi);
		u.physical_RealToSphereData(o_u);
		v.physical_RealToSphereData(o_v);
	}
};


#endif /* SRC_SWEREXI_SPH_HPP_ */
