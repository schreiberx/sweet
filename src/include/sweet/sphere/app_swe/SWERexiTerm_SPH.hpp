/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_SWEREXI_SPH_HPP_
#define SRC_SWEREXI_SPH_HPP_

#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>

#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>
#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalComplex.hpp>


/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWERexiTerm_SPH
{
	/// SPH configuration
	SphereDataConfig *sphereDataConfig;

	/// Solver for given alpha
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverPhi;
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverVel;

	SphereOperatorsComplex opComplex;


	/// scalar infront of RHS
	std::complex<double> rhs_scalar;

	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;

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

public:
	SWERexiTerm_SPH()	:
		sphereDataConfig(nullptr)
	{
	}


	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup(
			SphereDataConfig *i_sphereDataConfig,

			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,

			double i_radius,
			double i_coriolis_omega,
			double i_avg_geopotential,

			double i_timestep_size,

			bool i_include_coriolis_effect = true
	)
	{
		use_formulation_with_coriolis_effect = i_include_coriolis_effect;
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		inv_r = 1.0/r;

		coriolis_omega = i_coriolis_omega;

		two_omega = 2.0*coriolis_omega;
		avg_geopotential = i_avg_geopotential;

		sphereDataConfig = i_sphereDataConfig;

		opComplex.setup(sphereDataConfig, r);

		sphSolverPhi.setup(sphereDataConfig, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);

		if (use_formulation_with_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega*alpha*alpha, r);
			sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);
			sphSolverPhi.solver_component_rexi_z4(	-avg_geopotential*alpha*two_omega, r);
			sphSolverPhi.solver_component_rexi_z5(	avg_geopotential/alpha*two_omega*two_omega*two_omega, r);
			sphSolverPhi.solver_component_rexi_z6(	avg_geopotential*2.0*two_omega*two_omega, r);
		}
		sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential*alpha*alpha, r);
		if (use_formulation_with_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z8(	-avg_geopotential*two_omega*two_omega, r);
		}

		sphSolverVel.setup(sphereDataConfig, 2);
		sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
		if (use_formulation_with_coriolis_effect)
		{
			sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);
		}
	}




	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	inline
	void solve(
			const SphereData &i_phi0,
			const SphereData &i_u0,
			const SphereData &i_v0,

			SphereData &o_phi,
			SphereData &o_u,
			SphereData &o_v
	)
	{
		solve_complexRHS(
				Convert_SphereData_To_SphereDataComplex::physical_convert(i_phi0),
				Convert_SphereData_To_SphereDataComplex::physical_convert(i_u0),
				Convert_SphereData_To_SphereDataComplex::physical_convert(i_v0),

				o_phi,
				o_u,
				o_v
			);
	}



	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	void solve_complexRHS(
			const SphereDataComplex &i_phi0,
			const SphereDataComplex &i_u0,
			const SphereDataComplex &i_v0,

			SphereData &o_phi,
			SphereData &o_u,
			SphereData &o_v
	)
	{
#if 1
		// TODO: replace with spectral operation
		SphereDataComplex mu(i_phi0.sphereDataConfig);
		mu.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, std::complex<double> &o_data)
				{
					o_data = mu;
				}
			);
#endif

		const SphereDataComplex &phi0 = i_phi0;
		const SphereDataComplex &u0 = i_u0;
		const SphereDataComplex &v0 = i_v0;

		SphereDataComplex div0 = inv_r*opComplex.div(u0, v0);
		SphereDataComplex eta0 = inv_r*opComplex.vort(u0, v0);

		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex u(sphereDataConfig);
		SphereDataComplex v(sphereDataConfig);

		if (use_formulation_with_coriolis_effect)
		{

#if 0
			// only works for Robert formulation!
			SphereDataComplex tmp = (
					-(alpha*alpha*i_u0 - two_omega*two_omega*opComplex.mu2(i_u0)) +
					2.0*alpha*two_omega*opComplex.mu(i_v0)
				);

			SphereDataComplex Fc_k =	two_omega*ir*(tmp-opComplex.mu2(tmp));

#else

			SphereDataComplex Fc_k =	two_omega*inv_r*opComplex.grad_lat(mu)*(
										-(alpha*alpha*i_u0 - two_omega*two_omega*opComplex.mu2(i_u0)) +
										2.0*alpha*two_omega*opComplex.mu(i_v0)
									);

#endif

			SphereDataComplex foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*opComplex.mu(eta0)) +
									(alpha*i_phi0 + two_omega*two_omega*(1.0/alpha)*opComplex.mu2(i_phi0));

			SphereDataComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*opComplex.mu2(foo) +
									(avg_geopotential/alpha)*Fc_k;


			phi = sphSolverPhi.solve(rhs);

			SphereDataComplex a = i_u0 + inv_r*opComplex.grad_lon(phi);
			SphereDataComplex b = i_v0 + inv_r*opComplex.grad_lat(phi);

			SphereDataComplex rhsa = alpha*a - two_omega*opComplex.mu(b);
			SphereDataComplex rhsb = two_omega*opComplex.mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa);
			v = sphSolverVel.solve(rhsb);
		}
		else
		{
			SphereDataComplex rhs = avg_geopotential*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			u = (1.0/alpha) * (u0 + inv_r*opComplex.grad_lon(phi));
			v = (1.0/alpha) * (v0 + inv_r*opComplex.grad_lat(phi));
		}

		phi *= beta;
		u *= beta;
		v *= beta;

		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(phi);
		o_u = Convert_SphereDataComplex_To_SphereData::physical_convert(u);
		o_v = Convert_SphereDataComplex_To_SphereData::physical_convert(v);
	}
};


#endif /* SRC_SWEREXI_SPH_HPP_ */
