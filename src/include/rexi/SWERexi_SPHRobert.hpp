/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SWEREXI_SPHROBERT_HPP_
#define SRC_SWEREXI_SPHROBERT_HPP_

#include <complex>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>

#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>
#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/SphBandedMatrixPhysicalComplex.hpp>



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWERexi_SPHRobert
{
	/// SPH configuration
	SphereDataConfig *sphereDataConfig;

	/// SPH configuration
	SphereDataConfig *sphereDataConfigSolver;

	/// Solver for given alpha
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverPhi;
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverVel;

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
	SWERexi_SPHRobert()	:
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
			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,
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

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

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




	SphereDataComplex kappa(
			const SphereDataComplex &i_data
	)	const
	{
		return (alpha*alpha)*i_data + two_omega*two_omega*SphereOperatorsComplex::mu2(i_data);
	}


	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	void solve_complex(
			const SphereDataComplex &i_phi0,
			const SphereDataComplex &i_u0,
			const SphereDataComplex &i_v0,

			SphereDataComplex &o_phi,
			SphereDataComplex &o_u,
			SphereDataComplex &o_v
	)
	{
		const SphereDataComplex &phi0 = i_phi0;
		const SphereDataComplex &u0 = i_u0;
		const SphereDataComplex &v0 = i_v0;

		SphereDataComplex div0 = inv_r*SphereOperatorsComplex::robert_div(u0, v0);
		SphereDataComplex eta0 = inv_r*SphereOperatorsComplex::robert_vort(u0, v0);

		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex u(sphereDataConfig);
		SphereDataComplex v(sphereDataConfig);

		if (use_formulation_with_coriolis_effect)
		{
#if 1
			/**
			 * Both versions (this and the version below) results of similar accuracy
			 */
			// only valid for Robert formulation!
			SphereDataComplex Fc_k =	two_omega*inv_r*(
										-(alpha*alpha*u0 - two_omega*two_omega*SphereOperatorsComplex::mu2(u0)) +
										2.0*alpha*two_omega*SphereOperatorsComplex::mu(v0)
									);

			SphereDataComplex foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*SphereOperatorsComplex::mu(eta0)) +
										(alpha*phi0 + two_omega*two_omega*(1.0/alpha)*SphereOperatorsComplex::mu2(phi0));

			SphereDataComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*SphereOperatorsComplex::mu2(foo)
									- (avg_geopotential/alpha)*Fc_k;

#else

			double fj = inv_r*two_omega;
			double phi_bar = avg_geopotential;

			SphereDataComplex f(sphereDataConfig);
			f.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, std::complex<double> &o_data)
					{
						o_data = mu*two_omega;
					}
				);

			SphereDataComplex Fp_i = fj*(-(alpha*alpha-f*f));
			SphereDataComplex Fp_j = fj*(2.0*alpha*f);

			SphereDataComplex Fck = Fp_i*u0 + Fp_j*v0;

			SphereDataComplex rhs =
					kappa(
							phi_bar*(div0 - f*(1.0/alpha)*eta0)
							+ (alpha + f*f*(1.0/alpha))*phi0
					)
					- phi_bar/alpha*Fck;

#endif

			phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

			SphereDataComplex a = u0 + inv_r*SphereOperatorsComplex::robert_grad_lon(phi);
			SphereDataComplex b = v0 + inv_r*SphereOperatorsComplex::robert_grad_lat(phi);

			SphereDataComplex rhsa = alpha*a - two_omega*SphereOperatorsComplex::mu(b);
			SphereDataComplex rhsb = two_omega*SphereOperatorsComplex::mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
		}
		else
		{
			SphereDataComplex rhs = avg_geopotential*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			u = (1.0/alpha) * (u0 + inv_r*SphereOperatorsComplex::robert_grad_lon(phi));
			v = (1.0/alpha) * (v0 + inv_r*SphereOperatorsComplex::robert_grad_lat(phi));
		}

//		std::cout << beta << std::endl;

		o_phi = phi * beta;
		o_u = u * beta;
		o_v = v * beta;
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
		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex u(sphereDataConfig);
		SphereDataComplex v(sphereDataConfig);

		solve_complex(
				i_phi0,
				i_u0,
				i_v0,
				phi,
				u,
				v
			);

		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(phi);
		o_u = Convert_SphereDataComplex_To_SphereData::physical_convert(u);
		o_v = Convert_SphereDataComplex_To_SphereData::physical_convert(v);

#if 0
		static int i = 0;
		i++;

		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		char buffer[1024];

		sprintf(buffer, "data_phi_%04i_%f_%f.csv", i, alpha.real(), alpha.imag());
		(Convert_SphereDataComplex_To_SphereData::physical_convert(phi*alpha-i_phi0)).physical_file_write(buffer);

		sprintf(buffer, "data_u_%04i_%f_%f.csv", i, alpha.real(), alpha.imag());
		(Convert_SphereDataComplex_To_SphereData::physical_convert(u*alpha-i_u0)).physical_file_write(buffer);

		sprintf(buffer, "data_v_%04i_%f_%f.csv", i, alpha.real(), alpha.imag());
		(Convert_SphereDataComplex_To_SphereData::physical_convert(v*alpha-i_v0)).physical_file_write(buffer);
#endif

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

};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
