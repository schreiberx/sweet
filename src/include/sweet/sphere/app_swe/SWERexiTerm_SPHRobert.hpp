/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_SWEREXI_SPHROBERT_HPP_
#define SRC_SWEREXI_SPHROBERT_HPP_

#include <complex>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereDataPhysicalComplex.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>

#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>
#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalComplex.hpp>



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWERexiTerm_SPHRobert
{
	/// SPH configuration
	SphereDataConfig *sphereDataConfig;

	/// SPH configuration
	SphereDataConfig *sphereDataConfigSolver;

	/// Solver for given alpha
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverPhi;
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverVel;

	SphereOperators op;
	SphereOperatorsComplex opComplex;

	/// scalar infront of RHS
	std::complex<double> rhs_scalar;

	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;

	bool use_formulation_with_coriolis_effect;

	bool use_f_sphere;

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
	SWERexiTerm_SPHRobert()	:
		sphereDataConfig(nullptr),
		sphereDataConfigSolver(nullptr)
	{
	}


	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup_velocityformulation_progphiuv(
			SphereDataConfig *i_sphereDataConfig,
			SphereDataConfig *i_sphereDataConfigSolver,

			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,

			double i_radius,
			double i_coriolis_omega,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_use_formulation_with_coriolis_effect,
			bool i_use_f_sphere
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereDataConfigSolver = i_sphereDataConfigSolver;

		use_formulation_with_coriolis_effect = i_use_formulation_with_coriolis_effect;
		use_f_sphere = i_use_f_sphere;
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		inv_r = 1.0/r;

		coriolis_omega = i_coriolis_omega;

		two_omega = 2.0*coriolis_omega;
		avg_geopotential = i_avg_geopotential;

		op.setup(sphereDataConfig, r);
		opComplex.setup(sphereDataConfig, r);

		sphSolverPhi.setup(sphereDataConfigSolver, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);

		if (use_formulation_with_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega*alpha*alpha, r);
			sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);

			if (!use_f_sphere)
			{
				sphSolverPhi.solver_component_rexi_z4robert(	-avg_geopotential*alpha*two_omega, r);
				sphSolverPhi.solver_component_rexi_z5robert(	avg_geopotential/alpha*two_omega*two_omega*two_omega, r);
				sphSolverPhi.solver_component_rexi_z6robert(	avg_geopotential*2.0*two_omega*two_omega, r);
			}
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



	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup_vectorinvariant_progphivortdiv(
			SphereDataConfig *i_sphereDataConfig,
			SphereDataConfig *i_sphereDataConfigSolver,

			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,

			double i_radius,
			double i_coriolis_omega,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_use_formulation_with_coriolis_effect,
			bool i_use_f_sphere
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereDataConfigSolver = i_sphereDataConfigSolver;

		use_formulation_with_coriolis_effect = i_use_formulation_with_coriolis_effect;
		use_f_sphere = i_use_f_sphere;
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		inv_r = 1.0/r;

		coriolis_omega = i_coriolis_omega;

		two_omega = 2.0*coriolis_omega;
		avg_geopotential = i_avg_geopotential;

		double gh = i_avg_geopotential;

		op.setup(sphereDataConfig, r);
		opComplex.setup(sphereDataConfig, r);

		if (use_formulation_with_coriolis_effect && !use_f_sphere)
			sphSolverPhi.setup(sphereDataConfigSolver, 2);
		else
			sphSolverPhi.setup(sphereDataConfigSolver, 0);

		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha), r);

		if (use_formulation_with_coriolis_effect)
		{
			if (use_f_sphere)
				sphSolverPhi.solver_component_rexi_z1(	coriolis_omega*coriolis_omega, r);
			else
				sphSolverPhi.solver_component_rexi_z2(	two_omega*two_omega, r);
		}

		sphSolverPhi.solver_component_rexi_z7(	-gh, r);

//		std::cout << use_formulation_with_coriolis_effect << "\t" << use_f_sphere << "\t" << two_omega << "\t" << avg_geopotential << "\t" << r << "\t" << alpha << std::endl;

	}


	SphereDataComplex kappa(
			const SphereDataComplex &i_data
	)	const
	{
		return (alpha*alpha)*i_data + two_omega*two_omega*opComplex.mu2(i_data);
	}


	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	inline
	void solve_velocityformulation_progphiuv(
			const SphereData &i_phi0,
			const SphereData &i_u0,
			const SphereData &i_v0,

			SphereData &o_phi,
			SphereData &o_u,
			SphereData &o_v
	)
	{
		const SphereDataComplex phi0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_phi0);
		const SphereDataComplex u0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_u0);
		const SphereDataComplex v0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_v0);

#if 1
		SphereDataComplex div0 = inv_r*opComplex.robert_div(u0, v0);
		SphereDataComplex eta0 = inv_r*opComplex.robert_vort(u0, v0);
#else

		SphereData div0r(sphereDataConfig);
		SphereData eta0r(sphereDataConfig);
		op.robert_uv_to_vortdiv(i_u0.getSphereDataPhysical(), i_v0.getSphereDataPhysical(), eta0r, div0r);
		SphereDataComplex eta0 = Convert_SphereData_To_SphereDataComplex::physical_convert(eta0r);
		SphereDataComplex div0 = Convert_SphereData_To_SphereDataComplex::physical_convert(div0r);
#endif
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
			SphereDataComplex foo = avg_geopotential*(div0 - two_omega*(1.0/alpha)*opComplex.mu(eta0)) +
									(alpha*phi0 + two_omega*two_omega*(1.0/alpha)*opComplex.mu2(phi0));

			SphereDataComplex rhs(sphereDataConfig);

			if (use_f_sphere)
			{
				rhs =	alpha*alpha*foo +
						two_omega*two_omega*opComplex.mu2(foo);
			}
			else
			{
				SphereDataComplex Fc_k =	two_omega*inv_r*(
							-(alpha*alpha*u0 - two_omega*two_omega*opComplex.mu2(u0)) +
							2.0*alpha*two_omega*opComplex.mu(v0)
						);

				rhs =	alpha*alpha*foo +
						two_omega*two_omega*opComplex.mu2(foo)
						- (avg_geopotential/alpha)*Fc_k;
			}

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

			SphereDataComplex a = u0 + inv_r*opComplex.robert_grad_lon(phi);
			SphereDataComplex b = v0 + inv_r*opComplex.robert_grad_lat(phi);

			SphereDataComplex rhsa = alpha*a - two_omega*opComplex.mu(b);
			SphereDataComplex rhsb = two_omega*opComplex.mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
		}
		else
		{
			SphereDataComplex rhs = avg_geopotential*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			u = (1.0/alpha) * (u0 + inv_r*opComplex.robert_grad_lon(phi));
			v = (1.0/alpha) * (v0 + inv_r*opComplex.robert_grad_lat(phi));
		}

		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(phi * beta);
		o_u = Convert_SphereDataComplex_To_SphereData::physical_convert(u * beta);
		o_v = Convert_SphereDataComplex_To_SphereData::physical_convert(v * beta);
	}





	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	inline
	void solve_advection_progphiuv_ver2xy(
			const SphereData &i_phi0,
			const SphereData &i_u0,
			const SphereData &i_v0,

			SphereData &o_phi,
			SphereData &o_u,
			SphereData &o_v
	)
	{
		const SphereDataPhysicalComplex phi0 = i_phi0.getSphereDataPhysicalComplex();
		const SphereDataPhysicalComplex u0 = i_u0.getSphereDataPhysicalComplex();
		const SphereDataPhysicalComplex v0 = i_v0.getSphereDataPhysicalComplex();

		SphereDataPhysicalComplex fg(sphereDataConfig);
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, std::complex<double> &o_data)
			{
				o_data = mu;
			}
		);

#if 1

		SphereDataComplex div0r = inv_r*opComplex.robert_div((SphereDataComplex)u0, v0);
		SphereDataComplex eta0r = inv_r*opComplex.robert_vort(u0, v0);

#else
		SphereData eta0r(sphereDataConfig);
		SphereData div0r(sphereDataConfig);
		op.robert_uv_to_vortdiv(i_u0.getSphereDataPhysical(), i_v0.getSphereDataPhysical(), eta0r, div0r);
#endif
		SphereDataPhysicalComplex eta0 = eta0r.getSphereDataPhysicalComplex();
		SphereDataPhysicalComplex div0 = div0r.getSphereDataPhysicalComplex();

		SphereDataComplex phi;
		SphereDataComplex u;
		SphereDataComplex v;


		if (use_formulation_with_coriolis_effect)
		{
			/**
			 * Both versions (this and the version below) results of similar accuracy
			 */
			// only valid for Robert formulation!

#if 1
			SphereDataPhysicalComplex Fc_k = two_omega*inv_r*(
										-(alpha*alpha*u0 - two_omega*two_omega*fg*fg*u0) +
										2.0*alpha*two_omega*fg*v0
									);

			SphereDataPhysicalComplex foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*fg*eta0) +
										(alpha*phi0 + two_omega*two_omega*(1.0/alpha)*fg*fg*phi0);

			SphereDataPhysicalComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*fg*fg*foo
									- (avg_geopotential/alpha)*Fc_k;

#else
			SphereDataComplex Fc_k =	two_omega*inv_r*(
										-(alpha*alpha*(SphereDataComplex)u0 - two_omega*two_omega*opComplex.mu2(u0)) +
										2.0*alpha*two_omega*opComplex.mu(v0)
									);

			SphereDataComplex foo = 	avg_geopotential*((SphereDataComplex)div0 - two_omega*(1.0/alpha)*opComplex.mu(eta0)) +
										(alpha*(SphereDataComplex)phi0 + two_omega*two_omega*(1.0/alpha)*opComplex.mu2(phi0));

			SphereDataComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*opComplex.mu2(foo)
									- (avg_geopotential/alpha)*Fc_k;
#endif

			phi = sphSolverPhi.solve(((SphereDataComplex)rhs).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

#if 1
			SphereDataPhysicalComplex a = u0 + inv_r*(opComplex.robert_grad_lon(phi)).getSphereDataPhysicalComplex();
			SphereDataPhysicalComplex b = v0 + inv_r*(opComplex.robert_grad_lat(phi)).getSphereDataPhysicalComplex();

			SphereDataComplex rhsa = alpha*a - two_omega*fg*b;
			SphereDataComplex rhsb = two_omega*fg*a + alpha*b;

			u = sphSolverVel.solve(((SphereDataComplex)rhsa).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(((SphereDataComplex)rhsb).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

#else

			SphereDataComplex a = (SphereDataComplex)u0 + inv_r*opComplex.robert_grad_lon(phi);
			SphereDataComplex b = (SphereDataComplex)v0 + inv_r*opComplex.robert_grad_lat(phi);

			SphereDataComplex rhsa = alpha*a - two_omega*opComplex.mu(b);
			SphereDataComplex rhsb = two_omega*opComplex.mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
#endif
		}
		else
		{
			SphereDataComplex rhs = avg_geopotential*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			u = (1.0/alpha) * (u0 + inv_r*opComplex.robert_grad_lon(phi).getSphereDataPhysicalComplex());
			v = (1.0/alpha) * (v0 + inv_r*opComplex.robert_grad_lat(phi).getSphereDataPhysicalComplex());
		}

		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(phi * beta);
		o_u = Convert_SphereDataComplex_To_SphereData::physical_convert(u * beta);
		o_v = Convert_SphereDataComplex_To_SphereData::physical_convert(v * beta);
	}



	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	inline
	void solve_vectorinvariant_progphivortdiv(
			const SphereData &i_phi0,
			const SphereData &i_eta0,
			const SphereData &i_div0,

			SphereData &o_phi,
			SphereData &o_eta,
			SphereData &o_div
	)
	{
		const SphereDataComplex phi0c = Convert_SphereData_To_SphereDataComplex::physical_convert(i_phi0);
		const SphereDataComplex eta0c = Convert_SphereData_To_SphereDataComplex::physical_convert(i_eta0);
		const SphereDataComplex div0c = Convert_SphereData_To_SphereDataComplex::physical_convert(i_div0);


		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex eta(sphereDataConfig);
		SphereDataComplex div(sphereDataConfig);

		if (use_formulation_with_coriolis_effect)
		{
			SphereDataComplex rhs =
								avg_geopotential*div0c
								+ alpha*phi0c;

			if (use_formulation_with_coriolis_effect)
			{
				if (use_f_sphere)
				{
					rhs = rhs
						+ (coriolis_omega*coriolis_omega/alpha)*(phi0c)
						- (avg_geopotential*coriolis_omega/alpha)*(eta0c);
				}
				else
				{
					rhs = rhs
						+ (two_omega*two_omega/alpha)*opComplex.mu2(phi0c)
						- (avg_geopotential*two_omega/alpha)*opComplex.mu(eta0c);
				}
			}

			phi = sphSolverPhi.solve(
										((SphereDataComplex)rhs).spectral_returnWithDifferentModes(sphereDataConfigSolver)
									).spectral_returnWithDifferentModes(sphereDataConfig);

			div = -1.0/avg_geopotential*(phi0c - alpha*phi);


			if (use_formulation_with_coriolis_effect)
			{
				if (use_f_sphere)
					eta = (1.0/alpha)*(eta0c + coriolis_omega*(div));
				else
					eta = (1.0/alpha)*(eta0c + two_omega*opComplex.mu(div));
			}
		}
		else
		{
#if 1

			SphereDataComplex rhs =
								avg_geopotential*div0c
								+ alpha*phi0c;

			phi = sphSolverPhi.solve(
								((SphereDataComplex)rhs).spectral_returnWithDifferentModes(sphereDataConfigSolver)
							).spectral_returnWithDifferentModes(sphereDataConfig);


			div = -1.0/avg_geopotential*(phi0c - alpha*phi);
			eta = 1.0/alpha*(eta0c);

#else

			SphereDataComplex rhs = avg_geopotential*div0c + alpha*phi0c;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			div = -1.0/avg_geopotential*(phi0c - alpha*phi);

			eta = 1.0/alpha*(eta0c);

#endif
		}

		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(phi * beta);
		o_eta = Convert_SphereDataComplex_To_SphereData::physical_convert(eta * beta);
		o_div = Convert_SphereDataComplex_To_SphereData::physical_convert(div * beta);
	}


};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
