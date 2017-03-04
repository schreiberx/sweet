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

		op.setup(sphereDataConfig, r);
		opComplex.setup(sphereDataConfig, r);

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
		return (alpha*alpha)*i_data + two_omega*two_omega*opComplex.mu2(i_data);
	}


	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	inline
	void solve_advection_progphiuv(
			const SphereData &i_phi0,
			const SphereData &i_u0,
			const SphereData &i_v0,

			SphereData &o_phi,
			SphereData &o_u,
			SphereData &o_v
	)
	{
		const SphereDataComplex &phi0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_phi0);
		const SphereDataComplex &u0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_u0);
		const SphereDataComplex &v0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_v0);

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
			SphereDataComplex Fc_k =	two_omega*inv_r*(
										-(alpha*alpha*u0 - two_omega*two_omega*opComplex.mu2(u0)) +
										2.0*alpha*two_omega*opComplex.mu(v0)
									);

			SphereDataComplex foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*opComplex.mu(eta0)) +
										(alpha*phi0 + two_omega*two_omega*(1.0/alpha)*opComplex.mu2(phi0));

			SphereDataComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*opComplex.mu2(foo)
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
	void solve_advection_progphiuv_ver2(
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

#if 0
		{
			SphereDataPhysical u0gr(sphereDataConfig);
			SphereDataPhysical v0gr(sphereDataConfig);
			op.robert_vortdiv_to_uv(i_eta0, i_div0, u0gr, v0gr);

			SphereData u(sphereDataConfig);
			SphereData v(sphereDataConfig);

			solve_advection_progphiuv(
					i_phi0, (SphereData)u0gr, (SphereData)v0gr,
					o_phi, u, v
			);

			op.robert_uv_to_vortdiv(
					u.getSphereDataPhysical(),
					v.getSphereDataPhysical(),
					o_eta,
					o_div
				);
		}

		return;
#endif

		const SphereDataPhysicalComplex phi0g = i_phi0.getSphereDataPhysicalComplex();
		const SphereDataPhysicalComplex eta0g = i_eta0.getSphereDataPhysicalComplex();
		const SphereDataPhysicalComplex div0g = i_div0.getSphereDataPhysicalComplex();

		SphereDataPhysicalComplex fg(sphereDataConfig);
		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, std::complex<double> &o_data)
			{
				o_data = mu;
			}
		);


		SphereDataPhysical u0gr(sphereDataConfig);
		SphereDataPhysical v0gr(sphereDataConfig);
		op.robert_vortdiv_to_uv(i_eta0, i_div0, u0gr, v0gr);

		SphereDataPhysicalComplex u0g = u0gr;
		SphereDataPhysicalComplex v0g = v0gr;

		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex eta(sphereDataConfig);
		SphereDataComplex div(sphereDataConfig);

		SphereDataComplex u(sphereDataConfig);
		SphereDataComplex v(sphereDataConfig);

		if (use_formulation_with_coriolis_effect)
		{
			/**
			 * Both versions (this and the version below) results of similar accuracy
			 */
			// only valid for Robert formulation!

#if 1
			SphereDataPhysicalComplex Fc_k = two_omega*inv_r*(
										-(alpha*alpha*u0g - two_omega*two_omega*fg*fg*u0g) +
										2.0*alpha*two_omega*fg*v0g
									);

			SphereDataPhysicalComplex foo = 	avg_geopotential*(div0g - two_omega*(1.0/alpha)*fg*eta0g) +
										(alpha*phi0g + two_omega*two_omega*(1.0/alpha)*fg*fg*phi0g);

			SphereDataPhysicalComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*fg*fg*foo
									- (avg_geopotential/alpha)*Fc_k;

#else
			SphereDataComplex Fc_k =	two_omega*inv_r*(
										-(alpha*alpha*(SphereDataComplex)u0g - two_omega*two_omega*opComplex.mu2(u0g)) +
										2.0*alpha*two_omega*opComplex.mu(v0g)
									);

			SphereDataComplex foo = 	avg_geopotential*((SphereDataComplex)div0g - two_omega*(1.0/alpha)*opComplex.mu(eta0)) +
										(alpha*(SphereDataComplex)phi0g + two_omega*two_omega*(1.0/alpha)*opComplex.mu2(phi0g));

			SphereDataComplex rhs =	alpha*alpha*foo +
									two_omega*two_omega*opComplex.mu2(foo)
									- (avg_geopotential/alpha)*Fc_k;
#endif

			phi = sphSolverPhi.solve(((SphereDataComplex)rhs).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

#if 0

			SphereDataPhysicalComplex a = u0g + inv_r*(opComplex.robert_grad_lon(phi)).getSphereDataPhysicalComplex();
			SphereDataPhysicalComplex b = v0g + inv_r*(opComplex.robert_grad_lat(phi)).getSphereDataPhysicalComplex();

			SphereDataComplex rhsa = alpha*a - two_omega*fg*b;
			SphereDataComplex rhsb = two_omega*fg*a + alpha*b;

			SphereDataComplex u = sphSolverVel.solve(((SphereDataComplex)rhsa).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			SphereDataComplex v = sphSolverVel.solve(((SphereDataComplex)rhsb).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

			SphereDataPhysical ur = Convert_SphereDataComplex_To_SphereData::physical_convert(u).getSphereDataPhysical();
			SphereDataPhysical ui = Convert_SphereDataComplex_To_SphereData::physical_convert_imag(u).getSphereDataPhysical();

			SphereDataPhysical vr = Convert_SphereDataComplex_To_SphereData::physical_convert(v).getSphereDataPhysical();
			SphereDataPhysical vi = Convert_SphereDataComplex_To_SphereData::physical_convert_imag(v).getSphereDataPhysical();

			SphereData etar(sphereDataConfig);
			SphereData etai(sphereDataConfig);
			SphereData divr(sphereDataConfig);
			SphereData divi(sphereDataConfig);

			op.robert_uv_to_vortdiv(ur, vr, etar, divr);
			op.robert_uv_to_vortdiv(ui, vi, etai, divi);


			eta = Convert_SphereData_To_SphereDataComplex::physical_convert(etar) +
					Convert_SphereData_To_SphereDataComplex::physical_convert(etai)*std::complex<double>(0, 1);

			div = Convert_SphereData_To_SphereDataComplex::physical_convert(divr) +
					Convert_SphereData_To_SphereDataComplex::physical_convert(divi)*std::complex<double>(0, 1);

#else

			SphereDataComplex a = (SphereDataComplex)u0g + inv_r*opComplex.robert_grad_lon(phi);
			SphereDataComplex b = (SphereDataComplex)v0g + inv_r*opComplex.robert_grad_lat(phi);

			SphereDataComplex rhsa = alpha*a - two_omega*opComplex.mu(b);
			SphereDataComplex rhsb = two_omega*opComplex.mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
#endif
		}
		else
		{
			SphereDataComplex rhs = avg_geopotential*div0g + alpha*phi0g;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			SphereDataComplex u = (1.0/alpha) * (u0g + inv_r*opComplex.robert_grad_lon(phi).getSphereDataPhysicalComplex());
			SphereDataComplex v = (1.0/alpha) * (v0g + inv_r*opComplex.robert_grad_lat(phi).getSphereDataPhysicalComplex());

#if 1
			SphereDataPhysical ur = Convert_SphereDataComplex_To_SphereData::physical_convert(u).getSphereDataPhysical();
			SphereDataPhysical ui = Convert_SphereDataComplex_To_SphereData::physical_convert_imag(u).getSphereDataPhysical();

			SphereDataPhysical vr = Convert_SphereDataComplex_To_SphereData::physical_convert(v).getSphereDataPhysical();
			SphereDataPhysical vi = Convert_SphereDataComplex_To_SphereData::physical_convert_imag(v).getSphereDataPhysical();

			SphereData etar(sphereDataConfig);
			SphereData etai(sphereDataConfig);
			SphereData divr(sphereDataConfig);
			SphereData divi(sphereDataConfig);

			op.robert_uv_to_vortdiv(ur, vr, etar, divr);
			op.robert_uv_to_vortdiv(ui, vi, etai, divi);


			eta = Convert_SphereData_To_SphereDataComplex::physical_convert(etar) +
					Convert_SphereData_To_SphereDataComplex::physical_convert(etai)*std::complex<double>(0, 1);

			div = Convert_SphereData_To_SphereDataComplex::physical_convert(divr) +
					Convert_SphereData_To_SphereDataComplex::physical_convert(divi)*std::complex<double>(0, 1);
#else
			op.robert_uv_to_vortdiv(u, v, eta, div);
#endif
		}

		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert(phi * beta);
//		o_eta = Convert_SphereDataComplex_To_SphereData::physical_convert(eta * beta);
//		o_div = Convert_SphereDataComplex_To_SphereData::physical_convert(div * beta);

		SphereDataPhysical ur = Convert_SphereDataComplex_To_SphereData::physical_convert(u * beta).getSphereDataPhysical();
		SphereDataPhysical vr = Convert_SphereDataComplex_To_SphereData::physical_convert(v * beta).getSphereDataPhysical();

		op.robert_uv_to_vortdiv(ur, vr, o_eta, o_div);
	}


};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
