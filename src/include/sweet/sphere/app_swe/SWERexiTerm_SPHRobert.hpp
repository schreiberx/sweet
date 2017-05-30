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

	/// scalar in front of RHS
	std::complex<double> rhs_scalar;

	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;

	bool use_f_sphere;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// Coriolis omega
	double coriolis;

	/// Average geopotential
	double gh;

	int variant_id;

	SphereDataPhysicalComplex mug;

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
			double i_coriolis,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_use_f_sphere
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereDataConfigSolver = i_sphereDataConfigSolver;

		use_f_sphere = i_use_f_sphere;
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		inv_r = 1.0/r;

		if (use_f_sphere)
			coriolis = i_coriolis;
		else
			coriolis = 2.0*i_coriolis;

		gh = i_avg_geopotential;

		op.setup(sphereDataConfig, r);
		opComplex.setup(sphereDataConfig, r);

		if (coriolis != 0)
		{
			if (use_f_sphere)
			{
				// NOTHING TO DO HERE
			}
			else
			{
				sphSolverPhi.setup(sphereDataConfigSolver, 4);
				sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
				sphSolverPhi.solver_component_rexi_z2(	2.0*coriolis*coriolis*alpha*alpha, r);
				sphSolverPhi.solver_component_rexi_z3(	(coriolis*coriolis)*(coriolis*coriolis), r);
				sphSolverPhi.solver_component_rexi_z4robert(	-gh*alpha*coriolis, r);
				sphSolverPhi.solver_component_rexi_z5robert(	gh/alpha*coriolis*coriolis*coriolis, r);
				sphSolverPhi.solver_component_rexi_z6robert(	gh*2.0*coriolis*coriolis, r);
				sphSolverPhi.solver_component_rexi_z7(	-gh*alpha*alpha, r);
				sphSolverPhi.solver_component_rexi_z8(	-gh*coriolis*coriolis, r);

				sphSolverVel.setup(sphereDataConfigSolver, 2);
				sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
				sphSolverVel.solver_component_rexi_z2(	coriolis*coriolis, r);
			}
		}
		else
		{
			// NOTHING TO DO HERE
		}
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

#if 0
#if 1
		SphereDataComplex div0 = inv_r*opComplex.robert_div(u0, v0);
		SphereDataComplex eta0 = inv_r*opComplex.robert_vort(u0, v0);
#else

		SphereData div0r(sphereDataConfig);
		SphereData eta0r(sphereDataConfig);
		op.robert_uv_to_vortdiv(i_u0.getSphereDataPhysical(), i_v0.getSphereDataPhysical(), eta0r, div0r);
		SphereDataComplex eta0 = Convert_SphereData_To_SphereDataComplex::physical_convert_real(eta0r);
		SphereDataComplex div0 = Convert_SphereData_To_SphereDataComplex::physical_convert_real(div0r);
#endif
		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex u(sphereDataConfig);
		SphereDataComplex v(sphereDataConfig);

		if (coriolis != 0)
		{
#if 1
			/**
			 * Both versions (this and the version below) results of similar accuracy
			 */
			SphereDataComplex rhs(sphereDataConfig);

			if (use_f_sphere)
			{
				SphereDataComplex foo = gh*(div0 - coriolis*(1.0/alpha)*eta0) +
										(alpha*phi0 + coriolis*coriolis*(1.0/alpha)*phi0);

				SphereDataComplex Fc_k =	coriolis*inv_r*(
							-(alpha*alpha*u0 - coriolis*coriolis*(u0)) +
							2.0*alpha*coriolis*(v0)
						);

				rhs =	alpha*alpha*foo +
						coriolis*coriolis*foo;
			}
			else
			{
				SphereDataComplex foo = gh*(div0 - coriolis*(1.0/alpha)*opComplex.mu(eta0)) +
										(alpha*phi0 + coriolis*coriolis*(1.0/alpha)*opComplex.mu2(phi0));

				SphereDataComplex Fc_k =	coriolis*inv_r*(
							-(alpha*alpha*u0 - coriolis*coriolis*opComplex.mu2(u0)) +
							2.0*alpha*coriolis*opComplex.mu(v0)
						);

				rhs =	alpha*alpha*foo +
						coriolis*coriolis*opComplex.mu2(foo)
						- (gh/alpha)*Fc_k;
			}

#else

			double fj = inv_r*coriolis;
			double phi_bar = gh;

			SphereDataComplex f(sphereDataConfig);
			f.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, std::complex<double> &o_data)
					{
						o_data = mu*coriolis;
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

			SphereDataComplex rhsa(sphereDataConfig);
			SphereDataComplex rhsb(sphereDataConfig);

			if (use_f_sphere)
			{
				rhsa = alpha*a - coriolis*b;
				rhsb = coriolis*a + alpha*b;
			}
			else
			{
				rhsa = alpha*a - coriolis*opComplex.mu(b);
				rhsb = coriolis*opComplex.mu(a) + alpha*b;
			}

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
		}
		else
		{
			SphereDataComplex rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

			u = (1.0/alpha) * (u0 + inv_r*opComplex.robert_grad_lon(phi));
			v = (1.0/alpha) * (v0 + inv_r*opComplex.robert_grad_lat(phi));
		}
#endif


		SphereDataComplex div0 = inv_r*opComplex.robert_div(u0, v0);
		SphereDataComplex vort0 = inv_r*opComplex.robert_vort(u0, v0);

		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex u(sphereDataConfig);
		SphereDataComplex v(sphereDataConfig);

		if (coriolis != 0)
		{
			if (use_f_sphere)
			{
				SphereDataComplex rhs = gh*(div0 - coriolis/alpha*vort0) + (alpha+coriolis*coriolis/alpha)*phi0;
				phi = rhs.spectral_solve_helmholtz(alpha*alpha + coriolis*coriolis, -gh, r);

				std::complex<double> inv_k = 1.0/(coriolis*coriolis + alpha*alpha);
				SphereDataComplex a = u0 + inv_r*opComplex.robert_grad_lon(phi);
				SphereDataComplex b = v0 + inv_r*opComplex.robert_grad_lat(phi);
				u = (inv_k*alpha)*a - (inv_k*coriolis)*b;
				v = (inv_k*coriolis)*a + (inv_k*alpha)*b;
			}
			else
			{
				SphereDataComplex Fc_k =
						coriolis*inv_r*(
							-(alpha*alpha*u0 - coriolis*coriolis*opComplex.mu2(u0)) +
							2.0*alpha*coriolis*opComplex.mu(v0)
						);

				SphereDataComplex foo =
						gh*(div0 - coriolis*(1.0/alpha)*opComplex.mu(vort0)) +
						(alpha*phi0 + coriolis*coriolis*(1.0/alpha)*opComplex.mu2(phi0));

				SphereDataComplex rhs =	alpha*alpha*foo +
						coriolis*coriolis*opComplex.mu2(foo)
						- (gh/alpha)*Fc_k;


				phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

				SphereDataComplex a = u0 + inv_r*opComplex.robert_grad_lon(phi);
				SphereDataComplex b = v0 + inv_r*opComplex.robert_grad_lat(phi);

				SphereDataComplex rhsa(sphereDataConfig);
				SphereDataComplex rhsb(sphereDataConfig);

				rhsa = alpha*a - coriolis*opComplex.mu(b);
				rhsb = coriolis*opComplex.mu(a) + alpha*b;

				u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
				v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			}
		}
		else
		{
			SphereDataComplex rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

			u = (1.0/alpha) * (u0 + inv_r*opComplex.robert_grad_lon(phi));
			v = (1.0/alpha) * (v0 + inv_r*opComplex.robert_grad_lat(phi));
		}


		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert_real(phi * beta);
		o_u = Convert_SphereDataComplex_To_SphereData::physical_convert_real(u * beta);
		o_v = Convert_SphereDataComplex_To_SphereData::physical_convert_real(v * beta);
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
			double i_coriolis,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_use_f_sphere,
			int i_variant_id = 0
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereDataConfigSolver = i_sphereDataConfigSolver;

		use_f_sphere = i_use_f_sphere;
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		inv_r = 1.0/r;

		variant_id = i_variant_id;

		if (use_f_sphere)
			coriolis = i_coriolis;
		else
			coriolis = 2.0*i_coriolis;

		gh = i_avg_geopotential;

		op.setup(sphereDataConfig, r);
		opComplex.setup(sphereDataConfig, r);


		if (coriolis != 0)
		{
			if (use_f_sphere)
			{
				// Not needed
			}
			else
			{
				sphSolverPhi.setup(sphereDataConfigSolver, 4);
				sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
				sphSolverPhi.solver_component_rexi_z2(	2.0*coriolis*coriolis*alpha*alpha, r);
				sphSolverPhi.solver_component_rexi_z3(	(coriolis*coriolis)*(coriolis*coriolis), r);
				sphSolverPhi.solver_component_rexi_z4robert(	-gh*alpha*coriolis, r);
				sphSolverPhi.solver_component_rexi_z5robert(	gh/alpha*coriolis*coriolis*coriolis, r);
				sphSolverPhi.solver_component_rexi_z6robert(	gh*2.0*coriolis*coriolis, r);
				sphSolverPhi.solver_component_rexi_z7(	-gh*alpha*alpha, r);
				sphSolverPhi.solver_component_rexi_z8(	-gh*coriolis*coriolis, r);

				mug.setup(sphereDataConfig);
				mug.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, std::complex<double> &o_data)
					{
						o_data = mu;
					}
				);


				sphSolverVel.setup(sphereDataConfigSolver, 2);
				sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
				sphSolverVel.solver_component_rexi_z2(	coriolis*coriolis, r);
			}
		}
		else
		{
			// Not needed
		}
#if 0
		if (coriolis != 0 && !use_f_sphere)
			sphSolverPhi.setup(sphereDataConfigSolver, 2);
		else
			sphSolverPhi.setup(sphereDataConfigSolver, 0);

		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha), r);

		if (coriolis != 0)
		{
			if (use_f_sphere)
				sphSolverPhi.solver_component_rexi_z1(	coriolis*coriolis, r);
			else
				sphSolverPhi.solver_component_rexi_z2(	coriolis*coriolis, r);
		}

		sphSolverPhi.solver_component_rexi_z7(	-gh, r);
#endif
	}


#if 0
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
				o_data = mu*coriolis;
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


		if (coriolis != 0)
		{
			/**
			 * Both versions (this and the version below) results of similar accuracy
			 */
			// only valid for Robert formulation!

#if 1
			SphereDataPhysicalComplex Fc_k = coriolis*inv_r*(
										-(alpha*alpha*u0 - coriolis*coriolis*fg*fg*u0) +
										2.0*alpha*coriolis*fg*v0
									);

			SphereDataPhysicalComplex foo = 	gh*(div0 - coriolis*(1.0/alpha)*fg*eta0) +
										(alpha*phi0 + coriolis*coriolis*(1.0/alpha)*fg*fg*phi0);

			SphereDataPhysicalComplex rhs =	alpha*alpha*foo +
									coriolis*coriolis*fg*fg*foo
									- (gh/alpha)*Fc_k;

#else
			SphereDataComplex Fc_k =	coriolis*inv_r*(
										-(alpha*alpha*(SphereDataComplex)u0 - coriolis*coriolis*opComplex.mu2(u0)) +
										2.0*alpha*coriolis*opComplex.mu(v0)
									);

			SphereDataComplex foo = 	gh*((SphereDataComplex)div0 - coriolis*(1.0/alpha)*opComplex.mu(eta0)) +
										(alpha*(SphereDataComplex)phi0 + coriolis*coriolis*(1.0/alpha)*opComplex.mu2(phi0));

			SphereDataComplex rhs =	alpha*alpha*foo +
									coriolis*coriolis*opComplex.mu2(foo)
									- (gh/alpha)*Fc_k;
#endif

			phi = sphSolverPhi.solve(((SphereDataComplex)rhs).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

#if 1
			SphereDataPhysicalComplex a = u0 + inv_r*(opComplex.robert_grad_lon(phi)).getSphereDataPhysicalComplex();
			SphereDataPhysicalComplex b = v0 + inv_r*(opComplex.robert_grad_lat(phi)).getSphereDataPhysicalComplex();

			SphereDataComplex rhsa = alpha*a - coriolis*fg*b;
			SphereDataComplex rhsb = coriolis*fg*a + alpha*b;

			u = sphSolverVel.solve(((SphereDataComplex)rhsa).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(((SphereDataComplex)rhsb).spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

#else

			SphereDataComplex a = (SphereDataComplex)u0 + inv_r*opComplex.robert_grad_lon(phi);
			SphereDataComplex b = (SphereDataComplex)v0 + inv_r*opComplex.robert_grad_lat(phi);

			SphereDataComplex rhsa = alpha*a - coriolis*opComplex.mu(b);
			SphereDataComplex rhsb = coriolis*opComplex.mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
#endif
		}
		else
		{
			SphereDataComplex rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

			u = (1.0/alpha) * (u0 + inv_r*opComplex.robert_grad_lon(phi).getSphereDataPhysicalComplex());
			v = (1.0/alpha) * (v0 + inv_r*opComplex.robert_grad_lat(phi).getSphereDataPhysicalComplex());
		}

		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert_real(phi * beta);
		o_u = Convert_SphereDataComplex_To_SphereData::physical_convert_real(u * beta);
		o_v = Convert_SphereDataComplex_To_SphereData::physical_convert_real(v * beta);
	}
#endif



	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	inline
	void solve_vectorinvariant_progphivortdiv(
			const SphereData &i_phi0,
			const SphereData &i_vort0,
			const SphereData &i_div0,

			SphereData &o_phi,
			SphereData &o_vort,
			SphereData &o_div
	)
	{
		const SphereDataComplex phi0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_phi0);
		const SphereDataComplex vort0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_vort0);
		const SphereDataComplex div0 = Convert_SphereData_To_SphereDataComplex::physical_convert(i_div0);

		SphereDataComplex phi(sphereDataConfig);
		SphereDataComplex vort(sphereDataConfig);
		SphereDataComplex div(sphereDataConfig);

		if (coriolis != 0)
		{
			if (use_f_sphere)
			{
				SphereDataComplex rhs = gh*(div0 - coriolis/alpha*vort0) + (alpha+coriolis*coriolis/alpha)*phi0;
				phi = rhs.spectral_solve_helmholtz(alpha*alpha + coriolis*coriolis, -gh, r);

				div = -1.0/gh*(phi0 - alpha*phi);
				vort = (1.0/alpha)*(vort0 + coriolis*(div));
			}
			else
			{
				SphereDataPhysicalComplex u0g(sphereDataConfig);
				SphereDataPhysicalComplex v0g(sphereDataConfig);
				opComplex.robert_vortdiv_to_uv(vort0, div0, u0g, v0g, r);

				SphereDataComplex rhs(sphereDataConfig);

				if (variant_id == 0)
				{
					SphereDataPhysicalComplex phi0g = phi0.getSphereDataPhysicalComplex();

					SphereDataPhysicalComplex Fc_k =
							coriolis*inv_r*(
									-(-coriolis*coriolis*mug*mug + alpha*alpha)*u0g
									+ 2.0*alpha*coriolis*mug*v0g
							);

					SphereDataPhysicalComplex foo =
							(gh*(div0.getSphereDataPhysicalComplex() - (1.0/alpha)*coriolis*mug*vort0.getSphereDataPhysicalComplex())) +
							(alpha*phi0g + (1.0/alpha)*coriolis*coriolis*mug*mug*phi0g);

					SphereDataPhysicalComplex rhsg =
							alpha*alpha*foo +
							coriolis*coriolis*mug*mug*foo
							- (gh/alpha)*Fc_k;

					rhs = rhsg;
				}
				else
				{
// see normal_mode_exp_analysis_earth_scale
//#error "DONT USE THIS FORMULATION AS IT CREATES OUTLIERS IN THE DISPERSION!!!"

					SphereDataComplex Fc_k =
							coriolis*inv_r*(
									-(-coriolis*coriolis*mug*mug + alpha*alpha)*u0g
									+ 2.0*alpha*coriolis*mug*v0g
							);

					SphereDataComplex foo =
							(gh*(div0 - (1.0/alpha)*coriolis*opComplex.mu(vort0))) +
							(alpha*phi0 + (1.0/alpha)*coriolis*coriolis*opComplex.mu2(phi0));

					rhs =	alpha*alpha*foo +
							coriolis*coriolis*opComplex.mu2(foo)
							- (gh/alpha)*Fc_k;
				}

				phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

				/*
				 * Solve without inverting a matrix
				 */
				SphereDataPhysicalComplex u0(sphereDataConfig);
				SphereDataPhysicalComplex v0(sphereDataConfig);

				opComplex.robert_vortdiv_to_uv(vort0, div0, u0, v0, r);


				SphereDataPhysicalComplex a(sphereDataConfig);
				SphereDataPhysicalComplex b(sphereDataConfig);
				if (variant_id == 0)
				{
					SphereDataPhysicalComplex gradu(sphereDataConfig);
					SphereDataPhysicalComplex gradv(sphereDataConfig);
					opComplex.robert_grad_to_vec(phi, gradu, gradv, r);

					a = u0 + gradu;
					b = v0 + gradv;
				}
				else
				{
					// see normal_mode_exp_analysis_earth_scale
					//	#error "DONT USE THIS FORMULATION AS IT CREATES OUTLIERS IN DISPERSION!!!"

					a = u0 + inv_r*opComplex.robert_grad_lon(phi).getSphereDataPhysicalComplex();
					b = v0 + inv_r*opComplex.robert_grad_lat(phi).getSphereDataPhysicalComplex();
				}

				SphereDataPhysicalComplex k = (coriolis*coriolis*(mug*mug)+alpha*alpha);
				SphereDataPhysicalComplex u = (alpha*a - coriolis*mug*(b))/k;
				SphereDataPhysicalComplex v = (coriolis*mug*(a) + alpha*b)/k;

				opComplex.robert_uv_to_vortdiv(u, v, vort, div, r);
			}
		}
		else
		{
			SphereDataComplex rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

			div = -1.0/gh*(phi0 - alpha*phi);
			vort = 1.0/alpha*(vort0);
		}


		o_phi = Convert_SphereDataComplex_To_SphereData::physical_convert_real(phi * beta);
		o_vort = Convert_SphereDataComplex_To_SphereData::physical_convert_real(vort * beta);
		o_div = Convert_SphereDataComplex_To_SphereData::physical_convert_real(div * beta);
	}


};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
