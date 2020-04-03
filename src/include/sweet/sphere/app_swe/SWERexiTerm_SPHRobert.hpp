/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_SWEREXI_SPHROBERT_HPP_
#define SRC_SWEREXI_SPHROBERT_HPP_

#include <complex>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalComplex.hpp>
#include <sweet/sphere/Convert_SphereDataSpectral_to_SphereDataSpectralComplex.hpp>
#include <sweet/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_PhysicalComplex.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereOperators_SphereDataComplex.hpp>



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWERexiTerm_SPHRobert
{
	/// SPH configuration
//	SphereDataConfig *sphereDataConfig;

	/// SPH configuration
	const SphereData_Config *sphereDataConfigSolver;

	/// Solver for given alpha
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverPhi;
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverVel;

	SphereOperators_SphereData op;
	SphereOperators_SphereDataComplex opComplex;

	/// scalar in front of RHS
	std::complex<double> rhs_scalar;

	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;


	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// Coriolis omega
	double two_coriolis_omega;

	/// f0
	double f0;

	bool use_f_sphere;

	bool no_coriolis;

	/// Average geopotential
	double gh;

	SphereData_PhysicalComplex mug;

public:
	SWERexiTerm_SPHRobert()	:
//		sphereDataConfig(nullptr),
		sphereDataConfigSolver(nullptr)
	{
	}




	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup_vectorinvariant_progphivortdiv(
			const SphereData_Config *i_sphereDataConfigSolver,
			const SimulationVariables *i_simVars,

			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,

			double i_radius,
			double i_coriolis_omega,
			double i_f0,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_use_f_sphere,
			bool i_no_coriolis
	)
	{
		sphereDataConfigSolver = i_sphereDataConfigSolver;

		use_f_sphere = i_use_f_sphere;
		no_coriolis = i_no_coriolis;

		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		inv_r = 1.0/r;

		if (use_f_sphere)
			f0 = i_f0;
		else
			two_coriolis_omega = 2.0*i_coriolis_omega;

		gh = i_avg_geopotential;

		op.setup(sphereDataConfigSolver, &(i_simVars->sim));
		opComplex.setup(sphereDataConfigSolver, &(i_simVars->sim));


		if (!use_f_sphere)
		{
			sphSolverPhi.setup(sphereDataConfigSolver, 4);
			sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
			sphSolverPhi.solver_component_rexi_z2(	2.0*two_coriolis_omega*two_coriolis_omega*alpha*alpha, r);
			sphSolverPhi.solver_component_rexi_z3(	(two_coriolis_omega*two_coriolis_omega)*(two_coriolis_omega*two_coriolis_omega), r);
			sphSolverPhi.solver_component_rexi_z4robert(	-gh*alpha*two_coriolis_omega, r);
			sphSolverPhi.solver_component_rexi_z5robert(	gh/alpha*two_coriolis_omega*two_coriolis_omega*two_coriolis_omega, r);
			sphSolverPhi.solver_component_rexi_z6robert(	gh*2.0*two_coriolis_omega*two_coriolis_omega, r);
			sphSolverPhi.solver_component_rexi_z7(	-gh*alpha*alpha, r);
			sphSolverPhi.solver_component_rexi_z8(	-gh*two_coriolis_omega*two_coriolis_omega, r);

			mug.setup(sphereDataConfigSolver);
			mug.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, std::complex<double> &o_data)
				{
					o_data = mu;
				}
			);
		}
	}



	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	inline
	void solve_vectorinvariant_progphivortdiv(
			const SphereData_Spectral &i_phi0,
			const SphereData_Spectral &i_vort0,
			const SphereData_Spectral &i_div0,

			SphereData_Spectral &o_phi,
			SphereData_Spectral &o_vort,
			SphereData_Spectral &o_div
	)
	{
		const SphereData_SpectralComplex phi0 = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_phi0);
		const SphereData_SpectralComplex vort0 = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_vort0);
		const SphereData_SpectralComplex div0 = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_div0);

		SphereData_SpectralComplex phi(sphereDataConfigSolver);
		SphereData_SpectralComplex vort(sphereDataConfigSolver);
		SphereData_SpectralComplex div(sphereDataConfigSolver);

		if (no_coriolis)
		{
			SphereData_SpectralComplex rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

			vort = (1.0/alpha)*vort0;
			div = -1.0/gh*(phi0 - alpha*phi);
		}
		else if (use_f_sphere)
		{
			SphereData_SpectralComplex rhs = gh*(div0 - f0/alpha*vort0) + (alpha+f0*f0/alpha)*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha + f0*f0, -gh, r);

			vort = (1.0/alpha)*(vort0 + f0*(div));
			div = -1.0/gh*(phi0 - alpha*phi);
		}
		else
		{
			SphereData_PhysicalComplex u0g(sphereDataConfigSolver);
			SphereData_PhysicalComplex v0g(sphereDataConfigSolver);

			opComplex.robert_vortdiv_to_uv(vort0, div0, u0g, v0g, r);

			SphereData_SpectralComplex rhs(sphereDataConfigSolver);


			SphereData_PhysicalComplex phi0g = phi0.getSphereDataPhysicalComplex();

			SphereData_PhysicalComplex Fc_k =
					two_coriolis_omega*inv_r*(
							-(-two_coriolis_omega*two_coriolis_omega*mug*mug + alpha*alpha)*u0g
							+ 2.0*alpha*two_coriolis_omega*mug*v0g
					);

			SphereData_PhysicalComplex foo =
					(gh*(div0.getSphereDataPhysicalComplex() - (1.0/alpha)*two_coriolis_omega*mug*vort0.getSphereDataPhysicalComplex())) +
					(alpha*phi0g + (1.0/alpha)*two_coriolis_omega*two_coriolis_omega*mug*mug*phi0g);

			SphereData_PhysicalComplex rhsg =
					alpha*alpha*foo +
					two_coriolis_omega*two_coriolis_omega*mug*mug*foo
					- (gh/alpha)*Fc_k;
			// convert to spectral space
			rhs = rhsg;


			phi = sphSolverPhi.solve(rhs);

			/*
			 * Solve without inverting a matrix
			 */
			SphereData_PhysicalComplex u0(sphereDataConfigSolver);
			SphereData_PhysicalComplex v0(sphereDataConfigSolver);

			opComplex.robert_vortdiv_to_uv(vort0, div0, u0, v0, r);

			SphereData_PhysicalComplex a(sphereDataConfigSolver);
			SphereData_PhysicalComplex b(sphereDataConfigSolver);

#if 1

			SphereData_PhysicalComplex gradu(sphereDataConfigSolver);
			SphereData_PhysicalComplex gradv(sphereDataConfigSolver);

			opComplex.robert_grad_to_vec(phi, gradu, gradv, r);
			a = u0 + gradu;
			b = v0 + gradv;

#else
			// LEAVE THIS CODE HERE!!!
			// THIS code resulted in the requirement of a mode extension of 2 additional modes
			// for the REXI solver!!!
			// However, it turned out not to be necessary :-(
            a = u0 + inv_r*opComplex.robert_grad_lon(phi).getSphereDataPhysicalComplex();
            b = v0 + inv_r*opComplex.robert_grad_lat(phi).getSphereDataPhysicalComplex();

#endif

			SphereData_PhysicalComplex k = (two_coriolis_omega*two_coriolis_omega*(mug*mug)+alpha*alpha);
			SphereData_PhysicalComplex u = (alpha*a - two_coriolis_omega*mug*(b))/k;
			SphereData_PhysicalComplex v = (two_coriolis_omega*mug*(a) + alpha*b)/k;

			opComplex.robert_uv_to_vortdiv(u, v, vort, div, r);

		}

		o_phi = Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(phi * beta);
		o_vort = Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(vort * beta);
		o_div = Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(div * beta);
	}


};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
