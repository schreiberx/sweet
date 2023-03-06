/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWEREXI_SPH_NEW_HPP_
#define SRC_SWEREXI_SPH_NEW_HPP_

#include <complex>
#include <sweet/core/sphere/Convert_SphereDataSpectral_to_SphereDataSpectralComplex.hpp>
#include <sweet/core/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_PhysicalComplex.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperatorsComplex.hpp>
#include "../helpers/SWESphBandedMatrixPhysicalComplex.hpp"



/**
 * REXI solver for SWE
 */
class SWERexiTerm_SPH
{
	/// SPH configuration
	const sweet::SphereDataConfig *sphereDataConfig;

	/// Solver for given alpha
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolverComplexDiv;

	sweet::SphereOperatorsComplex opsComplex;

	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;


	/// (real) timestep size
	double timestep_size;

	/// earth radius
	double sphere_radius;

	/// inverse of earth radius
	double ir;

	/// Coriolis omega
	double two_coriolis_omega;

	/// f0
	double f0;

	bool use_f_sphere;

	bool no_coriolis;

	/// Average geopotential
	double gh0;

	/// pseudo timestep size to use implicit solver for REXI
	std::complex<double> dt_implicit;

public:
	SWERexiTerm_SPH()	:
		sphereDataConfig(nullptr)
	{
	}




	/**
	 * REXI term formulated as backward Euler
	 */
	inline
	void solve_vectorinvariant_progphivortdiv(
			const sweet::SphereData_Spectral &i_phi0,
			const sweet::SphereData_Spectral &i_vrt0,
			const sweet::SphereData_Spectral &i_div0,

			sweet::SphereData_Spectral &o_phi,
			sweet::SphereData_Spectral &o_vort,
			sweet::SphereData_Spectral &o_div
	)
	{
		sweet::SphereData_SpectralComplex phi1(sphereDataConfig);
		sweet::SphereData_SpectralComplex vrt1(sphereDataConfig);
		sweet::SphereData_SpectralComplex div1(sphereDataConfig);

		sweet::SphereData_SpectralComplex phi0 = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_phi0);
		sweet::SphereData_SpectralComplex vrt0 = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_vrt0);
		sweet::SphereData_SpectralComplex div0 = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_div0);

		/*
		 * Preprocessing:
		 * U0* = U0 * beta/alpha
		 * dt_implicit = -timestep_size/alpha
		 */
		phi0 *= beta/alpha;
		div0 *= beta/alpha;
		vrt0 *= beta/alpha;

		if (no_coriolis)
		{
#if 0

			sweet::SphereData_SpectralComplex rhs = div0 + opsComplex.implicit_L(phi0, dt_implicit);
			div1 = opsComplex.implicit_helmholtz(rhs, gh0*dt_implicit*dt_implicit, sphere_radius);

			phi1 = phi0 - dt_implicit*gh0*div1;
			vrt1 = vrt0;

#else

			std::complex<double> dt_two_omega = dt_implicit*0.0;

			sweet::SphereData_SpectralComplex rhs = div0 + opsComplex.implicit_FJinv(vrt0, dt_two_omega) + opsComplex.implicit_L(phi0, dt_implicit);
			div1 = sphSolverComplexDiv.solve(rhs);
			phi1 = phi0 - dt_implicit*gh0*div1;
			vrt1 = opsComplex.implicit_Jinv(vrt0 - opsComplex.implicit_F(div1, dt_two_omega), dt_two_omega);

#endif
		}
		else if (use_f_sphere)
		{
			SWEETError("TODO: This needs to be ported to the new formulation");
#if 0
			// TODO
			sweet::SphereData_SpectralComplex rhs = dt_implicit*gh0*(div0 - f0*vrt0) + (1.0+f0*f0)*phi0;
			phi1 = rhs.spectral_solve_helmholtz_one(-gh0, sphere_radius);

			vrt1 = (vrt0 + f0*div1);
			div1 = -1.0/(dt_implicit*gh0)*(phi0 - phi1);
#endif
		}
		else
		{
			std::complex<double> dt_two_omega = dt_implicit*two_coriolis_omega;

			sweet::SphereData_SpectralComplex rhs = div0 + opsComplex.implicit_FJinv(vrt0, dt_two_omega) + opsComplex.implicit_L(phi0, dt_implicit);
			div1 = sphSolverComplexDiv.solve(rhs);
			phi1 = phi0 - dt_implicit*gh0*div1;
			vrt1 = opsComplex.implicit_Jinv(vrt0 - opsComplex.implicit_F(div1, dt_two_omega), dt_two_omega);
		}

		o_phi = Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(phi1);
		o_vort = Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(vrt1);
		o_div = Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(div1);
	}



	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup_vectorinvariant_progphivortdiv(
			const sweet::SphereDataConfig *i_sphereDataConfigSolver,
			const sweet::ShackDictionary *i_shackDict,

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
		sphereDataConfig = i_sphereDataConfigSolver;

		sphere_radius = i_radius;
		gh0 = i_avg_geopotential;
		no_coriolis = i_no_coriolis;
		timestep_size = i_timestep_size;
		use_f_sphere = i_use_f_sphere;

		if (use_f_sphere)
			f0 = i_f0;
		else
			two_coriolis_omega = 2.0*i_coriolis_omega;

		alpha = i_alpha;
		beta = i_beta;

		opsComplex.setup(sphereDataConfig, &(i_shackDict->sim));

		dt_implicit = -timestep_size/alpha;

		if (!no_coriolis)
		{
			if (!use_f_sphere)
			{
				std::complex<double> dt_two_omega = dt_implicit*two_coriolis_omega;

				sphSolverComplexDiv.setup(sphereDataConfig, 4);
				sphSolverComplexDiv.solver_component_implicit_J(dt_two_omega);
				sphSolverComplexDiv.solver_component_implicit_FJinvF(dt_two_omega);
				sphSolverComplexDiv.solver_component_implicit_L(gh0*dt_implicit, dt_implicit, sphere_radius);
			}
		}
		else
		{
			std::complex<double> dt_two_omega = dt_implicit*0.0;

			sphSolverComplexDiv.setup(sphereDataConfig, 4);
			sphSolverComplexDiv.solver_component_implicit_J(dt_two_omega);
			sphSolverComplexDiv.solver_component_implicit_FJinvF(dt_two_omega);
			sphSolverComplexDiv.solver_component_implicit_L(gh0*dt_implicit, dt_implicit, sphere_radius);

		}
	}


};


#endif
