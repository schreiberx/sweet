
#ifndef TIMESTEPPING_IMPLICIT_RK1_HPP
#define TIMESTEPPING_IMPLICIT_RK1_HPP

#include <sweet/sphere/SphereData.hpp>
#include <limits>

class SWETimesteppingImplicitRK1
{
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
	SWETimesteppingImplicitRK1()	:
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

		r = i_radius;
		inv_r = 1.0/r;

		coriolis_omega = i_coriolis_omega;

		two_omega = 2.0*coriolis_omega;
		avg_geopotential = i_avg_geopotential;

		sphSolverPhi.setup(sphereDataConfigSolver, 4);
		sphSolverPhi.solver_component_rexi_z1(	I, r);

		if (use_formulation_with_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega, r);
			sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);

			sphSolverPhi.solver_component_rexi_z4robert(	-avg_geopotential*two_omega, r);
			sphSolverPhi.solver_component_rexi_z5robert(	avg_geopotential*two_omega*two_omega*two_omega, r);
			sphSolverPhi.solver_component_rexi_z6robert(	avg_geopotential*2.0*two_omega*two_omega, r);
		}

		sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential, r);
		if (use_formulation_with_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z8(	-avg_geopotential*two_omega*two_omega, r);
		}



		sphSolverVel.setup(sphereDataConfigSolver, 2);
		sphSolverVel.solver_component_rexi_z1(	I, r);
		if (use_formulation_with_coriolis_effect)
		{
			sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);
		}
	}




	SphereDataComplex kappa(
			const SphereDataComplex &i_data
	)	const
	{
		return i_data + two_omega*two_omega*SphereOperatorsComplex::mu2(i_data);
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
										-(u0 - two_omega*two_omega*SphereOperatorsComplex::mu2(u0)) +
										2.0*two_omega*SphereOperatorsComplex::mu(v0)
									);

			SphereDataComplex foo = 	avg_geopotential*(div0 - two_omega*SphereOperatorsComplex::mu(eta0)) +
										(alpha*phi0 + two_omega*two_omega*SphereOperatorsComplex::mu2(phi0));

			SphereDataComplex rhs =	foo +
									two_omega*two_omega*SphereOperatorsComplex::mu2(foo)
									- avg_geopotential*Fc_k;

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

			SphereDataComplex Fp_i = fj*(f*f);
			SphereDataComplex Fp_j = fj*(2.0*f);

			SphereDataComplex Fck = Fp_i*u0 + Fp_j*v0;

			SphereDataComplex rhs =
					kappa(
							phi_bar*(div0 - f*eta0)
							+ (I + f*f)*phi0
					)
					- phi_bar/alpha*Fck;

#endif

			phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

			SphereDataComplex a = u0 + inv_r*SphereOperatorsComplex::robert_grad_lon(phi);
			SphereDataComplex b = v0 + inv_r*SphereOperatorsComplex::robert_grad_lat(phi);

			SphereDataComplex rhsa = a - two_omega*SphereOperatorsComplex::mu(b);
			SphereDataComplex rhsb = two_omega*SphereOperatorsComplex::mu(a) + b;

			u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
		}
		else
		{
			SphereDataComplex rhs = avg_geopotential*div0 + phi0;
			phi = rhs.spectral_solve_helmholtz(I, -avg_geopotential, r);

			u = (u0 + inv_r*SphereOperatorsComplex::robert_grad_lon(phi));
			v = (v0 + inv_r*SphereOperatorsComplex::robert_grad_lat(phi));
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

#endif
