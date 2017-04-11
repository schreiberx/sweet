/*
 * SWEImplicit_SPHRobert.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
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

	/// Crank-Nicolson damping factor
	double crank_nicolson_damping_factor = 0.5;

	/// false if formulation with coriolis effect should NOT be used even if f != 0
//	bool use_formulation_with_coriolis_effect;

	bool use_f_sphere;

	// order of time stepping.
	// 1: backward Euler
	// 2: Crank Nicolson
	int timestepping_order;

	// PDE variant
	int pde_variant_id;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// set to f0 for f-sphere and 2*coriolis_omega for non-f-sphere
	double coriolis;

	/// Average geopotential
	double gh;

//	std::complex<double> I = 1;

	SphereDataPhysical fg;
	SphereDataPhysical mug;

public:
	SWEImplicit_SPHRobert(SphereOperators &i_op)	:
		op(i_op),
		sphereDataConfig(nullptr),
		sphereDataConfigSolver(nullptr)
	{
	}


	inline
	void setCrankNicolsonDampingFactor(
			double i_crank_nicolson_damping_factor
	)
	{
		crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;
	}



	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup_velocityformulation_progphiuv(
			SphereDataConfig *i_sphereDataConfig,
			SphereDataConfig *i_sphereDataConfigSolver,

			double i_radius,
			double i_coriolis_omega,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_f_sphere,
			int i_timestepping_order
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereDataConfigSolver = i_sphereDataConfigSolver;
		timestep_size = i_timestep_size;
		use_f_sphere = i_f_sphere;
		timestepping_order = i_timestepping_order;

		if (use_f_sphere)
			coriolis = i_coriolis_omega;
		else
			coriolis = 2.0*i_coriolis_omega;

		alpha = -1.0/timestep_size;
		beta = -1.0/timestep_size;

		if (i_timestepping_order == 2)
		{
			/*
			 * Crank-Nicolson method:
			 *
			 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
			 *
			 * with q the CN damping facor with no damping for q=0.5
			 */

			// scale dt by the damping factor to reuse solver structure
			alpha /= crank_nicolson_damping_factor;
			beta /= crank_nicolson_damping_factor;
		}

		r = i_radius;
		inv_r = 1.0/r;

		gh = i_avg_geopotential;

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
	void solve_velocityformulation_progphiuv(
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

		SphereData phi0 = i_phi0;
		SphereData u0 = i_u0;
		SphereData v0 = i_v0;


		/*
		 * Crank-Nicolson method:
		 *
		 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
		 *
		 * with q the CN damping facor with no damping for q=0.5
		 */

		if (timestepping_order == 2)
		{
			/**
			 * Compute explicit time step
			 */
			SphereData &o_phi_t = o_phi;
			SphereData &o_u_t = o_u;
			SphereData &o_v_t = o_v;

			if (coriolis != 0)
			{
				if (use_f_sphere)
				{
					o_phi_t = -(op.robert_div_lon(i_u0)+op.robert_div_lat(i_v0)) * (gh*inv_r);
					o_u_t = -op.robert_grad_lon(i_phi0) * inv_r + coriolis*(i_v0);
					o_v_t = -op.robert_grad_lat(i_phi0) * inv_r - coriolis*(i_u0);
				}
				else
				{
					o_phi_t = -(op.robert_div_lon(i_u0)+op.robert_div_lat(i_v0)) * (gh*inv_r);
					o_u_t = -op.robert_grad_lon(i_phi0) * inv_r + f(i_v0);
					o_v_t = -op.robert_grad_lat(i_phi0) * inv_r - f(i_u0);
				}
			}
			else
			{
				o_phi_t = -(op.robert_div_lon(i_u0)+op.robert_div_lat(i_v0)) * (gh*inv_r);
				o_u_t = -op.robert_grad_lon(i_phi0) * inv_r;
				o_v_t = -op.robert_grad_lat(i_phi0) * inv_r;
			}

			double fac = timestep_size*(1.0-crank_nicolson_damping_factor);
			// run single time step for rhs

			phi0 += fac*o_phi_t;
			u0 += fac*o_u_t;
			v0 += fac*o_v_t;
		}


		SphereData div0 = inv_r*op.robert_div(u0, v0);
		SphereData vort0 = inv_r*op.robert_vort(u0, v0);

		SphereData phi(sphereDataConfig);
		SphereData u(sphereDataConfig);
		SphereData v(sphereDataConfig);

		if (coriolis != 0)
		{
			if (use_f_sphere)
			{
				SphereData rhs = gh*(div0 - coriolis/alpha*vort0) + (alpha+coriolis*coriolis/alpha)*phi0;
				phi = rhs.spectral_solve_helmholtz(alpha*alpha + coriolis*coriolis, -gh, r);

				double inv_k = 1.0/(coriolis*coriolis + alpha*alpha);
				SphereData a = u0 + inv_r*op.robert_grad_lon(phi);
				SphereData b = v0 + inv_r*op.robert_grad_lat(phi);
				u = (inv_k*alpha)*a - (inv_k*coriolis)*b;
				v = (inv_k*coriolis)*a + (inv_k*alpha)*b;
			}
			else
			{
				SphereData Fc_k =
						coriolis*inv_r*(
							-(alpha*alpha*u0 - coriolis*coriolis*op.mu2(u0)) +
							2.0*alpha*coriolis*op.mu(v0)
						);

				SphereData foo =
						gh*(div0 - coriolis*(1.0/alpha)*op.mu(vort0)) +
						(alpha*phi0 + coriolis*coriolis*(1.0/alpha)*op.mu2(phi0));

				SphereData rhs =	alpha*alpha*foo +
						coriolis*coriolis*op.mu2(foo)
						- (gh/alpha)*Fc_k;


				phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);

				SphereData a = u0 + inv_r*op.robert_grad_lon(phi);
				SphereData b = v0 + inv_r*op.robert_grad_lat(phi);

				SphereData rhsa(sphereDataConfig);
				SphereData rhsb(sphereDataConfig);

				rhsa = alpha*a - coriolis*op.mu(b);
				rhsb = coriolis*op.mu(a) + alpha*b;

				u = sphSolverVel.solve(rhsa.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
				v = sphSolverVel.solve(rhsb.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
			}
		}
		else
		{
			SphereData rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

			u = (1.0/alpha) * (u0 + inv_r*op.robert_grad_lon(phi));
			v = (1.0/alpha) * (v0 + inv_r*op.robert_grad_lat(phi));
		}

		o_phi = phi * beta;
		o_u = u * beta;
		o_v = v * beta;
	}



	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup_vectorinvariantformulation_progphivortdiv(
			SphereDataConfig *i_sphereDataConfig,
			SphereDataConfig *i_sphereDataConfigSolver,

			double i_r,
			double i_coriolis_omega,
			double i_avg_geopotential,
			double i_timestep_size,

			bool i_f_sphere,
			int i_timestepping_order,

			int i_pde_variant_id
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereDataConfigSolver = i_sphereDataConfigSolver;
		timestep_size = i_timestep_size;
		use_f_sphere = i_f_sphere;
		timestepping_order = i_timestepping_order;

		pde_variant_id = i_pde_variant_id;

		if (use_f_sphere)
			coriolis = i_coriolis_omega;
		else
			coriolis = 2.0*i_coriolis_omega;

		alpha = -1.0/timestep_size;
		beta = -1.0/timestep_size;

		if (i_timestepping_order == 2)
		{
			/*
			 * Crank-Nicolson method:
			 *
			 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
			 *
			 * with q the CN damping facor with no damping for q=0.5
			 */

			// scale dt by the damping factor to reuse solver structure

			alpha /= crank_nicolson_damping_factor;
			beta /= crank_nicolson_damping_factor;
		}

		r = i_r;
		inv_r = 1.0/r;

		gh = i_avg_geopotential;

		op.setup(sphereDataConfig, r);

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

				fg.setup(sphereDataConfig);
				fg.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, double &o_data)
					{
						o_data = mu*coriolis;
					}
				);

				mug.setup(sphereDataConfig);
				mug.physical_update_lambda_gaussian_grid(
					[&](double lon, double mu, double &o_data)
					{
						o_data = mu;
					}
				);

#if 0
				sphSolverVel.setup(sphereDataConfigSolver, 2);
				sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
				sphSolverVel.solver_component_rexi_z2(	coriolis*coriolis, r);
#endif
			}
		}
		else
		{
			// Not needed
		}
	}



	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	void solve_vectorinvariantformulation_progphivortdiv(
			const SphereData &i_phi0,
			const SphereData &i_vort0,
			const SphereData &i_div0,

			SphereData &o_phi,
			SphereData &o_vort,
			SphereData &o_div,

			double &o_dt
	)
	{
		o_dt = timestep_size;

		SphereData phi0 = i_phi0;
		SphereData vort0 = i_vort0;
		SphereData div0 = i_div0;

		/*
		 * Crank-Nicolson method:
		 *
		 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
		 *
		 * with q the CN damping facor with no damping for q=0.5
		 */

		if (timestepping_order == 2)
		{
			SphereData o_phi_t(sphereDataConfig);
			SphereData o_vort_t(sphereDataConfig);
			SphereData o_div_t(sphereDataConfig);

			/*
			 * LINEAR
			 */
			if (coriolis != 0)
			{
				if (use_f_sphere)
				{
					o_phi_t = -gh*div0;
					o_div_t = -op.laplace(phi0) + coriolis*vort0;
					o_vort_t = -coriolis*div0;
				}
				else
				{
					SphereDataPhysical ug(sphereDataConfig);
					SphereDataPhysical vg(sphereDataConfig);

					op.robert_vortdiv_to_uv(vort0, div0, ug, vg);
					SphereDataPhysical phig = phi0.getSphereDataPhysical();

					SphereDataPhysical tmpg1 = ug*fg;
					SphereDataPhysical tmpg2 = vg*fg;

					op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

					o_vort_t *= -1.0;

					tmpg1 = ug*gh;
					tmpg2 = vg*gh;

					SphereData tmpspec(sphereDataConfig);
					op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

					o_phi_t *= -1.0;

					tmpspec = phig;
					tmpspec.request_data_spectral();
					o_div_t += -op.laplace(tmpspec);
				}
			}
			else
			{
				o_phi_t = -gh*div0;
				o_div_t = -op.laplace(phi0);
				o_vort_t.spectral_set_zero();
			}

			double fac = timestep_size*(1.0-crank_nicolson_damping_factor);
			// run single time step for rhs

			phi0 += fac*o_phi_t;
			vort0 += fac*o_vort_t;
			div0 += fac*o_div_t;
		}


		SphereData phi(sphereDataConfig);
		SphereData vort(sphereDataConfig);
		SphereData div(sphereDataConfig);

		if (coriolis != 0)
		{
			if (use_f_sphere)
			{
				SphereData rhs = gh*(div0 - coriolis/alpha*vort0) + (alpha+coriolis*coriolis/alpha)*phi0;
				phi = rhs.spectral_solve_helmholtz(alpha*alpha + coriolis*coriolis, -gh, r);

				div = -1.0/gh*(phi0 - alpha*phi);
				vort = (1.0/alpha)*(vort0 + coriolis*(div));
			
#if 0
				double inv_k = 1.0/(coriolis*coriolis + alpha*alpha);
				SphereData a = u0 + inv_r*op.robert_grad_lon(phi);
				SphereData b = v0 + inv_r*op.robert_grad_lat(phi);
				u = (inv_k*alpha)*a - (inv_k*coriolis)*b;
				v = (inv_k*coriolis)*a + (inv_k*alpha)*b;
#endif
			}
			else
			{
				SphereData rhs(sphereDataConfig);

				SphereDataPhysical u0g(sphereDataConfig);
				SphereDataPhysical v0g(sphereDataConfig);
				op.robert_vortdiv_to_uv(vort0, div0, u0g, v0g);

				if (pde_variant_id == 0)
				{
					SphereDataPhysical phi0g = phi0.getSphereDataPhysical();

					SphereDataPhysical Fc_k =
							coriolis*inv_r*(
									-(-coriolis*coriolis*mug*mug + alpha*alpha)*u0g
									+ 2.0*alpha*coriolis*mug*v0g
							);

					SphereDataPhysical foo =
							(gh*(div0.getSphereDataPhysical() - (1.0/alpha)*coriolis*mug*vort0.getSphereDataPhysical())) +
							(alpha*phi0g + (1.0/alpha)*coriolis*coriolis*mug*mug*phi0g);

					SphereDataPhysical rhsg =
							alpha*alpha*foo +
							coriolis*coriolis*mug*mug*foo
							- (gh/alpha)*Fc_k;

					rhs = rhsg;
				}
				else
				{
					SphereData Fc_k =
							coriolis*inv_r*(
									-(-coriolis*coriolis*mug*mug + alpha*alpha)*u0g
									+ 2.0*alpha*coriolis*mug*v0g
							);

					SphereData foo =
							gh*(div0 - coriolis*(1.0/alpha)*op.mu(vort0)) +
							(alpha*phi0 + coriolis*coriolis*(1.0/alpha)*op.mu2(phi0));

					rhs =	alpha*alpha*foo +
							coriolis*coriolis*op.mu2(foo)
							- (gh/alpha)*Fc_k;
				}

				phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);


				SphereDataPhysical u0(sphereDataConfig);
				SphereDataPhysical v0(sphereDataConfig);

				op.robert_vortdiv_to_uv(vort0, div0, u0, v0);

				SphereDataPhysical a(sphereDataConfig);
				SphereDataPhysical b(sphereDataConfig);

				if (pde_variant_id == 0)
				{
					SphereDataPhysical gradu(sphereDataConfig);
					SphereDataPhysical gradv(sphereDataConfig);
					op.robert_grad_to_vec(phi, gradu, gradv, r);

					a = u0 + gradu;
					b = v0 + gradv;
				}
				else
				{

					a = u0 + inv_r*op.robert_grad_lon(phi).getSphereDataPhysical();
					b = v0 + inv_r*op.robert_grad_lat(phi).getSphereDataPhysical();
				}


				SphereDataPhysical k = (coriolis*coriolis*mug*mug+alpha*alpha);
				SphereDataPhysical u = (alpha*a - coriolis*mug*(b))/k;
				SphereDataPhysical v = (coriolis*mug*(a) + alpha*b)/k;

				op.robert_uv_to_vortdiv(u, v, vort, div);
			}
		}
		else
		{
			SphereData rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);

			vort = 1.0/alpha*(vort0);
			div = -1.0/gh*(phi0 - alpha*phi);
		}

		o_phi = (phi * beta);
		o_vort = (vort * beta);
		o_div = (div * beta);
	}


	inline
	SphereData f(const SphereData &i_sphData)
	{
		return op.mu(i_sphData*coriolis);
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
		SphereData vort0 = inv_r*op.robert_vort(u0, v0);

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

			SphereData foo = 	avg_geopotential*(div0 - two_omega*(1.0/alpha)*op.mu(vort0)) +
										(alpha*phi0 + two_omega*two_omega*(1.0/alpha)*op.mu2(phi0));

			SphereData rhs =	alpha*alpha*foo +
									two_omega*two_omega*op.mu2(foo)
									- (avg_geopotential/alpha)*Fc_k;

			phi = sphSolverPhi.solve(rhs.spectral_returnWithDifferentModes(sphereDataConfigSolver)).spectral_returnWithDifferentModes(sphereDataConfig);
		}
		else
*/
		{
			SphereData rhs = gh*div0 + alpha*phi0;
			phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);
		}

		o_phi = phi * beta;
	}
};


#endif /* SRC_SWEREXI_SPHROBERT_HPP_ */
