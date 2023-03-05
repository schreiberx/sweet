/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Include these files and directory for compilation
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators/SWE_Sphere_TS_l_erk.cpp
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators/SWE_Sphere_TS_lg_erk.cpp
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators/SWE_Sphere_TS_l_irk.cpp
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 *
 * MULE_SCONS_OPTIONS: --lapack=enable
 * MULE_SCONS_OPTIONS: --fortran-source=enable
 */

#include <cmath>

#include "../programs/swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"

#include <sweet/core/SimulationVariables.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>

#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperatorsComplex.hpp>

#include <sweet/core/sphere/Convert_SphereDataSpectral_to_SphereDataSpectralComplex.hpp>
#include <sweet/core/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>


#include "../programs/swe_sphere_timeintegrators/helpers/SWESphBandedMatrixPhysicalComplex.hpp"
#include "../programs/swe_sphere_timeintegrators/helpers/SWESphBandedMatrixPhysicalReal.hpp"
#include "../programs/swe_sphere_timeintegrators/SWE_Sphere_TS_l_erk.hpp"
#include "../programs/swe_sphere_timeintegrators/SWE_Sphere_TS_l_irk.hpp"


SimulationVariables simVars;


template<
	typename phys_value_type = double,
	typename sphere_data_spec_type = SphereData_Spectral,
	typename sphere_data_phys_type = SphereData_Physical,
	typename sphere_operators_type = SphereOperators,
	typename sph_banded_solver_type = SphBandedMatrixPhysicalReal<std::complex<double>>
>
class Test
{
	double double_precision_digits = -1.0;

	void check_error(double err, double rel_value)
	{
		if (err > double_precision_digits*rel_value || std::isnan(err))
		{
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
			std::cout << "+ err: " << err << std::endl;
			std::cout << "+ max_err_threshold: " << double_precision_digits*rel_value << std::endl;
			std::cout << "+ rel_value: " << rel_value << std::endl;
			SWEETError("Error too high");
		}
	};


	static
	void benchmark_setup_geostrophic_balance(
			SphereOperators &ops,
			SphereData_Spectral &o_phi,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div
	)
	{
		BenchmarksSphereSWE benchmarks;
		benchmarks.setup(simVars, ops, "geostrophic_balance_linear_16");

		benchmarks.master->get_initial_state(o_phi, o_vrt, o_div);
	}


	static
	void benchmark_setup_geostrophic_balance(
			SphereOperators &ops,
			sweet::SphereData_SpectralComplex &o_phi,
			sweet::SphereData_SpectralComplex &o_vrt,
			sweet::SphereData_SpectralComplex &o_div
	)
	{

		SphereData_Spectral phi(ops.sphereDataConfig);
		SphereData_Spectral vrt(ops.sphereDataConfig);
		SphereData_Spectral div(ops.sphereDataConfig);

		benchmark_setup_geostrophic_balance(ops, phi, vrt, div);

		{
			// complex
			o_phi = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(phi);
			o_vrt = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(vrt);
			o_div = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(div);
		}

	}


	static
	void benchmark_setup_pvd(
			SphereOperators &opsReal,
			SphereData_Spectral &o_phi,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div
	)
	{
		BenchmarksSphereSWE benchmarks;
		benchmarks.setup(simVars, opsReal, "gaussian_bumps_pvd");

		benchmarks.master->get_initial_state(o_phi, o_vrt, o_div);
	}


	static
	void benchmark_setup_pvd(
			SphereOperators &opsReal,
			sweet::SphereData_SpectralComplex &o_phi,
			sweet::SphereData_SpectralComplex &o_vrt,
			sweet::SphereData_SpectralComplex &o_div
	)
	{

		SphereData_Spectral phi(opsReal.sphereDataConfig);
		SphereData_Spectral vrt(opsReal.sphereDataConfig);
		SphereData_Spectral div(opsReal.sphereDataConfig);

		benchmark_setup_pvd(opsReal, phi, vrt, div);

		{
			// complex
			o_phi = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(phi);
			o_vrt = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(vrt);
			o_div = Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(div);
		}
	}



public:
	void run_tests(
			sweet::SphereDataConfig *sphereDataConfig,
			phys_value_type &alpha
	)
	{
		simVars.outputConfig();

		if (sizeof(phys_value_type) == sizeof(double))
		{
			double_precision_digits = 1e-10*2.0*sphereDataConfig->spectral_modes_m_max;
		}
		else
		{
			double_precision_digits = 1e-7*2.0*sphereDataConfig->spectral_modes_m_max;
		}

		/*
		 * Operators
		 */
		SphereOperators opsReal(sphereDataConfig, &(simVars.sim));
		sphere_operators_type ops(sphereDataConfig, &(simVars.sim));

		/*
		 * Setup local simulation variables
		 */
		double gh0 = simVars.sim.gravitation*simVars.sim.h0;
		double sphere_radius = simVars.sim.sphere_radius;

		double dt_implicit;

		if (gh0 == 0)
			dt_implicit = 1.0/sphere_radius;
		else
			dt_implicit = 1.0/std::sqrt(gh0)*sphere_radius;

		dt_implicit /= 4.0*2.0*sphereDataConfig->spectral_modes_n_max;

		double dt_two_omega = dt_implicit*2.0*simVars.sim.sphere_rotating_coriolis_omega;

		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << "dt_implicit: " << dt_implicit << std::endl;
		std::cout << "h0: " << simVars.sim.h0 << std::endl;
		std::cout << "gh0: " << gh0 << std::endl;
		std::cout << "gravitation: " << simVars.sim.gravitation << std::endl;
		std::cout << "sphere_rotating_coriolis_omega: " << simVars.sim.sphere_rotating_coriolis_omega << std::endl;
		std::cout << "sphere_radius: " << simVars.sim.sphere_radius << std::endl;


		/*
		 * Tests for stationary solution
		 */
		if (simVars.sim.sphere_rotating_coriolis_omega == 0.0)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* SKIPPING 'STATIONARY SOLUTION TESTS' (coriolis = 0 => only trivial stationary case exists)" << std::endl;
			std::cout << "*********************************************************" << std::endl;
		}
		else
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* STATIONARY SOLUTION TESTS" << std::endl;
			std::cout << "*********************************************************" << std::endl;


			sphere_data_spec_type phi(sphereDataConfig);
			sphere_data_spec_type vrt(sphereDataConfig);
			sphere_data_spec_type div(sphereDataConfig);

			benchmark_setup_geostrophic_balance(opsReal, phi, vrt, div);

			sphere_data_spec_type foo, rhs;
			sph_banded_solver_type sphSolverTest;
			double dt = dt_implicit;
			double scalar = 1.234;
			double err = -1;

			double phi_max = phi.toPhys().physical_reduce_max_abs();
			double vrt_max = vrt.toPhys().physical_reduce_max_abs();
			double div_max = div.toPhys().physical_reduce_max_abs();

			if (phi_max == 0)
				phi_max = 1;

			if (vrt_max == 0)
				vrt_max = phi_max/(simVars.sim.sphere_radius*simVars.sim.sphere_radius);

			if (div_max == 0)
				div_max = phi_max/(simVars.sim.sphere_radius*simVars.sim.sphere_radius);

			std::cout << std::endl;
			std::cout << "phi_max: " << phi_max << std::endl;
			std::cout << "vrt_max: " << vrt_max << std::endl;
			std::cout << "div_max: " << div_max << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " J(x)-x" << std::endl;
			{
				foo = ops.implicit_J(vrt, dt_two_omega) - vrt;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				foo = ops.implicit_J(phi, dt_two_omega) -phi;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " x-Jinv(x)" << std::endl;
			{
				foo =vrt - ops.implicit_Jinv(vrt, dt_two_omega);
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << "vrt - Jinv(vrt) = 0 = " << err << std::endl;
				check_error(err, vrt_max);

				foo =phi - ops.implicit_Jinv(phi, dt_two_omega);
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << "phi - Jinv(phi) = 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " FJinv(x)-F(Jinv(x))" << std::endl;
			{
				foo = ops.implicit_FJinv(vrt, dt_two_omega) - ops.implicit_F(ops.implicit_Jinv(vrt, dt_two_omega), dt_two_omega);
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << "FJinv(vrt) - F(Jinv(vrt)) = 0 = " << err << std::endl;
				check_error(err, vrt_max);

				foo = ops.implicit_FJinv(phi, dt_two_omega) - ops.implicit_F(ops.implicit_Jinv(phi, dt_two_omega), dt_two_omega);
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << "FJinv(phi) - F(Jinv(phi)) = 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " -Linv(F(vrt)) - phi = 0" << std::endl;
			{
				foo = -ops.implicit_Linv(ops.implicit_F(vrt, dt_two_omega), dt) - phi;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << "-Linv(F*vrt) - phi = 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " -F(vrt) - L(phi) = 0" << std::endl;
			{
				foo = -ops.implicit_F(vrt, dt_two_omega) - ops.implicit_L(phi, dt);
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << "-F*vrt - L*phi = 0 = " << err << std::endl;
				check_error(err, div_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " invert(Jinv)" << std::endl;
			{
				sphSolverTest.setup(sphereDataConfig, 1);
				sphSolverTest.solver_component_implicit_J(dt_two_omega);

				rhs = ops.implicit_J(vrt, dt_two_omega);
				foo = sphSolverTest.solve(rhs) -vrt;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				rhs = ops.implicit_J(phi, dt_two_omega);
				foo = sphSolverTest.solve(rhs) -phi;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " invert(F+I)" << std::endl;
			{
				sphSolverTest.setup(sphereDataConfig, 2);
				sphSolverTest.solver_component_implicit_F(dt_two_omega);
				sphSolverTest.solver_component_implicit_I(scalar);

				rhs = ops.implicit_F(vrt, dt_two_omega) + scalar*vrt;
				foo = sphSolverTest.solve(rhs) -vrt;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				rhs = ops.implicit_F(phi, dt_two_omega) + scalar*phi;
				foo = sphSolverTest.solve(rhs) -phi;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " b = F(a)+J(a)" << std::endl;
			std::cout << " c = FJinv(b)+J(b)" << std::endl;
			std::cout << " b' = invert(FJinv+J)*c" << std::endl;
			std::cout << " a' = invert(F+J)*b'" << std::endl;
			{
				sph_banded_solver_type sphSolverTestF_J;
				sphSolverTestF_J.setup(sphereDataConfig, 4);
				sphSolverTestF_J.solver_component_implicit_F(dt_two_omega);
				sphSolverTestF_J.solver_component_implicit_J(dt_two_omega);

				sph_banded_solver_type sphSolverTestFJinv_J;
				sphSolverTestFJinv_J.setup(sphereDataConfig, 4);
				sphSolverTestFJinv_J.solver_component_implicit_FJinv(dt_two_omega);
				sphSolverTestFJinv_J.solver_component_implicit_J(dt_two_omega);

				sphere_data_spec_type a = vrt;
				sphere_data_spec_type b = ops.implicit_F(a, dt_two_omega) + ops.implicit_J(a, dt_two_omega);
				sphere_data_spec_type c = ops.implicit_FJinv(b, dt_two_omega) + ops.implicit_J(b, dt_two_omega);
				sphere_data_spec_type b_ = sphSolverTestFJinv_J.solve(c);
				sphere_data_spec_type a_ = sphSolverTestF_J.solve(b_);

				err = (b-b_).toPhys().physical_reduce_max_abs();
				std::cout << " + b_: 0 = " << err << std::endl;

				err = (b-b_).toPhys().physical_reduce_max_abs();
				std::cout << " + a_: 0 = " << err << std::endl;
				check_error(err, a.toPhys().physical_reduce_max_abs());
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " invert(FJinvF+J)" << std::endl;
			{
				sphSolverTest.setup(sphereDataConfig, 4);
				sphSolverTest.solver_component_implicit_FJinvF(dt_two_omega);
				sphSolverTest.solver_component_implicit_J(dt_two_omega);

				rhs = ops.implicit_FJinv(ops.implicit_F(vrt, dt_two_omega), dt_two_omega) + ops.implicit_J(vrt, dt_two_omega);
				foo = sphSolverTest.solve(rhs) - vrt;

				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				rhs = ops.implicit_FJinv(ops.implicit_F(phi, dt_two_omega), dt_two_omega) + ops.implicit_J(phi, dt_two_omega);
				foo = sphSolverTest.solve(rhs) - phi;
				err = foo.toPhys().physical_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " FIN" << std::endl;
			std::cout << "============================================" << std::endl;
		}



		/*
		 * Tests with time-varying solution
		 */
		for (int i = 0; i < 8; i++)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* TIME VARYING SOLUTION TESTS " << i << std::endl;
			std::cout << "*********************************************************" << std::endl;


			sphere_data_spec_type phi(sphereDataConfig);
			sphere_data_spec_type vrt(sphereDataConfig);
			sphere_data_spec_type div(sphereDataConfig);

			benchmark_setup_pvd(opsReal, phi, vrt, div);

			double phi_order = phi.toPhys().physical_reduce_max_abs();
			double vrt_order = vrt.toPhys().physical_reduce_max_abs();
			double div_order = div.toPhys().physical_reduce_max_abs();

			if ((i & 1) == 0)
			{
				phi.spectral_set_zero();
				phi_order = 1.0;
				std::cout << " + phi: set to zero" << std::endl;
			}
			else
			{
				std::cout << " + phi: " << phi_order << std::endl;
			}

			if ((i & 2) == 0)
			{
				vrt.spectral_set_zero();
				vrt_order = phi_order/(simVars.sim.sphere_radius*simVars.sim.sphere_radius);
				std::cout << " + vrt: set to zero" << std::endl;
			}
			else
			{
				std::cout << " + vrt: " << vrt_order << std::endl;
			}

			if ((i & 4) == 0)
			{
				div.spectral_set_zero();
				div_order = phi_order/(simVars.sim.sphere_radius*simVars.sim.sphere_radius);
				std::cout << " + div: set to zero" << std::endl;
			}
			else
			{
				std::cout << " + div: " << div_order << std::endl;
			}


			std::cout << std::endl;
			std::cout << "phi_order: " << phi_order << std::endl;
			std::cout << "vrt_order: " << vrt_order << std::endl;
			std::cout << "div_order: " << div_order << std::endl;


#if 0
			{
				// Run one forward Euler timestep just to get some valid div field
				SWE_Sphere_TS_l_erk l_erk(simVars, ops);
				l_erk.setup(1);
				l_erk.run_timestep(phi, vrt, div, dt_implicit);
			}
#endif


			double gh0 = simVars.sim.gravitation*simVars.sim.h0;
			double sphere_radius = simVars.sim.sphere_radius;

			double dt_two_omega = dt_implicit*2.0*simVars.sim.sphere_rotating_coriolis_omega;

			sph_banded_solver_type sphSolverTest;

			sphere_data_spec_type& S1 = vrt;
			sphere_data_spec_type& S2 = div;
			sphere_data_spec_type& S3 = phi;

			sph_banded_solver_type sphSolverDiv;
			sphSolverDiv.setup(sphereDataConfig, 4);
			sphSolverDiv.solver_component_implicit_J(dt_two_omega);
			sphSolverDiv.solver_component_implicit_FJinvF(dt_two_omega);
			sphSolverDiv.solver_component_implicit_L(gh0*dt_implicit, dt_implicit, sphere_radius);

			sphere_data_spec_type rhs = S2 + ops.implicit_FJinv(S1, dt_two_omega) + ops.implicit_L(S3, dt_implicit);
			sphere_data_spec_type div1 = sphSolverDiv.solve(rhs);

			sphere_data_spec_type phi1 = S3 - dt_implicit*gh0*div1;
			sphere_data_spec_type vrt1 = ops.implicit_Jinv(S1 - ops.implicit_F(div1, dt_two_omega), dt_two_omega);

			/*
			 * Testcase for correct solver
			 */
			{
				double gh0 = simVars.sim.gravitation*simVars.sim.h0;
				double dt_two_omega = dt_implicit*2.0*simVars.sim.sphere_rotating_coriolis_omega;

				sphere_data_spec_type foo, rhs;
				sph_banded_solver_type sphSolverTest;
				double err = -1;

				std::cout << "============================================" << std::endl;
				std::cout << " CURL equation" << std::endl;
				{
					foo = ops.implicit_J(vrt1, dt_two_omega) + ops.implicit_F(div1, dt_two_omega) - S1;
					err = foo.toPhys().physical_reduce_max_abs();
					std::cout << " + err: 0 = " << err << std::endl;
					check_error(err, vrt_order);
				}

				std::cout << "============================================" << std::endl;
				std::cout << " DIV equation" << std::endl;
				{
					foo = ops.implicit_J(div1, dt_two_omega) - ops.implicit_F(vrt1, dt_two_omega) - ops.implicit_L(phi1, dt_implicit) - S2;
					err = foo.toPhys().physical_reduce_max_abs();
					std::cout << " + err: 0 = " << err << std::endl;
					check_error(err, div_order);
				}

				std::cout << "============================================" << std::endl;
				std::cout << " GEOPOT equation" << std::endl;
				{
					foo = phi1 + dt_implicit*gh0*div1 - S3;
					err = foo.toPhys().physical_reduce_max_abs();
					std::cout << " + err: 0 = " << err << std::endl;
					check_error(err, phi_order);
				}
				std::cout << std::endl;
			}
			std::cout << "Tests successful" << std::endl;

		}



		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << " + dt_implicit: " << dt_implicit << std::endl;
		std::cout << " + h0: " << simVars.sim.h0 << std::endl;
		std::cout << " + gh0: " << gh0 << std::endl;
		std::cout << " + gravitation: " << simVars.sim.gravitation << std::endl;
		std::cout << " + sphere_rotating_coriolis_omega: " << simVars.sim.sphere_rotating_coriolis_omega << std::endl;
		std::cout << " + sphere_radius: " << simVars.sim.sphere_radius << std::endl;
		std::cout << std::endl;


		std::cout << "*********************************************************" << std::endl;
		std::cout << "All tests successful" << std::endl;
		std::cout << "*********************************************************" << std::endl;
	}
};




template<
	typename phys_value_type = double,
	typename sphere_data_spec_type = SphereData_Spectral,
	typename sphere_data_phys_type = SphereData_Physical,
	typename sphere_operators_type = SphereOperators,
	typename sph_banded_solver_type = double
>
class TestFB
{
	double double_precision_digits = -1.0;

	void check_error(double err, double rel_value)
	{
		if (err > double_precision_digits*rel_value || std::isnan(err))
		{
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
			std::cout << "+ err: " << err << std::endl;
			std::cout << "+ max_err_threshold: " << double_precision_digits*rel_value << std::endl;
			std::cout << "+ rel_value: " << rel_value << std::endl;
			SWEETError("Error too high");
		}
	};

public:
	void run_tests(
			sweet::SphereDataConfig *sphereDataConfig,
			phys_value_type &alpha
	)
	{
		simVars.outputConfig();

		double_precision_digits = 1e-10*2.0*sphereDataConfig->spectral_modes_m_max;

		/*
		 * Operators
		 */
		SphereOperators ops(sphereDataConfig, &(simVars.sim));

		/*
		 * Setup local simulation variables
		 */
		double gh0 = simVars.sim.gravitation*simVars.sim.h0;
		double sphere_radius = simVars.sim.sphere_radius;

		double dt_implicit;

		if (gh0 == 0)
			dt_implicit = 1.0/sphere_radius;
		else
			dt_implicit = 1.0/std::sqrt(gh0)*sphere_radius;

		dt_implicit /= 4.0*2.0*sphereDataConfig->spectral_modes_n_max;

		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << "dt_implicit: " << dt_implicit << std::endl;
		std::cout << "h0: " << simVars.sim.h0 << std::endl;
		std::cout << "gh0: " << gh0 << std::endl;
		std::cout << "gravitation: " << simVars.sim.gravitation << std::endl;
		std::cout << "sphere_rotating_coriolis_omega: " << simVars.sim.sphere_rotating_coriolis_omega << std::endl;
		std::cout << "sphere_radius: " << simVars.sim.sphere_radius << std::endl;

		/*
		 * Test forward/backward time stepping
		 *
		 * (I - dtexpl*L)^-1 (I + dtimpl*L) = I
		 */
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* FORWARD / BACKWARD TIME INTEGRATION" << std::endl;
			std::cout << "*********************************************************" << std::endl;

			BenchmarksSphereSWE benchmarks;
			benchmarks.setup(simVars, ops, "gaussian_bumps_pvd");

			SphereData_Spectral phi(sphereDataConfig);
			SphereData_Spectral vrt(sphereDataConfig);
			SphereData_Spectral div(sphereDataConfig);

			benchmarks.master->get_initial_state(phi, vrt, div);

			//phi.spectral_set_zero();
			vrt.spectral_set_zero();
			div.spectral_set_zero();

			double dt_implicit = 1.0/std::sqrt(gh0)*sphere_radius;
			dt_implicit = 0.5*sphereDataConfig->spectral_modes_n_max;
			std::cout << "dt_implicit: " << dt_implicit << std::endl;

			if (0)
			{
				// Run one forward Euler timestep just to get some valid div field
				SWE_Sphere_TS_l_erk l_erk(simVars, ops);
				l_erk.setup(1);
				l_erk.run_timestep(phi, vrt, div, dt_implicit);
			}

			double phi_max = phi.toPhys().physical_reduce_max_abs();
			double vrt_max = vrt.toPhys().physical_reduce_max_abs();
			double div_max = div.toPhys().physical_reduce_max_abs();

			if (phi_max == 0)
				phi_max = 1;

			if (vrt_max == 0)
				vrt_max = phi_max/(simVars.sim.sphere_radius*simVars.sim.sphere_radius);

			if (div_max == 0)
				div_max = phi_max/(simVars.sim.sphere_radius*simVars.sim.sphere_radius);

			std::cout << std::endl;
			std::cout << "phi_max: " << phi_max << std::endl;
			std::cout << "vrt_max: " << vrt_max << std::endl;
			std::cout << "div_max: " << div_max << std::endl;
			std::cout << std::endl;

			{
				SphereData_Spectral test_phi = phi;
				SphereData_Spectral test_vrt = vrt;
				SphereData_Spectral test_div = div;

				double dt_explicit = dt_implicit;
				double dt_implicit = -dt_explicit;		// Backward with flipped minus sigh

				SWE_Sphere_TS_l_erk l_erk(simVars, ops);
				l_erk.setup(1);
				l_erk.run_timestep(test_phi, test_vrt, test_div, dt_explicit);

				SWE_Sphere_TS_l_irk l_irk(simVars, ops);
				l_irk.setup(1, dt_implicit);
				l_irk.run_timestep(test_phi, test_vrt, test_div, dt_implicit);

				double err_phi = (test_phi-phi).toPhys().physical_reduce_max_abs();
				double err_vrt = (test_vrt-vrt).toPhys().physical_reduce_max_abs();
				double err_div = (test_div-div).toPhys().physical_reduce_max_abs();

				std::cout << "err_phi: " << err_phi << std::endl;
				check_error(err_phi, phi_max);

				std::cout << "err_vrt: " << err_vrt << std::endl;
				check_error(err_vrt, vrt_max);

				std::cout << "err_div: " << err_div << std::endl;
				check_error(err_div, div_max);

			}
		}

		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << " + dt_implicit: " << dt_implicit << std::endl;
		std::cout << " + h0: " << simVars.sim.h0 << std::endl;
		std::cout << " + gh0: " << gh0 << std::endl;
		std::cout << " + gravitation: " << simVars.sim.gravitation << std::endl;
		std::cout << " + sphere_rotating_coriolis_omega: " << simVars.sim.sphere_rotating_coriolis_omega << std::endl;
		std::cout << " + sphere_radius: " << simVars.sim.sphere_radius << std::endl;
		std::cout << std::endl;


		std::cout << "*********************************************************" << std::endl;
		std::cout << "All tests successful" << std::endl;
		std::cout << "*********************************************************" << std::endl;
	}
};




int main(
		int i_argc,
		char *const i_argv[]
)
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	// Stop overriding simulation variables for this test case
	simVars.benchmark.benchmark_override_simvars = false;

	if (simVars.disc.space_res_spectral[0] == 0)
		SWEETError("Set number of spectral modes to use SPH!");

	sweet::SphereDataConfig sphereDataConfig;
	sphereDataConfig.setupAutoPhysicalSpace(
					simVars.disc.space_res_spectral[0],
					simVars.disc.space_res_spectral[1],
					&simVars.disc.space_res_physical[0],
					&simVars.disc.space_res_physical[1],
					simVars.misc.reuse_spectral_transformation_plans,
					simVars.misc.verbosity,
					simVars.parallelization.num_threads_space

			);

#if 1
	{
		/*
		 * Real-valued solver
		 *
		 * E.g. used for backward Euler
		 */
		double alpha_real = 3.0;

		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "* REAL" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		Test<
			double,
			SphereData_Spectral,
			SphereData_Physical,
			SphereOperators,
			SphBandedMatrixPhysicalReal<std::complex<double>>
		>t_real;

		t_real.run_tests(&sphereDataConfig, alpha_real);

		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "* REAL FB" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		TestFB<
			double,
			SphereData_Spectral,
			SphereData_Physical,
			SphereOperators,
			SphBandedMatrixPhysicalReal<std::complex<double>>
		>t_real_fb;

		t_real_fb.run_tests(&sphereDataConfig, alpha_real);
	}
#endif

#if 1
	{
		/*
		 * Complex-valued solver
		 *
		 * E.g. used for REXI
		 */
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "* COMPLEX" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;

		std::complex<double> alpha_complex(1.0, 3.0);

		Test<
			std::complex<double>,
			SphereData_SpectralComplex,
			SphereData_PhysicalComplex,
			SphereOperatorsComplex,
			SphBandedMatrixPhysicalComplex<std::complex<double>>
		>t_complex;

		t_complex.run_tests(&sphereDataConfig, alpha_complex);
	}
#endif

}

