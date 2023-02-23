
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <sweet/Stopwatch.hpp>
#include "../include/sweet/plane/PlaneData_Spectral.hpp"
#include <sweet/SimulationVariables.hpp>
#include "../include/sweet/plane/PlaneOperators.hpp"

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;


SimulationVariables simVars;

#if SWEET_DEBUG
#include <fenv.h>
static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked.  */

  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif


int main(int i_argc, char *i_argv[])
{
#if SWEET_DEBUG
	trapfpe();
#endif

	SimulationVariables simVars;
	simVars.disc.space_use_spectral_basis_diffs = 1;

	if (simVars.disc.space_use_spectral_basis_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.space_res_physical[0];
	std::size_t res_y = simVars.disc.space_res_physical[1];

	std::size_t max_res = 2048;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	for (; res_x <= max_res && res_y <= max_res; res_x *= 2, res_y *= 2)
	{
		/*
		 * error tolerance for machine accuracy
		 *
		 * We assume 1e-12 for double precision
		 */
//		double eps = 1e-9*tolerance_increase;
//		double eps_conv = 1e-3*tolerance_increase;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};


		simVars.disc.space_res_physical[0] = res[0];
		simVars.disc.space_res_physical[1] = res[1];
		simVars.reset();

		planeDataConfigInstance.setupAutoSpectralSpaceFromPhysical(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);


		/*
		 * keep u,v in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		PlaneData_Spectral u_ana(planeDataConfig);
		PlaneData_Spectral v_ana(planeDataConfig);
		PlaneData_Spectral f(planeDataConfig);
		PlaneData_Spectral g(planeDataConfig);

		PlaneData_Physical u_ana_phys(planeDataConfig);
		PlaneData_Physical v_ana_phys(planeDataConfig);
		PlaneData_Physical f_phys(planeDataConfig);
		PlaneData_Physical g_phys(planeDataConfig);
		/**
		 * Test iterative solver for Helmholtz problem
		 * u - u_xx + u_yy = u_x + u_y + f
		 * v - v_xx + v_yy = v_x + v_y + g
		 * <=> (I-\grad^2) u = \grad u + Q
		 * with known right hand side. Sources f and g are calculated for
		 * u = sin(2*PI*x)*cos(2*PI*y)
		 * v = cos(2*PI*x)*sin(2*PI*y)
		 */
		{
			for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
			{
				for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.space_res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.space_res_physical[1];

					// u and v to reconstruct
					u_ana_phys.physical_set_value(
						j, i,
						sin(2*M_PI*x)*cos(2*M_PI*y)
					);

					v_ana_phys.physical_set_value(
						j, i,
						cos(2*M_PI*x)*sin(2*M_PI*y)
					);

					// sources for the right hand side to fulfill the equation system for given u and v
					f_phys.physical_set_value(
						j, i,
						(1+8*M_PI*M_PI)*sin(2*M_PI*x)*cos(2*M_PI*y)-2*M_PI*cos(2*M_PI*x)*cos(2*M_PI*y)+2*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)
					);

					g_phys.physical_set_value(
						j, i,
						(1+8*M_PI*M_PI)*cos(2*M_PI*x)*sin(2*M_PI*y)-2*M_PI*cos(2*M_PI*x)*cos(2*M_PI*y)+2*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)
					);
				}
			}

			u_ana.loadPlaneDataPhysical(u_ana_phys);
			v_ana.loadPlaneDataPhysical(v_ana_phys);
			f.loadPlaneDataPhysical(f_phys);
			g.loadPlaneDataPhysical(g_phys);

			PlaneData_Spectral u(planeDataConfig);
			PlaneData_Spectral v(planeDataConfig);
			PlaneData_Physical u_phys(planeDataConfig);
			PlaneData_Physical v_phys(planeDataConfig);
			u_phys.physical_set_zero();
			v_phys.physical_set_zero();
			u.loadPlaneDataPhysical(u_phys);
			v.loadPlaneDataPhysical(v_phys);

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			PlaneData_Spectral rhs_u = u_ana;
			PlaneData_Spectral rhs_v = v_ana;
			rhs_u = op.diff_c_x(rhs_u)+op.diff_c_y(rhs_u);
			rhs_u += f;
			rhs_v = op.diff_c_x(rhs_v)+op.diff_c_y(rhs_v);
			rhs_v += g;
			PlaneData_Spectral lhs = (-(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);

			Stopwatch watch;
			watch.start();

			u = rhs_u.spectral_div_element_wise(lhs);
			v = rhs_v.spectral_div_element_wise(lhs);

			watch.stop();

			// compare with analytical solution
			double error_analytical_u = (u-u_ana).spectral_reduce_rms();
			double error_analytical_v = (v-v_ana).spectral_reduce_rms();

			std::cout << "    Computed the solution in " << watch() << " seconds" << std::endl;
			std::cout << "                           RMS error in u (analytical): " << error_analytical_u << std::endl;
			std::cout << "                           RMS error in v (analytical): " << error_analytical_v << std::endl;

		}
		std::cout << "TEST HELMHOLTZ REAL: DONE" << std::endl;
	}


	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
