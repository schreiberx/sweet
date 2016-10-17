
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <sweet/Stopwatch.hpp>
#include "../include/sweet/plane/PlaneData.hpp"
#include <sweet/SimulationVariables.hpp>
#include "../include/sweet/plane/PlaneOperators.hpp"

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>


SimulationVariables simVars;

#if SWEET_DEBUG_MODE
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
#if SWEET_DEBUG_MODE
	trapfpe();
#endif

	SimulationVariables simVars;
	simVars.disc.use_spectral_basis_diffs = 1;

	const char *bogus_var_names[] =
	{
			nullptr
	};



	if (simVars.disc.use_spectral_basis_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.res_physical[0];
	std::size_t res_y = simVars.disc.res_physical[1];

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


		simVars.disc.res_physical[0] = res[0];
		simVars.disc.res_physical[1] = res[1];
		simVars.reset();


		/*
		 * keep u,v in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		PlaneData u_ana(res);
		PlaneData v_ana(res);
		PlaneData f(res);
		PlaneData g(res);

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
			for (std::size_t j = 0; j < simVars.disc.res_physical[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.res_physical[1];

					// u and v to reconstruct
					u_ana.set(
						j, i,
						sin(2*M_PIl*x)*cos(2*M_PIl*y)
					);

					v_ana.set(
						j, i,
						cos(2*M_PIl*x)*sin(2*M_PIl*y)
					);

					// sources for the right hand side to fulfill the equation system for given u and v
					f.set(
						j, i,
						(1+8*M_PIl*M_PIl)*sin(2*M_PIl*x)*cos(2*M_PIl*y)-2*M_PIl*cos(2*M_PIl*x)*cos(2*M_PIl*y)+2*M_PIl*sin(2*M_PIl*x)*sin(2*M_PIl*y)
					);

					g.set(
						j, i,
						(1+8*M_PIl*M_PIl)*cos(2*M_PIl*x)*sin(2*M_PIl*y)-2*M_PIl*cos(2*M_PIl*x)*cos(2*M_PIl*y)+2*M_PIl*sin(2*M_PIl*x)*sin(2*M_PIl*y)
					);
				}
			}


			PlaneData u(res);
			PlaneData v(res);
			u.set_all(0);
			v.set_all(0);

			PlaneOperators op(simVars.disc.res_phys, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

			PlaneData rhs_u = u_ana;
			PlaneData rhs_v = v_ana;
			f.requestDataInSpectralSpace();
			g.requestDataInSpectralSpace();
			rhs_u = op.diff_c_x(rhs_u)+op.diff_c_y(rhs_u);
			rhs_u += f;
			rhs_v = op.diff_c_x(rhs_v)+op.diff_c_y(rhs_v);
			rhs_v += g;
			PlaneData lhs = (-(op.diff2_c_x + op.diff2_c_y)).addScalar_Cart(1.0);

			Stopwatch watch;
			watch.start();

			u = rhs_u.spec_div_element_wise(lhs);
			v = rhs_v.spec_div_element_wise(lhs);

			watch.stop();

			// compare with analytical solution
			double error_analytical_u = (u-u_ana).reduce_rms();
			double error_analytical_v = (v-v_ana).reduce_rms();

			std::cout << "    Computed the solution in " << watch() << " seconds" << std::endl;
			std::cout << "                           RMS error in u (analytical): " << error_analytical_u << std::endl;
			std::cout << "                           RMS error in v (analytical): " << error_analytical_v << std::endl;

		}
		std::cout << "TEST HELMHOLTZ REAL: DONE" << std::endl;
	}


	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
