
//#if !SWEET_USE_PLANE_SPECTRAL_SPACE
//	#error "Spectral space not activated"
//#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif



#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SWEETError.hpp>

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

	// override flag
	SimulationVariables simVars;
	simVars.disc.space_use_spectral_basis_diffs = true;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.disc.space_use_spectral_basis_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;


	double prev_error_diff_x = 0;
	double prev_error_diff_y = 0;
	//double prev_error_diff_z = 0;

	double prev_error_diff2_x = 0;
	double prev_error_diff2_y = 0;

	double freq_x = 2.0;
	double freq_y = 4.0;

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
		 * upper bound for roundoff errors.
		 *
		 * The upper bound of roundoff errors for a 1D FFT is given by sqrt(N) using
		 * an RMS norm, see "Roundoff Error Analysis of the Fast Fourier Transform" by George U. Ramos
		 *
		 * Let
		 * d: truncated roundoff errors
		 * T: matrix for Fourier transformation (not the FT)
		 *
		 * Then, according to the work mentioned above
		 * |T d|_rms <= sqrt(N) |d|rms
		 *
		 * Rearranging to Euclidian norm, this yields
		 * sqrt(1/N) |T d|_2 <= sqrt(N) * sqrt(1/N) |d|_2
		 *
		 * Cancelling the ugly sqrt terms, we get
		 * |T d|_2 <= sqrt(N) |d|_2
		 *
		 * So the max. error is increasing linearly to sqrt(N) with an increasing
		 * problem size N to be transformed via FT transformed.
		 *
		 * For a 2D problem, we assume the error thresholds to be additive correlated
		 * due to the 2D problem expressed by two successively 1D FFTs, hence the
		 * upper error thresholds should be related in an additive way.
		 *
		 * We can assume, that the error for a 2D problem is increasing linearly with sqrt(N):
		 *
		 * err ~ sqrt(Nx) + sqrt(Ny)
		 */
		double tolerance_increase = sqrt(res_x) + sqrt(res_y);


		double max_aspect = simVars.sim.plane_domain_size[0] / simVars.sim.plane_domain_size[1];
		if (max_aspect < 1.0)
			max_aspect = 1.0/max_aspect;

		tolerance_increase *= max_aspect;

		/*
		 * error tolerance for machine accuracy
		 *
		 * We assume 1e-12 for double precision
		 */
		double eps = 1e-9*tolerance_increase;

		/*
		 * error tolerance for convergence
		 *
		 * Here, we are very patronizing due to flickering convergence for coarse solutions which
		 * are not really representable in the Fouerier space where the discretization errors
		 * are dominating.
		 */
		double eps_convergence = 1e-4*tolerance_increase;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		simVars.disc.space_res_physical[0] = res[0];
		simVars.disc.space_res_physical[1] = res[1];
		simVars.reset();

		planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);


		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		PlaneData h(planeDataConfig);


		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "> Resolution (" << res_x << "x" << res_y << ")" << std::endl;
			std::cout << "> Domain size (" << simVars.sim.plane_domain_size[0] << "x" << simVars.sim.plane_domain_size[1] << ")" << std::endl;
			std::cout << "**********************************************" << std::endl;
			std::cout << "error tol = " << eps << std::endl;
			std::cout << "**********************************************" << std::endl;

			PlaneData zero(planeDataConfig);
			PlaneData two(planeDataConfig);
			PlaneData five(planeDataConfig);

			zero.physical_set_all(0);
			two.physical_set_all(2);
			five.physical_set_all(5);
			h.physical_set_all(0);

			double res2 = (double)(res[0]*res[1]);

			double add_test_two = (zero+two).reduce_rms_quad();
			double add_test_ten2 = (3.0+(five+two)).reduce_rms_quad();
			double add_test_ten3 = (17.0-(five+two)).reduce_rms_quad();
			double error = 0;


			error = std::abs(add_test_two-2.0);
			std::cout << "Add test two ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				SWEETError("EXIT");
			}

			double add_test_seven = (five+two).reduce_rms_quad();
			error = std::abs(add_test_seven-7.0);
			std::cout << "Add test seven ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				SWEETError("EXIT");
			}

			double add_test_ten = ((five+two)+3.0).reduce_rms_quad();

			error = std::abs(add_test_ten-10.0);
			std::cout << "Add test ten ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				SWEETError("EXIT");
			}

			error = std::abs(add_test_ten2-10.0);
			std::cout << "Add test ten2 ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				SWEETError("EXIT");
			}

			error = std::abs(add_test_ten3-10.0);
			std::cout << "Add test ten3 ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				SWEETError("EXIT");
			}

			// create sinus curve
			for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
			{
				for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
				{
					double x = ((double)i)/(double)res[0];
					double y = ((double)j)/(double)res[1];

					h.p_physical_set(j, i, sin(2.0*M_PIl*x)*cos(2.0*M_PIl*y));
				}
			}

			// TEST summation
			// has to be zero, error threshold unknown
			error = h.reduce_sum_quad()/res2;
			std::cout << "Sin test zero ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				std::cout << "ERROR THRESHOLDS ARE UNKNOWN for summation without abs(), may depend on N!!!" << std::endl;
				SWEETError("EXIT");
			}

			double sin_test_six = (h+6.0).reduce_sum_quad()/res2;
			error = std::abs(sin_test_six-6.0);
			std::cout << "Sin test add six ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED Sin test add six ||_2 with error " << error << std::endl;
				std::cout << "FAILED with error " << sin_test_six << std::endl;
				SWEETError("EXIT");
			}

			double sin_test_zero_mul = (h*two).reduce_sum_quad()/res2;
			error = sin_test_zero_mul;
			std::cout << "Sin test times 2 ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				SWEETError("EXIT");
			}
		}

		std::cout << "TEST A: DONE" << std::endl;



		/**
		 * Tests for basic operators which are not amplifying the solution depending on the domain size
		 */
		{
			PlaneData u(planeDataConfig);
			PlaneData v(planeDataConfig);

			PlaneOperators op(
					planeDataConfig,
					simVars.sim.plane_domain_size,
					simVars.disc.space_use_spectral_basis_diffs
			);
			for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
			{
				for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.space_res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.space_res_physical[1];

#define FUN_ID	1

	#if FUN_ID==1
					u.p_physical_set(j, i, sin(freq_x*M_PIl*x));
					v.p_physical_set(j, i, cos(freq_y*M_PIl*y));
	#elif FUN_ID==2
					u.p_physical_set(j, i, sin(freq_x*M_PIl*x));
					v.p_physical_set(j, i, 1.0/(cos(freq_y*M_PIl*y)+2.0));
	#endif

					h.p_physical_set(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==2
						sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)
	#endif
					);

#undef FUN_ID
				}
			}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			// force forward/backward conversion
			u.request_data_spectral();
			u.physical_space_data_valid = false;

			// force forward/backward conversion
			v.request_data_spectral();
			v.physical_space_data_valid = false;
#endif

//			double err_z = (u*v-h).reduce_rms_quad();

			if (simVars.disc.space_use_spectral_basis_diffs)
			{
#if 0
//				std::cout << "error (mul*mul-fun) = " << err_z << std::endl;

#if FUN_ID == 1
				if (err_z > eps || std::isnan(err_z))
				{
					SWEETError("SPEC: Error threshold exceeded for err_z!");
					SWEETError("EXIT");
				}
#endif
#endif

#if SWEET_USE_PLANE_SPECTRAL_SPACE
				double err3_laplace =
					(
							h-
							(op.diff2_c_x(h)+op.diff2_c_y(h)).
								spectral_div_element_wise(op.diff2_c_x+op.diff2_c_y)
					).reduce_rms_quad();

				std::cout << "SPEC: Error threshold for Laplace and its inverse: " << err3_laplace << std::endl;
				if (err3_laplace > eps || std::isnan(err3_laplace))
					SWEETError("SPEC: Error threshold for Laplace too high for spectral differentiation!");
#endif
			}
			else
			{

#if 0
#if FUN_ID == 1
				if (err_z > eps || std::isnan(err_z))
					SWEETError("SPEC: Error threshold exceeded for err_z!");
#endif
#endif

				//prev_error_diff_z = err_z;
			}
		}

		std::cout << "TEST B: DONE" << std::endl;


		/**
		 * Tests for 1st order differential operator
		 */
		{
			PlaneData u(planeDataConfig);
			PlaneData v(planeDataConfig);
			PlaneData h_diff_x(planeDataConfig);
			PlaneData h_diff_y(planeDataConfig);

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
			{
				for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
				{
#define FUN_ID	1
					double x = ((double)i+0.5)/(double)simVars.disc.space_res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.space_res_physical[1];

	#if FUN_ID==1
					u.p_physical_set(j, i, sin(freq_x*M_PIl*x));
					v.p_physical_set(j, i, cos(freq_y*M_PIl*y));
	#elif FUN_ID==2
					u.p_physical_set(j, i, sin(freq_x*M_PIl*x));
					v.p_physical_set(j, i, 1.0/(cos(freq_y*M_PIl*y)+2.0));
	#endif

					h.p_physical_set(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==2
						sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)
	#endif
					);

					h_diff_x.p_physical_set(
						j, i,
	#if FUN_ID==1
						freq_x*M_PIl*cos(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(double)simVars.sim.plane_domain_size[0]
	#elif FUN_ID==2
						2.0*sin(freq_x*M_PIl*x)*std::pow(cos(freq_y*M_PIl*y),2.0)*freq_x*M_PIl*cos(freq_x*M_PIl*x)/(double)simVars.sim.plane_domain_size[0]
	#elif FUN_ID==3
						freq_x*M_PIl*cos(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)/(double)simVars.sim.plane_domain_size[0]
	#endif
					);

					h_diff_y.p_physical_set(
						j, i,
	#if FUN_ID==1
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/(double)simVars.sim.plane_domain_size[1]
	#elif FUN_ID==2
						-2.0*std::pow(std::sin(freq_x*M_PIl*x),2.0)*std::cos(freq_y*M_PIl*y)*freq_y*M_PIl*std::sin(freq_y*M_PIl*y)/(double)simVars.sim.plane_domain_size[1]
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/pow(cos(freq_y*M_PIl*y)+2.0, 2.0)/(double)simVars.sim.plane_domain_size[1]
	#endif
					);
#undef FUN_ID
				}
			}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			// force forward/backward conversion
			u.request_data_spectral();
			u.physical_space_data_valid = false;

			// force forward/backward conversion
			v.request_data_spectral();
			v.physical_space_data_valid = false;
#endif

			double res_normalization = sqrt(1.0/(simVars.disc.space_res_physical[0]*simVars.disc.space_res_physical[1]));

			// normalization for diff = 2 pi / L
			double err_x = (op.diff_c_x(h)-h_diff_x).reduce_norm2()*res_normalization*simVars.sim.plane_domain_size[0]/(2.0*M_PIl);
			double err_y = (op.diff_c_y(h)-h_diff_y).reduce_norm2()*res_normalization*simVars.sim.plane_domain_size[1]/(2.0*M_PIl);
			double err_z = (u*v-h).reduce_norm2()*res_normalization;

			if (simVars.disc.space_use_spectral_basis_diffs)
			{
				std::cout << "error diff x = " << err_x << std::endl;
				std::cout << "error diff y = " << err_y << std::endl;
				std::cout << "error (mul*mul-fun) = " << err_z << std::endl;
				std::cout << "err tol = " << eps << std::endl;

				if (err_x > eps || std::isnan(err_x))
					SWEETError("SPEC: Error threshold for diff-X too high for spectral differentiation!");

				if (err_y > eps || std::isnan(err_y))
					SWEETError("SPEC: Error threshold for diff-Y too high for spectral differentiation!");

#if FUN_ID == 1
				if (err_z > eps || std::isnan(err_z))
				{
					std::cerr << "SPEC: Error threshold exceeded for err_z, value = " << err_z << std::endl;
					SWEETError("EXIT");
				}
#endif

#if SWEET_USE_PLANE_SPECTRAL_SPACE
				double err_int_x = (h-h_diff_x.spectral_div_element_wise(op.diff_c_x)).reduce_norm2_quad()*res_normalization;
				std::cout << "Testing spectral inverse x " << err_int_x << std::endl;

				if (err_int_x > eps || std::isnan(err_int_x))
				{
					std::cout << err_int_x << std::endl;
					SWEETError(std::string("SPEC: Error threshold for integration in x too high for spectral integration! "));
				}

				double err_int_y = (h-h_diff_y.spectral_div_element_wise(op.diff_c_y)).reduce_norm2_quad()*res_normalization;
				std::cout << "Testing spectral inverse y " << err_int_y << std::endl;

				if (err_int_y > eps || std::isnan(err_int_y))
				{
					std::cout << err_int_y << std::endl;
					SWEETError(std::string("SPEC: Error threshold for integration in y too high for spectral integration! "));
				}
#endif
			}
			else
			{
				double conv_x = prev_error_diff_x/err_x;
				double conv_y = prev_error_diff_y/err_y;
				//double conv_z = prev_error_diff_z/err_z;
				std::cout << "error diff x = " << err_x << std::endl;
				std::cout << "error diff y = " << err_y << std::endl;
				//std::cout << "error diff z = " << err_z << std::endl;
				std::cout << "conv x = " << conv_x << std::endl;
				std::cout << "conv y = " << conv_y << std::endl;
				//std::cout << "conv z = " << conv_z << std::endl;

				if (conv_x != 0)
				if (abs(conv_x-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv_x, no convergence given!" << std::endl;
					SWEETError("EXIT");
				}

				if (conv_y != 0)
				if (abs(conv_y-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv_y, no convergence given!" << std::endl;
					SWEETError("EXIT");
				}
#if 0
				if (abs(err_z) > eps)
				{
					std::cerr << "Cart: Error threshold exceeded for err_z!" << std::endl;
					SWEETError("EXIT");
				}
#endif

				prev_error_diff_x = err_x;
				prev_error_diff_y = err_y;
			}
		}

		std::cout << "TEST C: DONE" << std::endl;



		/**
		 * 2nd order differential operator
		 *
		 * note, that the function on which the 2nd diff operator is computed on has
		 * to be scaled up be a factor of domain_size^2, since e.g.
		 *
		 * diff(sin(2 pi x / size), x, x) = 4.0 pi^2 sin(2 pi x / size) / size^2
		 */
		{
			PlaneData h_diff2_x(planeDataConfig);
			PlaneData h_diff2_y(planeDataConfig);

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
			{
				for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.space_res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.space_res_physical[1];
#define FUN_ID 1
					h.p_physical_set(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==2
						sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)
	#endif
					);

					h_diff2_x.p_physical_set(
						j, i,
	#if FUN_ID==1
						freq_x*freq_x*M_PIl*M_PIl*(-1.0)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(simVars.sim.plane_domain_size[0]*simVars.sim.plane_domain_size[0])
	#elif FUN_ID==2
	//					2.0*sin(freq_x*M_PIl*x)*std::pow(cos(freq_y*M_PIl*y),2.0)*freq_x*M_PIl*cos(freq_x*M_PIl*x)/(double)parameters.sim.domain_size[0]
	#elif FUN_ID==3
	//					freq_x*M_PIl*cos(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)/(double)parameters.sim.domain_size[0]
	#endif
					);

					h_diff2_y.p_physical_set(
						j, i,
	#if FUN_ID==1
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*freq_y*M_PIl*cos(freq_y*M_PIl*y)/(simVars.sim.plane_domain_size[1]*simVars.sim.plane_domain_size[1])
	#elif FUN_ID==2
	//					-2.0*std::pow(std::sin(freq_x*M_PIl*x),2.0)*std::cos(freq_y*M_PIl*y)*freq_y*M_PIl*std::sin(freq_y*M_PIl*y)/(double)parameters.sim.domain_size[1]
	#elif FUN_ID==3
	//					sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/pow(cos(freq_y*M_PIl*y)+2.0, 2.0)/(double)parameters.sim.domain_size[1]
	#endif
					);
#undef FUN_ID
				}
			}

			double normalization = sqrt(1.0/(simVars.disc.space_res_physical[0]*simVars.disc.space_res_physical[1]));

			// diff2 normalization = 4.0 pi^2 / L^2
			double err2_x = (op.diff2_c_x(h)-h_diff2_x).reduce_norm2_quad()*normalization*(simVars.sim.plane_domain_size[0]*simVars.sim.plane_domain_size[0])/(4.0*M_PIl*M_PIl);
			double err2_y = (op.diff2_c_y(h)-h_diff2_y).reduce_norm2_quad()*normalization*(simVars.sim.plane_domain_size[1]*simVars.sim.plane_domain_size[1])/(4.0*M_PIl*M_PIl);

			if (simVars.disc.space_use_spectral_basis_diffs)
			{
				std::cout << "error diff2 x = " << err2_x << std::endl;
				std::cout << "error diff2 y = " << err2_y << std::endl;

				if (err2_x > eps)
				{
					std::cerr << "SPEC: Error threshold for diff2-X too high for spectral differentiation!" << std::endl;
					SWEETError("EXIT");
				}

				if (err2_y > eps)
				{
					std::cerr << "SPEC: Error threshold for diff2-Y too high for spectral differentiation!" << std::endl;
					SWEETError("EXIT");
				}
			}
			else
			{
				double conv2_x = prev_error_diff2_x/err2_x;
				double conv2_y = prev_error_diff2_y/err2_y;
				std::cout << "error diff2 x = " << err2_x << std::endl;
				std::cout << "error diff2 y = " << err2_y << std::endl;
				std::cout << "conv2 x = " << conv2_x << std::endl;
				std::cout << "conv2 y = " << conv2_y << std::endl;

				if (conv2_x != 0)
				if (abs(conv2_x-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv2_x, no convergence given!" << std::endl;
					SWEETError("EXIT");
				}

				if (conv2_y != 0)
				if (abs(conv2_y-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv2_y, no convergence given!" << std::endl;
					SWEETError("EXIT");
				}

				prev_error_diff2_x = err2_x;
				prev_error_diff2_y = err2_y;
			}
		}

		std::cout << "TEST D: DONE" << std::endl;



		/**
		 * Tests for helmholtz solver
		 */
		{
			PlaneData u(planeDataConfig);
			PlaneData v(planeDataConfig);

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
			{
				for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.space_res_physical[0];
					double y = ((double)j+0.5)/(double)simVars.disc.space_res_physical[1];

#define FUN_ID	3
					h.p_physical_set(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==2
						10+sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)+10.0
	#endif
					);

#undef FUN_ID
				}
			}

			if (simVars.disc.space_use_spectral_basis_diffs)
			{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
				double kappa = 406.666;
				double a = 2.0;
				double b = 5.0;
				/**
				 * Solve
				 *   ( kappa*h - c * (diff2x(h) + diff2y(h))) =
				 *   ( kappa - c*diff2x - c*diff2y) * h = rhs;
				 */
				PlaneData helmholtz_operator = (-a*b*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(kappa);

				PlaneData rhs =  kappa*h - a*b*(op.diff2_c_x(h) + op.diff2_c_y(h));

				double err3_helmholtz =
					(
							h-rhs.spectral_div_element_wise(helmholtz_operator)
					).reduce_rms_quad();

				std::cout << "SPEC: Error threshold for Helmholtz operator with kappa = " << kappa << " and its inverse: " << err3_helmholtz << std::endl;
				if (err3_helmholtz > eps)
				{
					std::cerr << "SPEC: Error threshold for Laplace too high for spectral differentiation!" << std::endl;
					SWEETError("EXIT");
				}
#endif
			}
		}

		std::cout << "TEST E: DONE" << std::endl;

	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
