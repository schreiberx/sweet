
//#if !SWEET_USE_SPECTRAL_SPACE
//	#error "Spectral space not activated"
//#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#if SWEET_USE_SPECTRAL_DEALIASING
#	warning	"Aliasing not working"
#endif


#include <sweet/DataArray.hpp>
#include <sweet/SimulationParameters.hpp>
#include <sweet/Operators2D.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

SimulationParameters parameters;

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

	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);

	SimulationParameters parameters;
	parameters.use_spectral_diffs = 1;
	parameters.setup(i_argc, i_argv);

	if (parameters.use_spectral_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;


	double prev_error_diff_x = 0;
	double prev_error_diff_y = 0;
	double prev_error_diff_z = 0;

	double prev_error_diff2_x = 0;
	double prev_error_diff2_y = 0;

	double freq_x = 2.0;
	double freq_y = 4.0;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = parameters.res[0];
	std::size_t res_y = parameters.res[1];

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
		double eps_convergence = 1e-4;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		parameters.res[0] = res[0];
		parameters.res[1] = res[1];
		parameters.reset();

		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		DataArray<2> h(res);


		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "> Resolution (" << res_x << "x" << res_y << ")" << std::endl;
			std::cout << "> Domain size (" << parameters.sim_domain_size[0] << "x" << parameters.sim_domain_size[1] << ")" << std::endl;
			std::cout << "**********************************************" << std::endl;
			std::cout << "error tol = " << eps << std::endl;
			std::cout << "**********************************************" << std::endl;

			DataArray<2> zero(res);
			DataArray<2> two(res);
			DataArray<2> five(res);

			zero.setAll(0);
			two.setAll(2);
			five.setAll(5);
			h.setAll(0);

			double res2 = (double)(res[0]*res[1]);

			double add_test_two = (zero+two).reduce_rms_quad();
			double add_test_seven = (five+two).reduce_rms_quad();
			double add_test_ten = ((five+two)+3.0).reduce_rms_quad();
			double error = 0;

			error = std::abs(add_test_two-2.0);
			std::cout << "Add test two ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			error = std::abs(add_test_seven-7.0);
			std::cout << "Add test seven ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			error = std::abs(add_test_ten-10.0);
			std::cout << "Add test ten ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			// create sinus curve
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i)/(double)res[0];
					double y = ((double)j)/(double)res[1];

					h.set(j, i, sin(2.0*M_PIl*x)*cos(2.0*M_PIl*y));
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
				exit(-1);
			}

			double sin_test_six = (h+6.0).reduce_sum_quad()/res2;
			error = std::abs(sin_test_six-6.0);
			std::cout << "Sin test add six ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED Sin test add six ||_2 with error " << error << std::endl;
				std::cout << "FAILED with error " << sin_test_six << std::endl;
				exit(-1);
			}

			double sin_test_zero_mul = (h*two).reduce_sum_quad()/res2;
			error = sin_test_zero_mul;
			std::cout << "Sin test times 2 ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}
		}

		std::cout << "TEST A: DONE" << std::endl;



		/**
		 * Tests for basic operators which are not amplifying the solution depending on the domain size
		 */
		{
			DataArray<2> u(res);
			DataArray<2> v(res);

			Operators2D op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs);

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

#define FUN_ID	1

	#if FUN_ID==1
					u.set(j, i, sin(freq_x*M_PIl*x));
					v.set(j, i, cos(freq_y*M_PIl*y));
	#elif FUN_ID==2
					u.set(j, i, sin(freq_x*M_PIl*x));
					v.set(j, i, 1.0/(cos(freq_y*M_PIl*y)+2.0));
	#endif

					h.set(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==2
						sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)
	#endif
					);
				}
			}

			// force forward/backward conversion
			u.requestDataInSpectralSpace();
			u.array_data_cartesian_space_valid = false;

			// force forward/backward conversion
			v.requestDataInSpectralSpace();
			v.array_data_cartesian_space_valid = false;

			double err_z = (u*v-h).reduce_rms_quad();

			if (parameters.use_spectral_diffs)
			{
				std::cout << "error (mul*mul-fun) = " << err_z << std::endl;

#if FUN_ID == 1
				if (err_z > eps)
				{
					std::cerr << "SPEC: Error threshold exceeded for err_z!" << std::endl;
					exit(-1);
				}
#endif
				double err3_laplace =
					(
							h-
							(op.diff2_c_x(h)+op.diff2_c_y(h)).
								spec_div_element_wise(op.diff2_c_x+op.diff2_c_y)
					).reduce_rms_quad();

				std::cout << "SPEC: Error threshold for Laplace and its inverse: " << err3_laplace << std::endl;
				if (err3_laplace > eps)
				{
					std::cerr << "SPEC: Error threshold for Laplace too high for spectral differentiation!" << std::endl;
					exit(-1);
				}
			}
			else
			{

#if FUN_ID == 1
				if (err_z > eps)
				{
					std::cerr << "SPEC: Error threshold exceeded for err_z!" << std::endl;
					exit(-1);
				}
#endif

				prev_error_diff_z = err_z;
			}
		}

		std::cout << "TEST B: DONE" << std::endl;


		/**
		 * Tests for 1st order differential operator
		 */
		{
			DataArray<2> u(res);
			DataArray<2> v(res);
			DataArray<2> h_diff_x(res);
			DataArray<2> h_diff_y(res);

			Operators2D op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs);

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

	#if FUN_ID==1
					u.set(j, i, sin(freq_x*M_PIl*x));
					v.set(j, i, cos(freq_y*M_PIl*y));
	#elif FUN_ID==2
					u.set(j, i, sin(freq_x*M_PIl*x));
					v.set(j, i, 1.0/(cos(freq_y*M_PIl*y)+2.0));
	#endif

					h.set(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==2
						sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)
	#endif
					);

					h_diff_x.set(
						j, i,
	#if FUN_ID==1
						freq_x*M_PIl*cos(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(double)parameters.sim_domain_size[0]
	#elif FUN_ID==2
						2.0*sin(freq_x*M_PIl*x)*std::pow(cos(freq_y*M_PIl*y),2.0)*freq_x*M_PIl*cos(freq_x*M_PIl*x)/(double)parameters.sim_domain_size[0]
	#elif FUN_ID==3
						freq_x*M_PIl*cos(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)/(double)parameters.sim_domain_size[0]
	#endif
					);

					h_diff_y.set(
						j, i,
	#if FUN_ID==1
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/(double)parameters.sim_domain_size[1]
	#elif FUN_ID==2
						-2.0*std::pow(std::sin(freq_x*M_PIl*x),2.0)*std::cos(freq_y*M_PIl*y)*freq_y*M_PIl*std::sin(freq_y*M_PIl*y)/(double)parameters.sim_domain_size[1]
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/pow(cos(freq_y*M_PIl*y)+2.0, 2.0)/(double)parameters.sim_domain_size[1]
	#endif
					);
				}
			}

			// force forward/backward conversion
			u.requestDataInSpectralSpace();
			u.array_data_cartesian_space_valid = false;

			// force forward/backward conversion
			v.requestDataInSpectralSpace();
			v.array_data_cartesian_space_valid = false;

			double res_normalization = sqrt(1.0/parameters.res2_dbl);

			// normalization for diff = 2 pi / L
			double err_x = (op.diff_c_x(h)-h_diff_x).reduce_norm2()*res_normalization*parameters.sim_domain_size[0]/(2.0*M_PIl);
			double err_y = (op.diff_c_y(h)-h_diff_y).reduce_norm2()*res_normalization*parameters.sim_domain_size[1]/(2.0*M_PIl);
			double err_z = (u*v-h).reduce_norm2()*res_normalization;

			if (parameters.use_spectral_diffs)
			{
				std::cout << "error diff x = " << err_x << std::endl;
				std::cout << "error diff y = " << err_y << std::endl;
				std::cout << "error (mul*mul-fun) = " << err_z << std::endl;
				std::cout << "err tol = " << eps << std::endl;

				if (err_x > eps)
				{
					std::cerr << "SPEC: Error threshold for diff-X too high for spectral differentiation!" << std::endl;
					exit(-1);
				}

				if (err_y > eps)
				{
					std::cerr << "SPEC: Error threshold for diff-Y too high for spectral differentiation!" << std::endl;
					exit(-1);
				}

#if FUN_ID == 1
				if (err_z > eps)
				{
					std::cerr << "SPEC: Error threshold exceeded for err_z, value = " << err_z << std::endl;
					exit(-1);
				}
#endif

				double err_int_x = (h-h_diff_x.spec_div_element_wise(op.diff_c_x)).reduce_norm2_quad()*res_normalization;
				std::cout << "Testing spectral inverse x " << err_int_x << std::endl;

				if (err_int_x > eps)
				{
					std::cerr << "SPEC: Error threshold for integration in x too high for spectral integration!" << std::endl;
					std::cout << err_int_x << std::endl;
					exit(-1);
				}

				double err_int_y = (h-h_diff_y.spec_div_element_wise(op.diff_c_y)).reduce_norm2_quad()*res_normalization;
				std::cout << "Testing spectral inverse y " << err_int_y << std::endl;

				if (err_int_y > eps)
				{
					std::cout << err_int_y << std::endl;
					std::cerr << "SPEC: Error threshold for integration in y too high for spectral integration!" << std::endl;
					exit(-1);
				}
			}
			else
			{
				double conv_x = prev_error_diff_x/err_x;
				double conv_y = prev_error_diff_y/err_y;
				double conv_z = prev_error_diff_z/err_z;
				std::cout << "error diff x = " << err_x << std::endl;
				std::cout << "error diff y = " << err_y << std::endl;
				std::cout << "error diff z = " << err_z << std::endl;
				std::cout << "conv x = " << conv_x << std::endl;
				std::cout << "conv y = " << conv_y << std::endl;
				std::cout << "conv z = " << conv_z << std::endl;

				if (conv_x != 0)
				if (abs(conv_x-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv_x, no convergence given!" << std::endl;
					exit(-1);
				}

				if (conv_y != 0)
				if (abs(conv_y-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv_y, no convergence given!" << std::endl;
					exit(-1);
				}

				if (abs(err_z) > eps)
				{
					std::cerr << "Cart: Error threshold exceeded for err_z!" << std::endl;
					exit(-1);
				}

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
			DataArray<2> h_diff2_x(res);
			DataArray<2> h_diff2_y(res);

			Operators2D op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs);

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

					h.set(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==2
						sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)
	#elif FUN_ID==3
						sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)
	#endif
					);

					h_diff2_x.set(
						j, i,
	#if FUN_ID==1
						freq_x*freq_x*M_PIl*M_PIl*(-1.0)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(parameters.sim_domain_size[0]*parameters.sim_domain_size[0])
	#elif FUN_ID==2
	//					2.0*sin(freq_x*M_PIl*x)*std::pow(cos(freq_y*M_PIl*y),2.0)*freq_x*M_PIl*cos(freq_x*M_PIl*x)/(double)parameters.sim_domain_size[0]
	#elif FUN_ID==3
	//					freq_x*M_PIl*cos(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)/(double)parameters.sim_domain_size[0]
	#endif
					);

					h_diff2_y.set(
						j, i,
	#if FUN_ID==1
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*freq_y*M_PIl*cos(freq_y*M_PIl*y)/(parameters.sim_domain_size[1]*parameters.sim_domain_size[1])
	#elif FUN_ID==2
	//					-2.0*std::pow(std::sin(freq_x*M_PIl*x),2.0)*std::cos(freq_y*M_PIl*y)*freq_y*M_PIl*std::sin(freq_y*M_PIl*y)/(double)parameters.sim_domain_size[1]
	#elif FUN_ID==3
	//					sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/pow(cos(freq_y*M_PIl*y)+2.0, 2.0)/(double)parameters.sim_domain_size[1]
	#endif
					);
				}
			}

			double normalization = sqrt(1.0/parameters.res2_dbl);

			// diff2 normalization = 4.0 pi^2 / L^2
			double err2_x = (op.diff2_c_x(h)-h_diff2_x).reduce_norm2_quad()*normalization*(parameters.sim_domain_size[0]*parameters.sim_domain_size[0])/(4.0*M_PIl*M_PIl);
			double err2_y = (op.diff2_c_y(h)-h_diff2_y).reduce_norm2_quad()*normalization*(parameters.sim_domain_size[1]*parameters.sim_domain_size[1])/(4.0*M_PIl*M_PIl);

			if (parameters.use_spectral_diffs)
			{
				std::cout << "error diff2 x = " << err2_x << std::endl;
				std::cout << "error diff2 y = " << err2_y << std::endl;

				if (err2_x > eps)
				{
					std::cerr << "SPEC: Error threshold for diff2-X too high for spectral differentiation!" << std::endl;
					exit(-1);
				}

				if (err2_y > eps)
				{
					std::cerr << "SPEC: Error threshold for diff2-Y too high for spectral differentiation!" << std::endl;
					exit(-1);
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
					exit(-1);
				}

				if (conv2_y != 0)
				if (abs(conv2_y-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv2_y, no convergence given!" << std::endl;
					exit(-1);
				}

				prev_error_diff2_x = err2_x;
				prev_error_diff2_y = err2_y;
			}
		}

		std::cout << "TEST D: DONE" << std::endl;
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
