
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#if SWEET_USE_SPECTRAL_DEALIASINGc
#	warning	"Aliasing not working / supported here"
#endif


#include <sweet/DataArray.hpp>
#include <sweet/Complex2DArrayFFT.hpp>
#include <sweet/SimulationParameters.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>
#include <complex>

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

typedef std::complex<double> complex;


int main(int i_argc, char *i_argv[])
{
#if SWEET_DEBUG_MODE
	trapfpe();
#endif

	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);
	std::cout << std::setprecision(6);
	std::cerr << std::setprecision(6);

	SimulationParameters parameters;
	parameters.use_spectral_diffs = 1;
	parameters.setup(i_argc, i_argv);

	if (parameters.use_spectral_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	double freq_x = 2.0;
	double freq_y = 2.0;

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



		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};
//		std::size_t res_aliasing[2] = {res_x*2, res_y*2};


		parameters.res[0] = res[0];
		parameters.res[1] = res[1];
		parameters.reset();

		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		Complex2DArrayFFT h_cart(res);


		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "> Resolution (" << res_x << "x" << res_y << ")" << std::endl;
			std::cout << "> Domain size (" << parameters.sim_domain_size[0] << "x" << parameters.sim_domain_size[1] << ")" << std::endl;
			std::cout << "**********************************************" << std::endl;
			std::cout << "error tol = " << eps << std::endl;
			std::cout << "**********************************************" << std::endl;

			Complex2DArrayFFT zero_cart(res);
			Complex2DArrayFFT two_cart(res);
			Complex2DArrayFFT five_cart(res);

			zero_cart.setAll(0.0, 0.0);
			two_cart.setAll(2.0, 0.0);
			five_cart.setAll(5.0, 0.0);
			h_cart.setAll(0.0, 0.0);

			double res2 = (double)(res[0]*res[1]);

			double add_test_two = (zero_cart+two_cart).reduce_rms_quad();
			double add_test_seven = (five_cart+two_cart).reduce_rms_quad();

			double add_test_four = (five_cart+two_cart).subScalar_Cart(complex(3.0,0.0)).reduce_rms_quad();
			double add_test_four_spec = (five_cart.toSpec()+two_cart.toSpec()).subScalar_Spec(complex(3.0,0.0)).toCart().reduce_rms_quad();

			double add_test_seven_spec = (five_cart.toSpec()+two_cart.toSpec()).toCart().reduce_rms_quad();
			double add_test_ten = ((five_cart+two_cart).addScalar_Cart(complex(3.0, 0.0))).reduce_rms_quad();
			double add_test_ten_spec = ((five_cart.toSpec()+two_cart.toSpec()).addScalar_Spec(complex(3.0, 0.0))).toCart().reduce_rms_quad();
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

			error = std::abs(add_test_seven_spec-7.0);
			std::cout << "Add test seven cart ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			error = std::abs(add_test_four-4.0);
			std::cout << "Add test four ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			error = std::abs(add_test_four_spec-4.0);
			std::cout << "Add test four cart ||_2 = " << error << std::endl;
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

			error = std::abs(add_test_ten_spec-10.0);
			std::cout << "Add test ten spec ||_2 = " << error << std::endl;
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

					h_cart.set(j, i, sin(2.0*M_PIl*x)*cos(2.0*M_PIl*y), 0);
				}
			}

			// TEST summation
			// has to be zero, error threshold unknown
			error = h_cart.reduce_sum_re_quad()/res2;
			std::cout << "Sin test zero ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				std::cout << "ERROR THRESHOLDS ARE UNKNOWN for summation without abs(), may depend on N!!!" << std::endl;
				exit(-1);
			}

			double sin_test_six = (h_cart.addScalar_Cart(complex(6.0,0.0))).reduce_sum_re_quad()/res2;
			error = std::abs(sin_test_six-6.0);
			std::cout << "Sin test add six ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED Sin test add six ||_2 with error " << error << std::endl;
				std::cout << "FAILED with error " << sin_test_six << std::endl;
				exit(-1);
			}

			Complex2DArrayFFT h_spec = h_cart.toSpec();
			Complex2DArrayFFT two_spec = two_cart.toSpec();

			double sin_test_zero_mul = (h_spec*two_spec).reduce_sum_re_quad()/res2;
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
			Complex2DArrayFFT u(res);
			Complex2DArrayFFT v(res);

//			Operators2D op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs);

			Complex2DArrayFFT op_diff2_c_x(res);
			Complex2DArrayFFT op_diff2_c_y(res);
			op_diff2_c_x.op_setup_diff2_x(parameters.sim_domain_size);
			op_diff2_c_y.op_setup_diff2_y(parameters.sim_domain_size);

			Complex2DArrayFFT op_diff_c_x(res);
			op_diff_c_x.op_setup_diff_x(parameters.sim_domain_size);

			Complex2DArrayFFT op_diff_c_y(res);
			op_diff_c_y.op_setup_diff_y(parameters.sim_domain_size);


			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

					u.set(j, i, sin(freq_x*M_PIl*x), 0.0);
					v.set(j, i, cos(freq_y*M_PIl*y), 0.0);

					h_cart.set(
						j, i,
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y),
						0.0
					);
				}
			}

			// force forward/backward conversion
			double err_z = (u*v-h_cart).reduce_rms_quad();

			std::cout << "error (mul*mul-fun) = " << err_z << std::endl;

			if (err_z > eps)
			{
				std::cerr << "SPEC: Error threshold exceeded for err_z!" << std::endl;
				exit(-1);
			}

			double err3_laplace =
				(
						h_cart-
							((op_diff2_c_x+op_diff2_c_y)(h_cart.toSpec())).
							spec_div_element_wise(op_diff2_c_x+op_diff2_c_y).toCart()
				).reduce_rms_quad();

			std::cout << "SPEC: Error threshold for Laplace and its inverse: " << err3_laplace << std::endl;
			if (err3_laplace > eps)
			{
				std::cerr << "SPEC: Error threshold for Laplace too high for spectral differentiation!" << std::endl;
				exit(-1);
			}


			double err3_laplace_check =
				(
						h_cart-
							((op_diff_c_x*op_diff_c_x+op_diff_c_y*op_diff_c_y)(h_cart.toSpec())).
							spec_div_element_wise(op_diff2_c_x+op_diff2_c_y).toCart()
				).reduce_rms_quad();

			std::cout << "SPEC: Error threshold for Laplace and its inverse (check): " << err3_laplace << std::endl;
			if (err3_laplace_check > eps)
			{
				std::cerr << "SPEC: Error threshold for Laplace check too high for spectral differentiation!" << std::endl;
				exit(-1);
			}
		}

		std::cout << "TEST B: DONE" << std::endl;


		/**
		 * Tests for 1st order differential operator
		 */
		{
			Complex2DArrayFFT op_diff_c_x(res);
			op_diff_c_x.op_setup_diff_x(parameters.sim_domain_size);

			Complex2DArrayFFT op_diff_c_y(res);
			op_diff_c_y.op_setup_diff_y(parameters.sim_domain_size);

			Complex2DArrayFFT u(res);
			Complex2DArrayFFT v(res);
			Complex2DArrayFFT h_diff_x(res);
			Complex2DArrayFFT h_diff_y(res);

//			Operators2D op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs);

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];
//					double x = ((double)i)/(double)parameters.res[0];
//					double y = ((double)j)/(double)parameters.res[1];

					u.set(j, i, sin(freq_x*M_PIl*x), 0);
					v.set(j, i, cos(freq_y*M_PIl*y), 0);

					h_cart.set(
						j, i,
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y),
						0
					);

					h_diff_x.set(
						j, i,
						freq_x*M_PIl*cos(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(double)parameters.sim_domain_size[0],
						0
					);

					h_diff_y.set(
						j, i,
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/(double)parameters.sim_domain_size[1],
						0
					);
				}
			}


			double res_normalization = std::sqrt(1.0/parameters.res2_dbl);

			// normalization for diff = 2 pi / L
			double err_x = (op_diff_c_x(h_cart.toSpec()).toCart()-h_diff_x).reduce_norm2_quad()*res_normalization*parameters.sim_domain_size[0]/(2.0*M_PIl);
			double err_y = (op_diff_c_y(h_cart.toSpec()).toCart()-h_diff_y).reduce_norm2_quad()*res_normalization*parameters.sim_domain_size[1]/(2.0*M_PIl);
			double err_z = (u*v-h_cart).reduce_norm2_quad()*res_normalization;

			std::cout << "error diff x = " << err_x << std::endl;
			std::cout << "error diff y = " << err_y << std::endl;

			Complex2DArrayFFT h_spec = h_cart.toSpec();
			Complex2DArrayFFT h_diff_xy_spec = op_diff_c_x(h_spec) + op_diff_c_y(h_spec);
			Complex2DArrayFFT op_diff_c_x_y = op_diff_c_x + op_diff_c_y;
			Complex2DArrayFFT h_diff_xy_spec_split = op_diff_c_x_y(h_spec);

			double err_xy = (
								h_diff_xy_spec.toCart()
								-h_diff_x-h_diff_y
						).reduce_norm2_quad()*res_normalization/(2.0*M_PIl);

			double err_xy_split = (
								h_diff_xy_spec_split.toCart()
								-h_diff_x-h_diff_y
						).reduce_norm2_quad()*res_normalization/(2.0*M_PIl);

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

			if (err_xy > eps)
			{
				std::cerr << "SPEC: Error threshold for diff-X+Y too high for spectral differentiation!" << std::endl;
				exit(-1);
			}

			if (err_xy_split > eps)
			{
				std::cerr << "SPEC: Error threshold for diff-X+Y split too high for spectral differentiation!" << std::endl;
				exit(-1);
			}

			if (err_z > eps)
			{
				std::cerr << "SPEC: Error threshold exceeded for err_z, value = " << err_z << std::endl;
				exit(-1);
			}

			double err_int_x = (h_cart-h_diff_x.toSpec().spec_div_element_wise(op_diff_c_x).toCart()).reduce_norm2_quad()*res_normalization;
			std::cout << "Testing spectral inverse x " << err_int_x << std::endl;

			if (err_int_x > eps)
			{
				std::cerr << "SPEC: Error threshold for integration in x too high for spectral integration!" << std::endl;
				std::cout << err_int_x << std::endl;
				exit(-1);
			}

			double err_int_y = (h_cart-h_diff_y.toSpec().spec_div_element_wise(op_diff_c_y).toCart()).reduce_norm2_quad()*res_normalization;
			std::cout << "Testing spectral inverse y " << err_int_y << std::endl;

			if (err_int_y > eps)
			{
				std::cout << err_int_y << std::endl;
				std::cerr << "SPEC: Error threshold for integration in y too high for spectral integration!" << std::endl;
				exit(-1);
			}
		}

		std::cout << "TEST C: DONE" << std::endl;


#if 1

		/**
		 * 2nd order differential operator
		 *
		 * note, that the function on which the 2nd diff operator is computed on has
		 * to be scaled up be a factor of domain_size^2, since e.g.
		 *
		 * diff(sin(2 pi x / size), x, x) = 4.0 pi^2 sin(2 pi x / size) / size^2
		 */
		{
			Complex2DArrayFFT h_diff2_x(res);
			Complex2DArrayFFT h_diff2_y(res);

//			Operators2D op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs);

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

					h_cart.set(
						j, i,
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y),
						0
					);

					h_diff2_x.set(
						j, i,
						freq_x*freq_x*M_PIl*M_PIl*(-1.0)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(parameters.sim_domain_size[0]*parameters.sim_domain_size[0]),
						0
					);

					h_diff2_y.set(
						j, i,
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*freq_y*M_PIl*cos(freq_y*M_PIl*y)/(parameters.sim_domain_size[1]*parameters.sim_domain_size[1]),
						0
					);
				}
			}

			double normalization = sqrt(1.0/parameters.res2_dbl);

			Complex2DArrayFFT op_diff2_c_x(res);
			Complex2DArrayFFT op_diff2_c_y(res);
			op_diff2_c_x.op_setup_diff2_x(parameters.sim_domain_size);
			op_diff2_c_y.op_setup_diff2_y(parameters.sim_domain_size);

			// diff2 normalization = 4.0 pi^2 / L^2
			double err2_x = (op_diff2_c_x(h_cart.toSpec()).toCart()-h_diff2_x).reduce_norm2_quad()*normalization*(parameters.sim_domain_size[0]*parameters.sim_domain_size[0])/(4.0*M_PIl*M_PIl);
			double err2_y = (op_diff2_c_y(h_cart.toSpec()).toCart()-h_diff2_y).reduce_norm2_quad()*normalization*(parameters.sim_domain_size[1]*parameters.sim_domain_size[1])/(4.0*M_PIl*M_PIl);

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
#endif

		std::cout << "TEST D: DONE" << std::endl;
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
