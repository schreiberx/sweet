
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <sweet/DataArray.hpp>
#include <sweet/Complex2DArrayFFT.hpp>
#include <sweet/SimulationVariables.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>
#include <complex>

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

typedef std::complex<double> complex;


int main(int i_argc, char *i_argv[])
{
#if SWEET_DEBUG_MODE
	trapfpe();
#endif

	SimulationVariables simVars;
	simVars.disc.use_spectral_diffs = 1;

	const char *bogus_var_names[] = {
			"use-fd-for-complex-array",	/// use finite differences for complex array
			nullptr
	};

	simVars.bogus.var[0] = 0;	// don't use FD per default for complex array

	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "		--use-fd-for-complex-array=[0/1]	Use finite-differences for derivatives in spectral space" << std::endl;
		std::cout << std::endl;
		return -1;
	}

	/*
	 * use finite differences for differential operators in complex array
	 */
	bool use_finite_differences_for_complex_array = simVars.bogus.var[0];

	if (use_finite_differences_for_complex_array)
	{
		std::cout << "********************************************************" << std::endl;
		std::cout << "*** Using finite-differences for complex array" << std::endl;
		std::cout << "********************************************************" << std::endl;
	}

	if (simVars.disc.use_spectral_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	double freq_x = 2.0;
	double freq_y = 2.0;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.res[0];
	std::size_t res_y = simVars.disc.res[1];

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
		double eps_conv = 1e-3*tolerance_increase;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};


		simVars.disc.res[0] = res[0];
		simVars.disc.res[1] = res[1];
		simVars.reset();


		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		Complex2DArrayFFT h_cart(res);


		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "> Resolution (" << res_x << "x" << res_y << ")" << std::endl;
			std::cout << "> Domain size (" << simVars.sim.domain_size[0] << "x" << simVars.sim.domain_size[1] << ")" << std::endl;
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

			double add_test_three_imaginary = ((five_cart+two_cart).addScalar_Cart(complex(-7.0, 3.0))).reduce_rms_quad();
			std::cout << "add_test_three_imaginary ||_2 = " << add_test_three_imaginary << std::endl;
			error = std::abs(add_test_three_imaginary-3.0);
			if (error > eps)
			{
				std::cout << "FAILED add_test_three_imaginary with error " << error;
				exit(-1);
			}

			double add_test_two_two_spec = ((five_cart.toSpec()+two_cart.toSpec()).addScalar_Spec(complex(-7.0-3.0, 4.0))).toCart().reduce_rms_quad();
			std::cout << "add_test_two_two_spec ||_2 = " << add_test_two_two_spec << std::endl;
			error = std::abs(add_test_two_two_spec-5.0);
			if (error > eps)
			{
				std::cout << "FAILED add_test_two_two_spec with error " << error;
				exit(-1);
			}

			double mul_test_two_times_two_spec = ((two_cart.toSpec()*complex(2.0, 0.0))).toCart().reduce_rms_quad();
			std::cout << "mul_test_two_times_two_spec ||_2 = " << mul_test_two_times_two_spec << std::endl;
			error = std::abs(mul_test_two_times_two_spec-4.0);
			if (error > eps)
			{
				std::cout << "FAILED mul_test_two_times_two_spec with error " << error;
				exit(-1);
			}




			// create sinus curve
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
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

//			Operators2D op(parameters.discretization.res, parameters.sim.domain_size, parameters.disc.use_spectral_diffs);

			Complex2DArrayFFT op_diff2_c_x(res);
			Complex2DArrayFFT op_diff2_c_y(res);
			op_diff2_c_x.op_setup_diff2_x(simVars.sim.domain_size, use_finite_differences_for_complex_array);
			op_diff2_c_y.op_setup_diff2_y(simVars.sim.domain_size, use_finite_differences_for_complex_array);

			Complex2DArrayFFT op_diff_c_x(res);
			op_diff_c_x.op_setup_diff_x(simVars.sim.domain_size, use_finite_differences_for_complex_array);

			Complex2DArrayFFT op_diff_c_y(res);
			op_diff_c_y.op_setup_diff_y(simVars.sim.domain_size, use_finite_differences_for_complex_array);


			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.res[0];
					double y = ((double)j+0.5)/(double)simVars.disc.res[1];

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
							(
								((op_diff2_c_x+op_diff2_c_y)(h_cart.toSpec())).
								spec_div_element_wise(op_diff2_c_x+op_diff2_c_y)
							).toCart()
				).reduce_rms_quad();

			std::cout << "SPEC: Error threshold for Laplace and its inverse: " << err3_laplace << std::endl;
			if (err3_laplace > eps)
			{
				std::cerr << "SPEC: Error threshold for Laplace too high for spectral differentiation!" << std::endl;
				exit(-1);
			}

			static double err3_laplace_check_prev = -1;
			double err3_laplace_check =
				(
						h_cart-
							((op_diff_c_x*op_diff_c_x+op_diff_c_y*op_diff_c_y)(h_cart.toSpec())).
							spec_div_element_wise(op_diff2_c_x+op_diff2_c_y).toCart()
				).reduce_rms_quad();

			std::cout << "Error for Laplace (diff*diff()) and its inverse (check): " << err3_laplace_check << std::endl;

			if (use_finite_differences_for_complex_array)
			{
				if (err3_laplace_check_prev != -1)
				{
					double conv = err3_laplace_check_prev/err3_laplace_check;
					std::cout << " + Convergence: " << conv << std::endl;
					if (std::abs(conv-4.0) > eps_conv)
					{
						std::cerr << "Convergence of laplace operator expected to be 4 ... aborting" << std::endl;
						exit(-1);
					}
				}
				err3_laplace_check_prev = err3_laplace_check;
			}
			else
			{
				if (err3_laplace_check > eps)
				{
					std::cerr << "SPEC: Error threshold for Laplace check too high for spectral differentiation!" << std::endl;
					exit(-1);
				}
			}
		}

		std::cout << "TEST B: DONE" << std::endl;


		/**
		 * Tests for 1st order differential operator
		 */
		{
			Complex2DArrayFFT op_diff_c_x(res);
			op_diff_c_x.op_setup_diff_x(simVars.sim.domain_size, use_finite_differences_for_complex_array);

			Complex2DArrayFFT op_diff_c_y(res);
			op_diff_c_y.op_setup_diff_y(simVars.sim.domain_size, use_finite_differences_for_complex_array);

			Complex2DArrayFFT u(res);
			Complex2DArrayFFT v(res);
			Complex2DArrayFFT h_diff_x(res);
			Complex2DArrayFFT h_diff_y(res);

//			Operators2D op(parameters.discretization.res, parameters.sim.domain_size, parameters.disc.use_spectral_diffs);

			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.res[0];
					double y = ((double)j+0.5)/(double)simVars.disc.res[1];
//					double x = ((double)i)/(double)parameters.discretization.res[0];
//					double y = ((double)j)/(double)parameters.discretization.res[1];

					u.set(j, i, sin(freq_x*M_PIl*x), 0);
					v.set(j, i, cos(freq_y*M_PIl*y), 0);

					h_cart.set(
						j, i,
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y),
						0
					);

					h_diff_x.set(
						j, i,
						freq_x*M_PIl*cos(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(double)simVars.sim.domain_size[0],
						0
					);

					h_diff_y.set(
						j, i,
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/(double)simVars.sim.domain_size[1],
						0
					);
				}
			}


			double res_normalization = std::sqrt(1.0/(simVars.disc.res[0]*simVars.disc.res[1]));

			// normalization for diff = 2 pi / L
			double err_x = (op_diff_c_x(h_cart.toSpec()).toCart()-h_diff_x).reduce_norm2_quad()*res_normalization*simVars.sim.domain_size[0]/(2.0*M_PIl);
			double err_y = (op_diff_c_y(h_cart.toSpec()).toCart()-h_diff_y).reduce_norm2_quad()*res_normalization*simVars.sim.domain_size[1]/(2.0*M_PIl);
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

			if (use_finite_differences_for_complex_array)
			{
				static double err_x_prev = -1;
				if (err_x_prev != -1)
				{
					double conv = err_x_prev/err_x;
					std::cout << " + Convergence: " << conv << std::endl;
					if (std::abs(conv-4.0) > eps_conv)
					{
						std::cerr << "Convergence of diff-x operator expected to be 4 ... aborting" << std::endl;
						exit(-1);
					}
				}
				err_x_prev = err_x;

				static double err_y_prev = -1;
				if (err_y_prev != -1)
				{
					double conv = err_y_prev/err_y;
					std::cout << " + Convergence: " << conv << std::endl;
					if (std::abs(conv-4.0) > eps_conv)
					{
						std::cerr << "Convergence of diff-y operator expected to be 4 ... aborting" << std::endl;
						exit(-1);
					}
				}
				err_y_prev = err_y;

				static double err_xy_prev = -1;
				if (err_xy_prev != -1)
				{
					double conv = err_xy_prev/err_xy;
					std::cout << " + Convergence: " << conv << std::endl;
					if (std::abs(conv-4.0) > eps_conv)
					{
						std::cerr << "Convergence of diff-xy operator expected to be 4 ... aborting" << std::endl;
						exit(-1);
					}
				}
				err_xy_prev = err_xy;

			}
			else
			{
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
			Complex2DArrayFFT h_diff2_x(res);
			Complex2DArrayFFT h_diff2_y(res);

//			Operators2D op(parameters.discretization.res, parameters.sim.domain_size, parameters.disc.use_spectral_diffs);

			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)simVars.disc.res[0];
					double y = ((double)j+0.5)/(double)simVars.disc.res[1];

					h_cart.set(
						j, i,
						sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y),
						0
					);

					h_diff2_x.set(
						j, i,
						freq_x*freq_x*M_PIl*M_PIl*(-1.0)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(simVars.sim.domain_size[0]*simVars.sim.domain_size[0]),
						0
					);

					h_diff2_y.set(
						j, i,
						-sin(freq_x*M_PIl*x)*freq_y*M_PIl*freq_y*M_PIl*cos(freq_y*M_PIl*y)/(simVars.sim.domain_size[1]*simVars.sim.domain_size[1]),
						0
					);
				}
			}

			double normalization = sqrt(1.0/(simVars.disc.res[0]*simVars.disc.res[1]));

			Complex2DArrayFFT op_diff2_c_x(res);
			Complex2DArrayFFT op_diff2_c_y(res);
			op_diff2_c_x.op_setup_diff2_x(simVars.sim.domain_size, use_finite_differences_for_complex_array);
			op_diff2_c_y.op_setup_diff2_y(simVars.sim.domain_size, use_finite_differences_for_complex_array);

			// diff2 normalization = 4.0 pi^2 / L^2
			double err2_x = (op_diff2_c_x(h_cart.toSpec()).toCart()-h_diff2_x).reduce_norm2_quad()*normalization*(simVars.sim.domain_size[0]*simVars.sim.domain_size[0])/(4.0*M_PIl*M_PIl);
			double err2_y = (op_diff2_c_y(h_cart.toSpec()).toCart()-h_diff2_y).reduce_norm2_quad()*normalization*(simVars.sim.domain_size[1]*simVars.sim.domain_size[1])/(4.0*M_PIl*M_PIl);

			std::cout << "error diff2 x = " << err2_x << std::endl;
			std::cout << "error diff2 y = " << err2_y << std::endl;


			if (use_finite_differences_for_complex_array)
			{
				static double err2_x_prev = -1;
				if (err2_x_prev != -1)
				{
					double conv = err2_x_prev/err2_x;
					std::cout << " + Convergence: " << conv << std::endl;
					if (std::abs(conv-4.0) > eps_conv)
					{
						std::cerr << "Convergence of diff2-x operator expected to be 4 ... aborting" << std::endl;
						exit(-1);
					}
				}
				err2_x_prev = err2_x;


				static double err2_y_prev = -1;
				if (err2_y_prev != -1)
				{
					double conv = err2_y_prev/err2_y;
					std::cout << " + Convergence: " << conv << std::endl;
					if (std::abs(conv-4.0) > eps_conv)
					{
						std::cerr << "Convergence of diff2-y operator expected to be 4 ... aborting" << std::endl;
						exit(-1);
					}
				}
				err2_y_prev = err2_y;

			}
			else
			{
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
		}

		std::cout << "TEST D: DONE" << std::endl;
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
