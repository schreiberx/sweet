/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <sweet/core/defaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/core/plane/PlaneComplex.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/ProgramArguments.hpp>

#include <ostream>
#include <cmath>


typedef std::complex<double> complex;


int main(
		int i_argc,
		char *i_argv[]
)
{
	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackPlaneDataOps *shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.printShackData();

	double freq_x = 2.0;
	double freq_y = 2.0;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = shackPlaneDataOps->space_res_physical[0];
	std::size_t res_y = shackPlaneDataOps->space_res_physical[1];

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
		//double eps_conv = 1e-3*tolerance_increase;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		shackPlaneDataOps->space_res_physical[0] = res[0];
		shackPlaneDataOps->space_res_physical[1] = res[1];
		shackPlaneDataOps->space_res_spectral[0] = 0;
		shackPlaneDataOps->space_res_spectral[1] = 0;

		sweet::PlaneDataConfig planeDataConfig;
		planeDataConfig.setupAuto(shackPlaneDataOps);
		ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(planeDataConfig);

		sweet::PlaneOperatorsComplex ops;
		ops.setup(planeDataConfig, shackPlaneDataOps);
		ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(ops);


		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		sweet::PlaneData_PhysicalComplex h_cart(planeDataConfig);


		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "> Resolution (" << res_x << "x" << res_y << ")" << std::endl;
			std::cout << "> Domain size (" << shackPlaneDataOps->plane_domain_size[0] << "x" << shackPlaneDataOps->plane_domain_size[1] << ")" << std::endl;
			std::cout << "**********************************************" << std::endl;
			std::cout << "error tol = " << eps << std::endl;
			std::cout << "**********************************************" << std::endl;

			sweet::PlaneData_PhysicalComplex zero_cart(planeDataConfig);
			sweet::PlaneData_PhysicalComplex two_cart(planeDataConfig);
			sweet::PlaneData_PhysicalComplex five_cart(planeDataConfig);

			zero_cart.physical_set_all_value(0.0, 0.0);
			two_cart.physical_set_all_value(2.0, 0.0);
			five_cart.physical_set_all_value(5.0, 0.0);
			h_cart.physical_set_all_value(0.0, 0.0);

			double res2 = (double)(res[0]*res[1]);

			double add_test_two = (zero_cart+two_cart).physical_reduce_rms_quad();
			double add_test_seven = (five_cart+two_cart).physical_reduce_rms_quad();

			double add_test_four = ((five_cart+two_cart)-complex(3.0,0.0)).physical_reduce_rms_quad();
			double add_test_four_spec = ((five_cart+two_cart)-complex(3.0,0.0)).physical_reduce_rms_quad();

			double add_test_seven_spec = (five_cart+two_cart).physical_reduce_rms_quad();
			double add_test_ten = ((five_cart+two_cart)+complex(3.0, 0.0)).physical_reduce_rms_quad();
			double add_test_ten_spec = ((five_cart+two_cart)+complex(3.0, 0.0)).physical_reduce_rms_quad();
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
			std::cout << "Add test seven spec ||_2 = " << error << std::endl;
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

			double add_test_three_imaginary = ((five_cart+two_cart)+complex(-7.0, 3.0)).physical_reduce_rms_quad();
			std::cout << "add_test_three_imaginary ||_2 = " << add_test_three_imaginary << std::endl;
			error = std::abs(add_test_three_imaginary-3.0);
			if (error > eps)
			{
				std::cout << "FAILED add_test_three_imaginary with error " << error;
				exit(-1);
			}

			double add_test_two_two_spec = ((five_cart+two_cart)+complex(-7.0-3.0, 4.0)).physical_reduce_rms_quad();
			std::cout << "add_test_two_two_spec ||_2 = " << add_test_two_two_spec << std::endl;
			error = std::abs(add_test_two_two_spec-5.0);
			if (error > eps)
			{
				std::cout << "FAILED add_test_two_two_spec with error " << error;
				exit(-1);
			}

			double mul_test_two_times_two_spec = ((two_cart*complex(2.0, 0.0))).physical_reduce_rms_quad();
			std::cout << "mul_test_two_times_two_spec ||_2 = " << mul_test_two_times_two_spec << std::endl;
			error = std::abs(mul_test_two_times_two_spec-4.0);
			if (error > eps)
			{
				std::cout << "FAILED mul_test_two_times_two_spec with error " << error;
				exit(-1);
			}




			// create sinus curve
			for (int j = 0; j < shackPlaneDataOps->space_res_physical[1]; j++)
			{
				for (int i = 0; i < shackPlaneDataOps->space_res_physical[0]; i++)
				{
					double x = ((double)i)/(double)res[0];
					double y = ((double)j)/(double)res[1];

					h_cart.physical_set_value(j, i, (double)(sin(2.0*M_PI*x)*cos(2.0*M_PI*y)), 0.0);
				}
			}

			// TEST summation
			// has to be zero, error threshold unknown
			error = h_cart.physical_reduce_sum_re_quad()/res2;
			std::cout << "Sin test zero ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				std::cout << "ERROR THRESHOLDS ARE UNKNOWN for summation without abs(), may depend on N!!!" << std::endl;
				exit(-1);
			}

			double sin_test_six = (h_cart+complex(6.0,0.0)).physical_reduce_sum_re_quad()/res2;
			error = std::abs(sin_test_six-6.0);
			std::cout << "Sin test add six ||_2 = " << error << std::endl;
			if (error > eps)
			{
				std::cout << "FAILED Sin test add six ||_2 with error " << error << std::endl;
				std::cout << "FAILED with error " << sin_test_six << std::endl;
				exit(-1);
			}

			sweet::PlaneData_PhysicalComplex h_spec = h_cart;
			sweet::PlaneData_PhysicalComplex two_spec = two_cart;

			double sin_test_zero_mul = (h_spec*two_spec).physical_reduce_sum_re_quad()/res2;
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
			sweet::PlaneData_PhysicalComplex u(planeDataConfig);
			sweet::PlaneData_PhysicalComplex v(planeDataConfig);

			for (int j = 0; j < shackPlaneDataOps->space_res_physical[1]; j++)
			{
				for (int i = 0; i < shackPlaneDataOps->space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)shackPlaneDataOps->space_res_physical[0];
					double y = ((double)j+0.5)/(double)shackPlaneDataOps->space_res_physical[1];

					u.physical_set_value(j, i, sin(freq_x*M_PI*x), 0.0);
					v.physical_set_value(j, i, cos(freq_y*M_PI*y), 0.0);

					h_cart.physical_set_value(
						j, i,
						sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y),
						0.0
					);
				}
			}

			// force forward/backward conversion
			double err_z = (u*v-h_cart).physical_reduce_rms_quad();

			std::cout << "error (mul*mul-fun) = " << err_z << std::endl;

			if (err_z > eps)
			{
				std::cerr << "SPEC: Error threshold exceeded for err_z!" << std::endl;
				exit(-1);
			}
			sweet::PlaneData_SpectralComplex h_spec(h_cart.planeDataConfig);
			h_spec.loadPlaneDataPhysical(h_cart);
			double err3_laplace =
				(
						h_cart-
							(
								//((ops.diff2_c_x+ops.diff2_c_y)(h_spec)).
								//spectral_div_element_wise(ops.diff2_c_x+ops.diff2_c_y)
								( ((ops.diff2_c_x+ops.diff2_c_y)(h_spec)).spectral_div_element_wise(ops.diff2_c_x+ops.diff2_c_y) ).toPhys()
							)
				).physical_reduce_rms_quad();

			std::cout << "SPEC: Error threshold for Laplace and its inverse: " << err3_laplace << std::endl;
			if (err3_laplace > eps)
			{
				std::cerr << "SPEC: Error threshold for Laplace too high for spectral differentiation!" << std::endl;
				exit(-1);
			}

			double err3_laplace_check =
				(
						h_cart-
							//((ops.diff_c_x.spectral_mul_element_wise(ops.diff_c_x)+ops.diff_c_y.spectral_mul_element_wise(ops.diff_c_y))(h_cart)).
							//spectral_div_element_wise(ops.diff2_c_x+ops.diff2_c_y)
							( ((ops.diff_c_x * ops.diff_c_x + ops.diff_c_y * ops.diff_c_y)(h_spec)) /
							(ops.diff2_c_x+ops.diff2_c_y) ).toPhys()
				).physical_reduce_rms_quad();

			std::cout << "Error for Laplace (diff*diff()) and its inverse (check): " << err3_laplace_check << std::endl;

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
			sweet::PlaneData_PhysicalComplex u(planeDataConfig);
			sweet::PlaneData_PhysicalComplex v(planeDataConfig);
			sweet::PlaneData_PhysicalComplex h_diff_x(planeDataConfig);
			sweet::PlaneData_PhysicalComplex h_diff_y(planeDataConfig);


//			Operators2D op(parameters.discretization.res, parameters.sim.domain_size, parameters.disc.use_spectral_diffs);

			for (int j = 0; j < shackPlaneDataOps->space_res_physical[1]; j++)
			{
				for (int i = 0; i < shackPlaneDataOps->space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)shackPlaneDataOps->space_res_physical[0];
					double y = ((double)j+0.5)/(double)shackPlaneDataOps->space_res_physical[1];
//					double x = ((double)i)/(double)parameters.discretization.res[0];
//					double y = ((double)j)/(double)parameters.discretization.res[1];

					u.physical_set_value(j, i, sin(freq_x*M_PI*x), 0);
					v.physical_set_value(j, i, cos(freq_y*M_PI*y), 0);

					h_cart.physical_set_value(
						j, i,
						sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y),
						0
					);

					h_diff_x.physical_set_value(
						j, i,
						freq_x*M_PI*cos(freq_x*M_PI*x)*cos(freq_y*M_PI*y)/(double)shackPlaneDataOps->plane_domain_size[0],
						0
					);

					h_diff_y.physical_set_value(
						j, i,
						-sin(freq_x*M_PI*x)*freq_y*M_PI*sin(freq_y*M_PI*y)/(double)shackPlaneDataOps->plane_domain_size[1],
						0
					);
				}
			}


			double res_normalization = std::sqrt(1.0/(shackPlaneDataOps->space_res_physical[0]*shackPlaneDataOps->space_res_physical[1]));

			// normalization for diff = 2 pi / L
			double err_x = (ops.diff_c_x(h_cart).toPhys()-h_diff_x).physical_reduce_norm2_quad()*res_normalization*shackPlaneDataOps->plane_domain_size[0]/(2.0*M_PI);
			double err_y = (ops.diff_c_y(h_cart).toPhys()-h_diff_y).physical_reduce_norm2_quad()*res_normalization*shackPlaneDataOps->plane_domain_size[1]/(2.0*M_PI);
			double err_z = (u*v-h_cart).physical_reduce_norm2_quad()*res_normalization;

			std::cout << "error diff x = " << err_x << std::endl;
			std::cout << "error diff y = " << err_y << std::endl;

			sweet::PlaneData_SpectralComplex h_spec(h_cart.planeDataConfig);
			h_spec.loadPlaneDataPhysical(h_cart);
			sweet::PlaneData_SpectralComplex h_diff_xy_spec  = ops.diff_c_x(h_spec) + ops.diff_c_y(h_spec);
			sweet::PlaneData_SpectralComplex h_diff_xy_spec_split  = (ops.diff_c_x + ops.diff_c_y)(h_spec);

			///PlaneDataComplex h_spec = h_cart;
			///PlaneDataComplex h_diff_xy_spec = ops.diff_c_x(h_spec) + ops.diff_c_y(h_spec);
			///PlaneDataComplex h_diff_xy_spec_split = (ops.diff_c_x + ops.diff_c_y)(h_spec);

			double err_xy = (
								h_diff_xy_spec
								-h_diff_x-h_diff_y
						).toPhys().physical_reduce_norm2_quad()*res_normalization/(2.0*M_PI);

			double err_xy_split = (
								h_diff_xy_spec_split
								-h_diff_x-h_diff_y
						).toPhys().physical_reduce_norm2_quad()*res_normalization/(2.0*M_PI);

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

			//double err_int_x = (h_cart-h_diff_x.spectral_div_element_wise(ops.diff_c_x)).reduce_norm2_quad()*res_normalization;
			sweet::PlaneData_SpectralComplex h_cart_spec(h_cart);
			sweet::PlaneData_SpectralComplex h_diff_x_spec(h_diff_x);
			double err_int_x = (h_cart_spec-h_diff_x_spec.spectral_div_element_wise(ops.diff_c_x)).toPhys().physical_reduce_norm2_quad()*res_normalization;
			std::cout << "Testing spectral inverse x " << err_int_x << std::endl;

			if (err_int_x > eps)
			{
				std::cerr << "SPEC: Error threshold for integration in x too high for spectral integration!" << std::endl;
				std::cout << err_int_x << std::endl;
				exit(-1);
			}

			//double err_int_y = (h_cart-h_diff_y.spectral_div_element_wise(ops.diff_c_y)).reduce_norm2_quad()*res_normalization;
			sweet::PlaneData_SpectralComplex h_diff_y_spec(h_diff_y);
			double err_int_y = (h_cart_spec-h_diff_y_spec.spectral_div_element_wise(ops.diff_c_y)).toPhys().physical_reduce_norm2_quad()*res_normalization;
			std::cout << "Testing spectral inverse y " << err_int_y << std::endl;

			if (err_int_y > eps)
			{
				std::cout << err_int_y << std::endl;
				std::cerr << "SPEC: Error threshold for integration in y too high for spectral integration!" << std::endl;
				exit(-1);
			}
		}

		std::cout << "TEST C: DONE" << std::endl;

	}


	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
