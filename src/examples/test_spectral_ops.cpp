
//#if !SWEET_USE_SPECTRAL_SPACE
//	#error "Spectral space not activated"
//#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif



#include <sweet/DataArray.hpp>
#include <sweet/SimulationParameters.hpp>
#include <sweet/Operators2D.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>

SimulationParameters parameters;


int main(int i_argc, char *i_argv[])
{
	// error tolerance
	double eps = 1e-10;
	double eps_convergence = 1e-4;

	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);

	SimulationParameters parameters;
	parameters.use_spectral_diffs = 1;
	parameters.setup(i_argc, i_argv);

	if (parameters.use_spectral_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	double prev_error_x = 0;
	double prev_error_y = 0;

	double freq_x = 2.0;
	double freq_y = 4.0;


	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = parameters.res[0];
	std::size_t res_y = parameters.res[1];

	for (; res_x <= 4096 && res_y <= 2048; res_x *= 2, res_y *= 2)
	{
		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		parameters.res[0] = res[0];
		parameters.res[1] = res[1];
		parameters.reset();

		DataArray<2> h(res);
		DataArray<2> u(res);
		DataArray<2> v(res);
		DataArray<2> h_diff_x(res);
		DataArray<2> h_diff_y(res);

		{
			DataArray<2> zero(res);
			DataArray<2> two(res);
			DataArray<2> five(res);
			DataArray<2> h(res);

			zero.setAll(0);
			two.setAll(2);
			five.setAll(5);

			double add_test_two = (zero+two).reduce_sumAbs_quad()/(double)(res[0]*res[1]);
			double add_test_seven = (five+two).reduce_sumAbs_quad()/(double)(res[0]*res[1]);
			double add_test_ten = ((five+two)+3.0).reduce_sumAbs_quad()/(double)(res[0]*res[1]);
			double error = 0;

			std::cout << "Add test two: " << add_test_two << std::endl;
			error = std::abs(add_test_two-2.0);
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			std::cout << "Add test seven: " << add_test_seven << std::endl;
			error = std::abs(add_test_seven-7.0);
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			std::cout << "Add test ten: " << add_test_seven << std::endl;
			error = std::abs(add_test_ten-10.0);
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

			// has to be zero
			double sin_test_zero = h.reduce_sum_quad()/(double)(res[0]*res[1]);
			std::cout << "Sin test zero: " << sin_test_zero << std::endl;
			error = std::abs(sin_test_zero);
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}

			double sin_test_six = (h+6.0).reduce_sum_quad()/(double)(res[0]*res[1]);
			std::cout << "Sin test add six: " << sin_test_six << std::endl;
			error = std::abs(sin_test_six-6.0);
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}


			double sin_test_zero_mul = (h*two).reduce_sum_quad()/(double)(res[0]*res[1]);
			std::cout << "Sin test add zero: " << sin_test_zero_mul << std::endl;
			error = std::abs(sin_test_zero_mul);
			if (error > eps)
			{
				std::cout << "FAILED with error " << error;
				exit(-1);
			}
		}

		Operators2D op(parameters.res, parameters.sim_domain_length, parameters.use_spectral_diffs);

//		double scale = ((double)parameters.sim_domain_length[0]*(double)parameters.sim_domain_length[1]);
		double scale = 1.0;

		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				double x = ((double)i)/(double)parameters.res[0];
				double y = ((double)j)/(double)parameters.res[1];

#define FUN_ID	1

#if FUN_ID==1
				u.set(j, i, sin(freq_x*M_PIl*x));
				v.set(j, i, cos(freq_y*M_PIl*y)*scale);
#elif FUN_ID==2
				u.set(j, i, sin(freq_x*M_PIl*x));
				v.set(j, i, 1.0/(cos(freq_y*M_PIl*y)+2.0)*scale);
#endif

				h.set(
					j, i,
#if FUN_ID==1
					sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*scale
#elif FUN_ID==2
					sin(freq_x*M_PIl*x)*sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*cos(freq_y*M_PIl*y)*scale
#elif FUN_ID==3
					sin(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)*scale
#endif
				);

				h_diff_x.set(
					j, i,
#if FUN_ID==1
					freq_x*M_PIl*cos(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(double)parameters.sim_domain_length[0]*scale
#elif FUN_ID==2
					2.0*sin(freq_x*M_PIl*x)*std::pow(cos(freq_y*M_PIl*y),2.0)*freq_x*M_PIl*cos(freq_x*M_PIl*x)/(double)parameters.sim_domain_length[0]*scale
#elif FUN_ID==3
					freq_x*M_PIl*cos(freq_x*M_PIl*x)/(cos(freq_y*M_PIl*y)+2.0)/(double)parameters.sim_domain_length[0]*scale
#endif
				);

				h_diff_y.set(
					j, i,
#if FUN_ID==1
					-sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/(double)parameters.sim_domain_length[1]*scale
#elif FUN_ID==2
					-2.0*std::pow(std::sin(freq_x*M_PIl*x),2.0)*std::cos(freq_y*M_PIl*y)*freq_y*M_PIl*std::sin(freq_y*M_PIl*y)/(double)parameters.sim_domain_length[1]*scale
#elif FUN_ID==3
					sin(freq_x*M_PIl*x)*freq_y*M_PIl*sin(freq_y*M_PIl*y)/pow(cos(freq_y*M_PIl*y)+2.0, 2.0)/(double)parameters.sim_domain_length[1]*scale
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

		double normalization = 1.0/((double)res[0]*(double)res[1]);

		double err_x = (op.diff_c_x(h)-h_diff_x).reduce_sumAbs_quad()*normalization;
		double err_y = (op.diff_c_y(h)-h_diff_y).reduce_sumAbs_quad()*normalization;
		double err_z = (u*v-h).reduce_sumAbs_quad()*normalization;

		if (parameters.use_spectral_diffs)
		{
			std::cout << "(" << res_x << "x" << res_y << ")\terr_x=" << err_x << "\terr_y=" << err_y << "\terr_z=" << err_z << std::endl;

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
			if (abs(err_z) > eps)
			{
				std::cerr << "SPEC: Error threshold exceeded for err_z!" << std::endl;
				exit(-1);
			}
#endif
		}
		else
		{
			double conv_x = prev_error_x/err_x;
			double conv_y = prev_error_y/err_y;
			std::cout <<
					"(" << res_x << "x" << res_y << ")" <<
					"\terr_x=" << err_x << "\terr_y=" << err_y << "\terr_z=" << err_z <<
					"\tconv_x=" << conv_x << "\tconv_y=" << conv_y  <<
					std::endl;


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

			prev_error_x = err_x;
			prev_error_y = err_y;
		}
	}


	return 0;
}
