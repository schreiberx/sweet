
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

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
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);
	std::cout << std::setprecision(4);
	std::cerr << std::setprecision(4);

	SimulationParameters parameters;
	parameters.use_spectral_diffs = 1;
	parameters.setup(i_argc, i_argv);

	if (parameters.use_spectral_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	double prev_error_x = 0;
	double prev_error_y = 0;

	double freq_x = 4.0;
	double freq_y = 4.0;

	for (std::size_t res_x = 32; res_x <= 4096; res_x *= 2)
//	for (std::size_t res_y = 32; res_y <= 4096; res_y *= 2)
	{
		std::size_t res_y = res_x;
		std::size_t res[2] = {res_x,res_y};

		parameters.res[0] = res[0];
		parameters.res[1] = res[1];
		parameters.reset();

		DataArray<2> h(res);
		DataArray<2> h_diff_x(res);
		DataArray<2> h_diff_y(res);

		Operators2D op(parameters.res, parameters.sim_domain_length, parameters.use_spectral_diffs);

		// scale factor to avoid very very tiny values
		double scale = ((double)parameters.sim_domain_length[0]*(double)parameters.sim_domain_length[1]);
//		scale *= scale;
//		double scale = 1.0;

		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				double x = ((double)i)/(double)parameters.res[0];
				double y = ((double)j)/(double)parameters.res[1];

				h.set(
					j, i,
					sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*scale
				);

				h_diff_x.set(
					j, i,
					freq_x*M_PIl*cos(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(double)parameters.sim_domain_length[0]*scale
				);

				h_diff_y.set(
					j, i,
					sin(freq_x*M_PIl*x)*freq_y*M_PIl*(-sin(freq_y*M_PIl*y))/(double)parameters.sim_domain_length[1]*scale
				);
			}
		}

		double normalization = 1.0/((double)res[0]*(double)res[1]);
		double err_x = (op.diff_c_x(h)-h_diff_x).reduce_sumAbs()*normalization;
		double err_y = (op.diff_c_y(h)-h_diff_y).reduce_sumAbs()*normalization;

		if (parameters.use_spectral_diffs)
		{
			std::cout << res_x << "\t" << res_y << "\t" << err_x << "\t" << err_y << std::endl;

			if (err_x > 10e-10)
				std::cerr << "Error threshold for diff-X too high for spectral differentiation!" << std::endl;

			if (err_y > 10e-10)
				std::cerr << "Error threshold for diff-Y too high for spectral differentiation!" << std::endl;
		}
		else
		{
			std::cout << res_x << "\t" << res_y << "\t" << err_x << "\t" << err_y << "\t" << prev_error_x/err_x << "\t" << prev_error_y/err_y << std::endl;

			prev_error_x = err_x;
			prev_error_y = err_y;
		}
	}


	return 1;
}
