/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/expIntegration/ExpFunctions.hpp>


int main()
{
	double d_re = M_PI*1e-2;
	double d_im = M_PI*1e-5;
	typedef double T;
	typedef std::complex<T> CT;

	sweet::ExpFunctions<T> fun;
	fun.setup("phi0");

	CT x, y;

	CT prev_y = 0;

	int exit_value = EXIT_SUCCESS;
	for (int phi_or_ups_n = 0; phi_or_ups_n < 4; phi_or_ups_n++)
	{
		for (double x_re = -2; x_re < 2; x_re+=d_re)
		{
			prev_y = 0;
			for (double x_im = -2; x_im < 2; x_im+=d_im)
			{
				x = CT(x_re, x_im);

				double err;

				if (phi_or_ups_n < 4)
				{
					y = fun.phiN(phi_or_ups_n, x);
					err = std::abs(y - fun.phiN(phi_or_ups_n, x, 100));
				}
				else
				{
					y = fun.upsN(phi_or_ups_n-3, x);
					err = std::abs(y - fun.phiN(phi_or_ups_n-3, x, 100));
				}

				if (err > 1e-15)
				{
					std::cerr << "phi_or_ups " << phi_or_ups_n << " Error too high: " << err << std::endl;
					std::cerr << "x: " << x << std::endl;
					std::cerr << "y: " << y << std::endl;
					exit_value = EXIT_FAILURE;
					exit(1);
				}

				if (std::abs(prev_y) != 0)
				{
					if (std::abs(prev_y - y) > d_im*8)
					{
						std::cerr << "y delta too large for phi " << phi_or_ups_n << std::endl;
						std::cerr << "x: " << x << std::endl;
						std::cerr << "y: " << y << std::endl;
						std::cerr << "prev_y: " << prev_y << std::endl;
						exit_value = EXIT_FAILURE;
						exit(1);
					}
				}

				prev_y = y;
			}

			if (exit_value != EXIT_SUCCESS)
			{
				std::cerr << "FAILURE" << std::endl;
				return exit_value;
			}
		}

		std::cout << x << "\t" << y << std::endl;
	}

	return 0;

}
