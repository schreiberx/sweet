/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */



#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#include <ostream>
#include <cmath>


PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

SimulationVariables simVars;



int main(
		int i_argc,
		char *i_argv[]
)
{
	// override flag
	SimulationVariables simVars;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.space_res_physical[0];
	std::size_t res_y = simVars.disc.space_res_physical[1];

	std::cout << "*************************************************************" << std::endl;
	std::cout << "Testing aliasing pattern with resolution " << res_x << " x " << res_y << std::endl;
	std::cout << "*************************************************************" << std::endl;
	std::size_t res[2] = {res_x, res_y};

	simVars.disc.space_res_physical[0] = res[0];
	simVars.disc.space_res_physical[1] = res[1];
	simVars.reset();

	planeDataConfigInstance.setupAuto(
			simVars.disc.space_res_physical,
			simVars.disc.space_res_spectral,
			simVars.misc.reuse_spectral_transformation_plans
		);

	std::size_t test_max_freqx = planeDataConfig->spectral_real_modes[0];
	std::size_t test_max_freqy = planeDataConfig->spectral_real_modes[1];

	std::cout << "*************************************************************" << std::endl;
	planeDataConfig->printInformation();
	std::cout << "*************************************************************" << std::endl;
	std::cout << std::endl;

	PlaneData_Physical h1_x(planeDataConfig);
	PlaneData_Physical h2_x(planeDataConfig);

	PlaneData_Physical h1_y(planeDataConfig);
	PlaneData_Physical h2_y(planeDataConfig);

	PlaneData_Physical h1h2_numerical_x(planeDataConfig);
	PlaneData_Physical h1h2_numerical(planeDataConfig);
	PlaneData_Physical h1h2_numerical_y(planeDataConfig);

	PlaneData_Physical h1h2_analytical_x_high(planeDataConfig);
	PlaneData_Physical h1h2_analytical_x_low(planeDataConfig);
	PlaneData_Physical h1h2_analytical_x(planeDataConfig);


	PlaneData_Physical h1h2_analytical_y_high(planeDataConfig);
	PlaneData_Physical h1h2_analytical_y_low(planeDataConfig);
	PlaneData_Physical h1h2_analytical_y(planeDataConfig);

	PlaneData_Physical h1h2_analytical(planeDataConfig);

	double epsilon = 1e-10;

	std::size_t max_modes_x = planeDataConfig->spectral_real_modes[0];
	std::size_t max_modes_y = planeDataConfig->spectral_real_modes[1];

	/*
	 * Use cos(x)*cos(x) instead of sin*sin?
	 */
//	int coscos_y = 1;
//	int coscos_x = 1;
	for (int coscos_y = 0; coscos_y < 2; coscos_y++)
	for (int coscos_x = 0; coscos_x < 2; coscos_x++)
	{
		//std::size_t freq1_y = 0;
		for (std::size_t freq1_y = 0; freq1_y < test_max_freqy; freq1_y++)
		{
			/*
			 * Setup frequency
			 * 	sin(fx1*pi*x)
			 */
			h1_y.physical_update_lambda_unit_coordinates_corner_centered(
				[&](double x, double y, double &io_data)
				{
					if (coscos_y)
						io_data = std::cos((double)freq1_y*2.0*M_PI*y);
					else
						io_data = std::sin((double)freq1_y*2.0*M_PI*y);
				}
			);


			//std::size_t freq2_y = 0;
			for (std::size_t freq2_y = 0; freq2_y < test_max_freqy; freq2_y++)
			{
				h2_y.physical_update_lambda_unit_coordinates_corner_centered(
					[&](double x, double y, double &io_data)
					{
						if (coscos_y)
							io_data = std::cos((double)freq2_y*2.0*M_PI*y);
						else
							io_data = std::sin((double)freq2_y*2.0*M_PI*y);
					}
				);

				h1h2_numerical_y = h1_y*h2_y;

				/*
				 * Compute analytical solution
				 *
				 * Setup lower frequency
				 * \frac{1}{2}\cos\left(\left(k_{1}-k_{2}\right)x\pi\right)
				 */
				h1h2_analytical_y_low.physical_update_lambda_unit_coordinates_corner_centered(
						[&](double x, double y, double &io_data)
						{
							io_data = 0.5*std::cos((double)((int)freq1_y-(int)freq2_y)*2.0*M_PI*y);
						}
					);

				if (freq1_y+freq2_y < max_modes_y)
				{
					/*
					 * Setup higher frequency
					 * -\frac{1}{2}\cos\left(\left(k_{1}+k_{2}\right)y\pi\right)
					 */
					h1h2_analytical_y_high.physical_update_lambda_unit_coordinates_corner_centered(
							[&](double x, double y, double &io_data)
							{
								io_data = -0.5*std::cos((double)((int)freq1_y+(int)freq2_y)*2.0*M_PI*y);

								if (coscos_y)
									io_data = -io_data;
							}
						);
				}
				else
				{
					h1h2_analytical_y_high.physical_set_all_value(0.);
					std::cout << "Higher mode truncated in this case!" << std::endl;
				}


				/*
				 * Merge solutions
				 */
				h1h2_analytical_y = h1h2_analytical_y_low + h1h2_analytical_y_high;

				for (std::size_t freq1_x = 0; freq1_x < test_max_freqx; freq1_x++)
				{
					/*
					 * Setup frequency
					 * 	sin(fx1*pi*x)
					 */
					h1_x.physical_update_lambda_unit_coordinates_corner_centered(
						[&](double x, double y, double &io_data)
						{
							if (coscos_x)
								io_data = std::cos((double)freq1_x*2.0*M_PI*x);
							else
								io_data = std::sin((double)freq1_x*2.0*M_PI*x);
						}
					);


					std::size_t freq2_x = freq1_x;
					//for (std::size_t freq_x2 = 0; freq_x2 < test_max_freqx; freq_x2++)
					{
						std::cout << "***********************************************************" << std::endl;
						std::cout << "Resolution physical space: (" << planeDataConfig->physical_res[0] << ", " << planeDataConfig->physical_res[1] << ")" << std::endl;
						std::cout << "           spectral space: (" << planeDataConfig->spectral_real_modes[0] << ", " << planeDataConfig->spectral_real_modes[1] << ")" << std::endl;
						std::cout << "Testing for frequency fx1=" << freq1_x << " and fx2=" << freq2_x << ", max frequency x=" << test_max_freqx << std::endl;
						std::cout << "                      fy1=" << freq1_y << " and fy2=" << freq2_y << ", max frequency y=" << test_max_freqy << std::endl;
						std::cout << "                      Test functions: " << (coscos_x ? "coscos_x" : "sinsin_x") << ", " << (coscos_y ? "coscos_y" : "sinsin_y") << std::endl;
						/*
						 * Setup frequency
						 * 	sin(fx1*pi*x)
						 */
						h2_x.physical_update_lambda_unit_coordinates_corner_centered(
							[&](double x, double y, double &io_data)
							{
								if (coscos_x)
									io_data = std::cos((double)freq2_x*2.0*M_PI*x);
								else
									io_data = std::sin((double)freq2_x*2.0*M_PI*x);
							}
						);


						/*
						 * NOTE
						 *
						 * See also
						 * ./doc/software_development_discussions/antialiasing/antialiasing_rule_and_tests.pdf
						 */
						h1h2_numerical_x = h1_x*h2_x;

						/*
						 * Setup lower frequency
						 * \frac{1}{2}\cos\left(\left(k_{1}-k_{2}\right)x\pi\right)
						 */
						h1h2_analytical_x_low.physical_update_lambda_unit_coordinates_cell_centered(
								[&](double x, double y, double &io_data)
								{
									io_data = 0.5*std::cos((double)((int)freq1_x-(int)freq2_x)*2.0*M_PI*x);
								}
							);


						if (freq1_x+freq2_x < max_modes_x)
						{
							/*
							 * Setup higher frequency
							 * -\frac{1}{2}\cos\left(\left(k_{1}+k_{2}\right)x\pi\right)
							 */
							h1h2_analytical_x_high.physical_update_lambda_unit_coordinates_corner_centered(
									[&](double x, double y, double &io_data)
									{
										io_data = -0.5*std::cos((double)((int)freq1_x+(int)freq2_x)*2.0*M_PI*x);

										if (coscos_x)
											io_data = -io_data;
									}
								);
						}
						else
						{
							h1h2_analytical_x_high.physical_set_all_value(0.);
							std::cout << "Higher mode truncated in this case!" << std::endl;
						}


						h1h2_analytical_x = h1h2_analytical_x_low + h1h2_analytical_x_high;

						{
							h1h2_numerical = h1h2_numerical_x + h1h2_numerical_y;
							h1h2_analytical = h1h2_analytical_x + h1h2_analytical_y;

							double max_error = (h1h2_analytical - h1h2_numerical).physical_reduce_max_abs();
							std::cout << "Max error: " << max_error << std::endl;

							if (max_error > epsilon)
							{
								PlaneData_Spectral tmp_spec(h1_y.planeDataConfig);

								std::cout << "H1_y:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1_y);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "H2_y:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h2_y);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "H1_x:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1_x);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "H2_x:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h2_x);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Numerical spectrum:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1h2_numerical);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum x:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1h2_analytical_x);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum y (low):" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1h2_analytical_y_low);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum y (high):" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1h2_analytical_y_high);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum y:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1h2_analytical_y);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Analytical spectrum:" << std::endl;
								tmp_spec.loadPlaneDataPhysical(h1h2_analytical);
								tmp_spec.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								if (max_error > epsilon)
									SWEETError("Error too high!");
							}
						}
					}
				}
			}
		}
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
