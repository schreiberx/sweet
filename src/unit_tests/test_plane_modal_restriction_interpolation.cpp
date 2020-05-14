#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
//	#error "Aliasing activated! Deactivate it to use this unit test"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif



#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>

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


void setupDataFreq(
		PlaneData &io_data,
		int fx,	///< frequency x
		int fy	///< frequency y
)
{
	int res_x = io_data.planeDataConfig->physical_res[0];
	int res_y = io_data.planeDataConfig->physical_res[1];

	// shift by half a cell to generate exactly this mode in spectral space
	double phase_shift = 0.0;

	io_data.physical_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = 0.0;

				if (fx >= 0)
					o_data += std::cos(((double)x+phase_shift)*(double)fx*M_PI*2.0/(double)res_x);

				if (fy >= 0)
					o_data += std::cos(((double)y+phase_shift)*(double)fy*M_PI*2.0/(double)res_y);
			}
	);
}


void setupData123(
		PlaneData &io_data
)
{

	io_data.physical_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = (x+1.0)+(y+3.0)*y;
			}
	);
}

int main(int i_argc, char *i_argv[])
{
	// override flag
	SimulationVariables simVars;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */

	if (simVars.disc.space_res_spectral[0] <= 0)
	{
		simVars.disc.space_res_spectral[0] = simVars.disc.space_res_physical[0];
		simVars.disc.space_res_spectral[1] = simVars.disc.space_res_physical[1];
	}
	else
	{
		simVars.disc.space_res_physical[0] = simVars.disc.space_res_spectral[0];
		simVars.disc.space_res_physical[1] = simVars.disc.space_res_spectral[1];
	}


	int max_res = 64;

	double epsilon = 1e-12;

	for (int res = simVars.disc.space_res_spectral[0]; res <= max_res; res+=2)
	{
		simVars.disc.space_res_physical[0] = res;
		simVars.disc.space_res_physical[1] = res;
		simVars.disc.space_res_spectral[0] = 0;
		simVars.disc.space_res_spectral[1] = 0;
		simVars.reset();

		/**
		 * Here we enforce the same physical and spectral resolution
		 */
		planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);
		planeDataConfigInstance.printInformation();
		std::cout << std::endl;

		if ((res & 1) == 1)
			SWEETError("Only even resolutions supported");

		PlaneData a(planeDataConfig);

		// relative spectral resolution deltas to test restriction / interpolation
		for (int spec_delta = -4; spec_delta <= 4; spec_delta+=2)
		{
			int dst_res_physical[2] = {
					(int)simVars.disc.space_res_physical[0]+spec_delta,
					(int)simVars.disc.space_res_physical[1]+spec_delta
			};

			if (dst_res_physical[0] < 4 || dst_res_physical[1] < 4)
				continue;

			int dst_res_spectral[2] = {0, 0};

			PlaneDataConfig planeDataConfigInstanceDst;
			planeDataConfigInstanceDst.setupAuto(dst_res_physical, dst_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

			PlaneDataConfig *planeDataConfigDst = &planeDataConfigInstanceDst;

			/*
			 * Iterate over relative frequencies
			 * Only iterate up to real_modes/2 frequency since higher frequencies would be only
			 * representable as aliased ones.
			 */
			for (int freq_x = 0; freq_x < (int)planeDataConfig->spectral_real_modes[0]; freq_x += 1)
			{
				for (int freq_y = 0; freq_y < (int)planeDataConfig->spectral_real_modes[1]; freq_y += 1)
				{
					std::cout << std::endl;
					std::cout << std::endl;
					std::cout << "*************************************************************" << std::endl;
					std::cout << "Testing (" << res << ", " << res << ") -> (" << planeDataConfigDst->spectral_modes[0] << ", " << planeDataConfigDst->spectral_modes[1] << ")"  << std::endl;
					std::cout << " + testing frequency in source data: " << freq_x << ", " << freq_y << std::endl;

					int test_freq_x = freq_x;
					if (planeDataConfigDst->spectral_real_modes[0] <= (std::size_t)freq_x)
						test_freq_x = -1;

					int test_freq_y = freq_y;
					if (planeDataConfigDst->spectral_real_modes[1] <= (std::size_t)freq_y)
						test_freq_y = -1;

					std::cout << " + testing frequency in destination data: " << test_freq_x << ", " << test_freq_y << std::endl;
					std::cout << "*************************************************************" << std::endl;

					double error = 0;

					/*
					 * Setup data with highest possible frequency
					 */

					{
						setupDataFreq(a, freq_x, freq_y);

#if 0
						std::cout << "A (physical):" << std::endl;
						a.print_physicalArrayData();
						std::cout << std::endl;

						std::cout << "A (spectral):" << std::endl;
						a.print_spectralData_zeroNumZero();
						std::cout << std::endl;
#endif

						{
							/*
							 * Test for conserving high frequencies in Fourier transformations
							 */
							PlaneData tmp = a;
							tmp.request_data_spectral();
							tmp.request_data_physical();

							error = (a-tmp).reduce_maxAbs();
							if (error > epsilon)
							{
								std::cout << "Error: " << error << std::endl;
								SWEETError("Test for conserving high frequencies failed! Results should be identical!");
							}
						}

						PlaneData b = a.spectral_returnWithDifferentModes(planeDataConfigDst);

						{
							PlaneData test(planeDataConfigDst);
							setupDataFreq(test, test_freq_x, test_freq_y);

#if 0
							std::cout << "Spectral data of b:" << std::endl;
							b.print_spectralData_zeroNumZero();
							std::cout << std::endl;

							std::cout << "Spectral data of test:" << std::endl;
							test.print_spectralData_zeroNumZero();
							std::cout << std::endl;
#endif

							error = (b-test).reduce_maxAbs();

							if (error > epsilon)
							{
								a.print_spectralData_zeroNumZero();
								std::cout << "**************************************************" << std::endl;
								std::cout << "* ERROR" << std::endl;
								std::cout << "* a = freq(fx,fy)" << std::endl;
								std::cout << "* b = interp_restrict(a)" << std::endl;
								std::cout << "* test = expected_interp_restrict(fx,fy)" << std::endl;
								std::cout << "**************************************************" << std::endl;
								std::cout << "Spectral data of A:" << std::endl;
								a.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Physical data of A:" << std::endl;
								a.print_physicalData_zeroNumZero();
								std::cout << std::endl;


								std::cout << "Spectral data of b:" << std::endl;
								b.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Physical data of b:" << std::endl;
								b.print_physicalData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Spectral data of test:" << std::endl;
								test.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Physical data of test:" << std::endl;
								test.print_physicalData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Error: " << error << std::endl;
								SWEETError("No modes changed! Results should be identical!");
							}

							std::cout << "PASSED (freq test) with error of " << error << std::endl;
						}
					}

					if (spec_delta >= 0)
					{
						std::cout << "TESTING for conservation of modes" << std::endl;

						setupData123(a);

						PlaneData b = a.spectral_returnWithDifferentModes(planeDataConfigDst);
						PlaneData test = b.spectral_returnWithDifferentModes(planeDataConfig);

						error = (a-test).reduce_maxAbs();
						if (error > epsilon)
						{
							std::cout << "**************************************************" << std::endl;
							std::cout << "* ERROR" << std::endl;
							std::cout << "**************************************************" << std::endl;
							std::cout << "Spectral data of A:" << std::endl;
							a.print_spectralData_zeroNumZero();

							std::cout << "Spectral data of b:" << std::endl;
							b.print_spectralData_zeroNumZero();

							std::cout << "Spectral data of test:" << std::endl;
							test.print_spectralData_zeroNumZero();

							std::cout << "Error: " << error << std::endl;
							SWEETError("Mode extension requested, but interpolation followed by restriction does not return identical results!");
						}

						std::cout << "PASSED (setup 123) with error of " << error << std::endl;
					}
				}
			}
		}
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
