#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
	#error "Aliasing activated! Deactivate it to use this unit test"
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
	int res_x = io_data.planeDataConfig->spectral_modes[0];
	int res_y = io_data.planeDataConfig->spectral_modes[1];

	/*
	 * Phase shift doesn't work with the highest representable frequency
	 */
//	double phase_shift = M_PI*0.3;
	double phase_shift = 0;

	io_data.physical_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				if (fx < 0 || fy < 0)
				{
					o_data = 0.0;
						return;
				}

				o_data = 1.0;

				if (fx >= 0)
					o_data *= std::cos(((double)x)*(double)fx*M_PI*2.0/(double)res_x + phase_shift);

				if (fy >= 0)
					o_data *= std::cos(((double)y)*(double)fy*M_PI*2.0/(double)res_y + phase_shift);
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

	if (simVars.disc.res_spectral[0] <= 0)
	{
		simVars.disc.res_spectral[0] = simVars.disc.res_physical[0];
		simVars.disc.res_spectral[1] = simVars.disc.res_physical[1];
	}
	else
	{
		simVars.disc.res_physical[0] = simVars.disc.res_spectral[0];
		simVars.disc.res_physical[1] = simVars.disc.res_spectral[1];
	}


	int max_res = 64;

	double epsilon = 1e-12;

	for (int res = simVars.disc.res_spectral[0]; res <= max_res; res+=2)
	{
		simVars.disc.res_physical[0] = res;
		simVars.disc.res_physical[1] = res;
		simVars.disc.res_spectral[0] = res;
		simVars.disc.res_spectral[1] = res;
		simVars.reset();

		/**
		 * Here we enforce the same physical and spectral resolution
		 */
		planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

		std::cout << "PHYS RES: " << planeDataConfig->physical_res[0] << " x " << planeDataConfig->physical_res[1] << std::endl;
		std::cout << "PHYS SIZE: " << planeDataConfig->physical_data_size[0] << " x " << planeDataConfig->physical_data_size[1] << std::endl;
		std::cout << "SPEC MODES: " << planeDataConfig->spectral_modes[0] << " x " << planeDataConfig->spectral_modes[1] << std::endl;
		std::cout << "SPEC SIZE: " << planeDataConfig->spectral_data_size[0] << " x " << planeDataConfig->spectral_data_size[1] << std::endl;
		std::cout << std::endl;

//		if (res_x != res_y)
//			FatalError("Only same resolutions in x and y dimension are currently supported!");

		if ((res & 1) == 1)
			FatalError("Only even resolutions supported");

		PlaneData a(planeDataConfig);

		// relative spectral resolution deltas to test restriction / interpolation
		//for (int spec_delta = -4; spec_delta <= 4; spec_delta+=2)
		for (int spec_delta = -4; spec_delta <= 4; spec_delta+=2)
		{
			int dst_res = (int)simVars.disc.res_physical[0]+spec_delta;

			PlaneDataConfig planeDataConfigInstanceDst;
			planeDataConfigInstanceDst.setup(
					dst_res, dst_res,
					dst_res, dst_res
				);

			PlaneDataConfig *planeDataConfigDst = &planeDataConfigInstanceDst;

			/*
			 * Iterate over relative frequencies
			 * Only iterate up to res_x/2 frequency since higher frequencies would be only
			 * representable as aliased ones.
			 */
			for (int freq_x = 0; freq_x <= res/2; freq_x += 1)
			{
				for (int freq_y = 0; freq_y <= res/2; freq_y += 1)
//				int freq_y = 2;
				{
					std::cout << "*************************************************************" << std::endl;
					std::cout << "Testing (" << res << ", " << res << ") -> (" << planeDataConfigDst->spectral_modes[0] << ", " << planeDataConfigDst->spectral_modes[1] << ")"  << std::endl;
					std::cout << "freq: " << freq_x << ", " << freq_y << std::endl;
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

						std::cout << "A (spectral):" << std::endl;
						a.print_spectralData_zeroNumZero();
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
								FatalError("Test for conserving high frequencies failed! Results should be identical!");
							}
						}

						PlaneData b = a.spectral_returnWithDifferentModes(planeDataConfigDst);

						std::cout << "Frequency in original resolution: " << freq_x << ", " << freq_y << std::endl;
						int test_freq_x = freq_x;
						if (dst_res/2 < freq_x)
							test_freq_x = -1;

						int test_freq_y = freq_y;
						if (dst_res/2 < freq_y)
							test_freq_y = -1;

						{
							std::cout << "Test for existing frequencies " << test_freq_x << ", " << test_freq_y << std::endl;

							PlaneData test(planeDataConfigDst);
							setupDataFreq(test, test_freq_x, test_freq_y);
#if 0
							std::cout << "Spectral data of b:" << std::endl;
							b.print_spectralData_zeroNumZero();

							std::cout << "Spectral data of test:" << std::endl;
							test.print_spectralData_zeroNumZero();
#endif
							error = (b-test).reduce_maxAbs();
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
								FatalError("No modes changed! Results should be identical!");
							}

							std::cout << "PASSED (freq test) with error of " << error << std::endl;
						}
					}

					if (spec_delta >= 0)
					{
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
							FatalError("Mode extension requested, but interpolation followed by restriction does not return identical results!");
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
