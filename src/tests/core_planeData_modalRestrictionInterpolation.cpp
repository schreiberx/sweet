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




void setupDataFreq(
		sweet::PlaneData_Spectral &io_data,
		int fx,	///< frequency x
		int fy	///< frequency y
)
{
	int res_x = io_data.planeDataConfig->physical_res[0];
	int res_y = io_data.planeDataConfig->physical_res[1];

	// shift by half a cell to generate exactly this mode in spectral space
	double phase_shift = 0.0;

	sweet::PlaneData_Physical tmp(io_data.planeDataConfig);

	tmp.physical_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = 0.0;

				if (fx >= 0)
					o_data += std::cos(((double)x+phase_shift)*(double)fx*M_PI*2.0/(double)res_x);

				if (fy >= 0)
					o_data += std::cos(((double)y+phase_shift)*(double)fy*M_PI*2.0/(double)res_y);
			}
	);

	io_data.loadPlaneDataPhysical(tmp);
}


void setupData123(
		sweet::PlaneData_Spectral &io_data
)
{

	sweet::PlaneData_Physical tmp(io_data.planeDataConfig);

	tmp.physical_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = (x+1.0)+(y+3.0)*y;
			}
	);

	io_data.loadPlaneDataPhysical(tmp);
}

int main(int i_argc, char *i_argv[])
{
	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackPlaneDataOps *shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.printShackData();

	if (shackPlaneDataOps->space_res_spectral[0] <= 0)
	{
		shackPlaneDataOps->space_res_spectral[0] = shackPlaneDataOps->space_res_physical[0];
		shackPlaneDataOps->space_res_spectral[1] = shackPlaneDataOps->space_res_physical[1];
	}
	else
	{
		shackPlaneDataOps->space_res_physical[0] = shackPlaneDataOps->space_res_spectral[0];
		shackPlaneDataOps->space_res_physical[1] = shackPlaneDataOps->space_res_spectral[1];
	}


	int max_res = 64;

	double epsilon = 1e-12;

	for (int res = shackPlaneDataOps->space_res_spectral[0]; res <= max_res; res+=2)
	{
		shackPlaneDataOps->space_res_physical[0] = res;
		shackPlaneDataOps->space_res_physical[1] = res;
		shackPlaneDataOps->space_res_spectral[0] = 0;
		shackPlaneDataOps->space_res_spectral[1] = 0;

		/**
		 * Here we enforce the same physical and spectral resolution
		 */
		sweet::PlaneData_Config planeDataConfig;
		planeDataConfig.setupAuto(shackPlaneDataOps);
		planeDataConfig.printInformation();
		std::cout << std::endl;

		if ((res & 1) == 1)
			SWEETError("Only even resolutions supported");

		sweet::PlaneData_Spectral a(planeDataConfig);

		// relative spectral resolution deltas to test restriction / interpolation
		for (int spec_delta = -4; spec_delta <= 4; spec_delta+=2)
		{
			int dst_res_physical[2] = {
					(int)shackPlaneDataOps->space_res_physical[0]+spec_delta,
					(int)shackPlaneDataOps->space_res_physical[1]+spec_delta
			};

			if (dst_res_physical[0] < 4 || dst_res_physical[1] < 4)
				continue;

			int dst_res_spectral[2] = {0, 0};

			sweet::PlaneData_Config planeDataConfigDst;
			planeDataConfigDst.setupAuto(dst_res_physical, dst_res_spectral, shackPlaneDataOps->reuse_spectral_transformation_plans);

			/*
			 * Iterate over relative frequencies
			 * Only iterate up to real_modes/2 frequency since higher frequencies would be only
			 * representable as aliased ones.
			 */
			for (int freq_x = 0; freq_x < (int)planeDataConfig.spectral_real_modes[0]; freq_x += 1)
			{
				for (int freq_y = 0; freq_y < (int)planeDataConfig.spectral_real_modes[1]; freq_y += 1)
				{
					std::cout << std::endl;
					std::cout << std::endl;
					std::cout << "*************************************************************" << std::endl;
					std::cout << "Testing (" << res << ", " << res << ") -> (" << planeDataConfigDst.spectral_modes[0] << ", " << planeDataConfigDst.spectral_modes[1] << ")"  << std::endl;
					std::cout << " + testing frequency in source data: " << freq_x << ", " << freq_y << std::endl;

					int test_freq_x = freq_x;
					if (planeDataConfigDst.spectral_real_modes[0] <= (std::size_t)freq_x)
						test_freq_x = -1;

					int test_freq_y = freq_y;
					if (planeDataConfigDst.spectral_real_modes[1] <= (std::size_t)freq_y)
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
							sweet::PlaneData_Spectral tmp = a;

							error = (a-tmp).spectral_reduce_max_abs();
							if (error > epsilon)
							{
								std::cout << "Error: " << error << std::endl;
								SWEETError("Test for conserving high frequencies failed! Results should be identical!");
							}
						}

						sweet::PlaneData_Spectral b = a.spectral_returnWithDifferentModes(planeDataConfigDst);

						{
							sweet::PlaneData_Spectral test(planeDataConfigDst);
							setupDataFreq(test, test_freq_x, test_freq_y);

#if 0
							std::cout << "Spectral data of b:" << std::endl;
							b.print_spectralData_zeroNumZero();
							std::cout << std::endl;

							std::cout << "Spectral data of test:" << std::endl;
							test.print_spectralData_zeroNumZero();
							std::cout << std::endl;
#endif

							error = (b-test).spectral_reduce_max_abs();

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
								a.toPhys().print_physicalData_zeroNumZero();
								std::cout << std::endl;


								std::cout << "Spectral data of b:" << std::endl;
								b.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Physical data of b:" << std::endl;
								b.toPhys().print_physicalData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Spectral data of test:" << std::endl;
								test.print_spectralData_zeroNumZero();
								std::cout << std::endl;

								std::cout << "Physical data of test:" << std::endl;
								test.toPhys().print_physicalData_zeroNumZero();
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

						sweet::PlaneData_Spectral b = a.spectral_returnWithDifferentModes(planeDataConfigDst);
						sweet::PlaneData_Spectral test = b.spectral_returnWithDifferentModes(planeDataConfig);

						error = (a-test).spectral_reduce_max_abs();
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
