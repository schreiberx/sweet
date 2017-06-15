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


void setupData(
		PlaneData &io_data,
		int fx,	///< frequency x
		int fy	///< frequency y
)
{
	int res_x = io_data.planeDataConfig->spectral_modes[0];
//	int res_y = io_data.planeDataConfig->spectral_modes[1];

	/*
	 * Phase shift doesn't work with the highest representable frequency
	 */
//	double phase_shift = M_PI*0.3;
	double phase_shift = 0;

	io_data.physical_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = std::cos(((double)x)*(double)fx*M_PI*2.0/(double)res_x + phase_shift);
			}
	);

#if 0
	a.physical_update_lambda_array_indices(
			[&](int x, int y, double &o_data)
			{
				o_data = (x+1)*(y*y+1);
			}
	);
#endif
#if 0
	a.spectral_update_lambda_array_indices(
			[&](int idx, std::complex<double> &o_data)
			{
				o_data = idx;
			}
	);
#endif
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


	int max_res = 128;

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
		//	for (int spec_delta = -4; spec_delta <= 4; spec_delta+=2)
		int spec_delta = 2;
		{
			double dst_res = simVars.disc.res_physical[0]+spec_delta;

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
			for (int freq_delta = -(int)res; freq_delta <= -res/2; freq_delta += 1)
			{
				int freq = res+freq_delta;

				std::cout << "*************************************************************" << std::endl;
				std::cout << "Testing (" << res << ", " << res << ") -> (" << planeDataConfigDst->spectral_modes[0] << ", " << planeDataConfigDst->spectral_modes[1] << ")"  << std::endl;
				std::cout << "freq: " << freq << std::endl;
				std::cout << "*************************************************************" << std::endl;

				double error = 0;

				/*
				 * Setup data with highest possible frequency
				 */

				setupData(a, freq, freq);

				std::cout << "A (physical):" << std::endl;
				a.print_physicalArrayData();

				std::cout << "A (spectral):" << std::endl;
				a.print_spectralData_zeroNumZero();

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

				if (dst_res/2 >= freq)
				{
					std::cout << "dst_res/2 >= freq: Current frequency is representable in new array" << std::endl;


					PlaneData test(planeDataConfigDst);
					setupData(test, freq, freq);

					std::cout << "Spectral data of b:" << std::endl;
					b.print_spectralData_zeroNumZero();

					std::cout << "Spectral data of test:" << std::endl;
					test.print_spectralData_zeroNumZero();

					error = (b-test).reduce_maxAbs();
					if (error > epsilon)
					{
						std::cout << "Spectral data of A:" << std::endl;
						a.print_spectralData_zeroNumZero();

						std::cout << "Spectral data of b:" << std::endl;
						b.print_spectralData_zeroNumZero();

						std::cout << "Spectral data of test:" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << "Error: " << error << std::endl;
						FatalError("No modes changed! Results should be identical!");
					}

					std::cout << "PASSED with error of " << error << std::endl;
				}
				else
				{
					std::cout << "dst_res/2 < freq: Modes should be truncated => zero" << std::endl;

					/*
					 * Reduced resolution => All modes should be truncated and the result should be zero
					 */
					error = b.reduce_maxAbs();
					if (error > epsilon)
					{
						std::cout << "Spectral data of A:" << std::endl;
						a.print_spectralData_zeroNumZero();

						std::cout << "Spectral data of B:" << std::endl;
						b.print_spectralData_zeroNumZero();

						std::cout << "Error: " << error << std::endl;
						FatalError("Reduction and truncation of modes! Resulting array should be zero!");
					}

					std::cout << "PASSED with error of " << error << std::endl;
				}
			}
		}
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
