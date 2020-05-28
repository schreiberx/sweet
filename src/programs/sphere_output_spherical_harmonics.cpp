/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_OUTPUTSPHERICALHARMONICS_HPP_
#define SRC_OUTPUTSPHERICALHARMONICS_HPP_


#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <cassert>
#include <sweet/SimulationVariables.hpp>



SimulationVariables simVars;


class AppOutputSphericalHarmonics
{
public:
	void run(SphereData_Config *sphereConfig)
	{
		SphereData_Spectral h(sphereConfig);

		int counter = 0;
		// iterate over modes
		for (int n = 0; n <= sphereConfig->spectral_modes_n_max; n++)
		{
			for (int m = 0; m <= std::min((int)sphereConfig->spectral_modes_m_max, n); m++)
			{
				h.spectral_update_lambda(
						[&](int i_n, int i_m, std::complex<double> &o_data)
						{
							if (i_n == n && i_m == m)
								o_data = 1;
							else
								o_data = 0;
						}
					);

				char buffer[1024];
				sprintf(buffer, "SPH_n%i_m%i.csv", n, m);

				h.toPhys().physical_file_write(buffer);

				counter++;
			}
		}

		assert(counter == sphereConfig->spectral_array_data_number_of_elements);
	}


};




int main(int i_argc, char *i_argv[])
{
	MemBlockAlloc::setup();

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			nullptr
	};

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
		return -1;

	SphereData_Config sphereDataConfig;

	if (simVars.disc.space_res_physical[0] <= 0)
	{
		sphereDataConfig.setupAutoPhysicalSpace(
				simVars.disc.space_res_spectral[0], simVars.disc.space_res_spectral[1],		// spectral (lon/lat)
				&simVars.disc.space_res_physical[0], &simVars.disc.space_res_physical[1],	// physical (lon/lat)
				simVars.misc.reuse_spectral_transformation_plans
		);
	}
	else
	{
		sphereDataConfig.setup(
				simVars.disc.space_res_spectral[0], simVars.disc.space_res_spectral[1],		// spectral (lon/lat)
				simVars.disc.space_res_physical[0], simVars.disc.space_res_physical[1],	// physical (lon/lat)
				simVars.misc.reuse_spectral_transformation_plans
			);
	}


	AppOutputSphericalHarmonics sph;
	sph.run(&sphereDataConfig);

	return 1;
}


#endif /* SRC_OUTPUTSPHERICALHARMONICS_HPP_ */
