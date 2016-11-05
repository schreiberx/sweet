/*
 * OutputSphericalHarmonics.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_OUTPUTSPHERICALHARMONICS_HPP_
#define SRC_OUTPUTSPHERICALHARMONICS_HPP_


#include <cassert>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/SimulationVariables.hpp>



SimulationVariables simVars;


class AppOutputSphericalHarmonics
{
public:
	void run(SphereDataConfig *sphereConfig)
	{
		SphereData h(sphereConfig);

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

				h.physical_file_write(buffer);

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

	SphereDataConfig sphereDataConfig;

	if (simVars.disc.res_physical[0] <= 0)
	{
		sphereDataConfig.setupAutoPhysicalSpace(
				simVars.disc.res_spectral[0], simVars.disc.res_spectral[1],		// spectral (lon/lat)
				&simVars.disc.res_physical[0], &simVars.disc.res_physical[1]	// physical (lon/lat)
		);
	}
	else
	{
		sphereDataConfig.setup(
				simVars.disc.res_spectral[0], simVars.disc.res_spectral[1],		// spectral (lon/lat)
				simVars.disc.res_physical[0], simVars.disc.res_physical[1]	// physical (lon/lat)
			);
	}


	AppOutputSphericalHarmonics sph;
	sph.run(&sphereDataConfig);

	return 1;
}


#endif /* SRC_OUTPUTSPHERICALHARMONICS_HPP_ */
