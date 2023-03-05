/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_OUTPUTSPHERICALHARMONICS_HPP_
#define SRC_OUTPUTSPHERICALHARMONICS_HPP_


#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <cassert>
#include <sweet/core/shacks/ShackDictionary.hpp>




class AppOutputSphericalHarmonics
{
public:
	void run(sweet::SphereDataConfig *sphereConfig)
	{
		sweet::SphereData_Spectral h(sphereConfig);

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
	sweet::ShackDictionary shackDict;

	if (!shackDict.setupFromMainParameters(i_argc, i_argv))
		return -1;

	sweet::SphereDataConfig sphereDataConfig;

	if (shackDict.disc.space_res_physical[0] <= 0)
	{
		sphereDataConfig.setupAutoPhysicalSpace(
				shackDict.disc.space_res_spectral[0], shackDict.disc.space_res_spectral[1],		// spectral (lon/lat)
				&shackDict.disc.space_res_physical[0], &shackDict.disc.space_res_physical[1],	// physical (lon/lat)
				shackDict.misc.reuse_spectral_transformation_plans,
				shackDict.misc.verbosity,
				shackDict.parallelization.num_threads_space
		);
	}
	else
	{
		sphereDataConfig.setup(
				shackDict.disc.space_res_spectral[0], shackDict.disc.space_res_spectral[1],		// spectral (lon/lat)
				shackDict.disc.space_res_physical[0], shackDict.disc.space_res_physical[1],	// physical (lon/lat)
				shackDict.misc.reuse_spectral_transformation_plans,
				shackDict.misc.verbosity,
				shackDict.parallelization.num_threads_space
			);
	}


	AppOutputSphericalHarmonics sph;
	sph.run(&sphereDataConfig);

	return 1;
}


#endif
