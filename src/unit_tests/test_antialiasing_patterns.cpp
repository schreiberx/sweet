
//#if !SWEET_USE_PLANE_SPECTRAL_SPACE
//	#error "Spectral space not activated"
//#endif

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


int main(int i_argc, char *i_argv[])
{
	// override flag
	SimulationVariables simVars;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.disc.use_spectral_basis_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.res_physical[0];
	std::size_t res_y = simVars.disc.res_physical[1];

	std::size_t max_res = 64;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	for (; res_x <= max_res && res_y <= max_res; res_x+=2, res_y+=2)
	{
		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing aliasing pattern with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		simVars.disc.res_physical[0] = res[0];
		simVars.disc.res_physical[1] = res[1];
		simVars.reset();

		planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical);

		planeDataConfigInstance.printInformation();

		std::cout << "PHYS RES: " << planeDataConfig->physical_res[0] << " x " << planeDataConfig->physical_res[1] << std::endl;
		std::cout << "PHYS SIZE: " << planeDataConfig->physical_data_size[0] << " x " << planeDataConfig->physical_data_size[1] << std::endl;
		std::cout << "SPEC MODES: " << planeDataConfig->spectral_modes[0] << " x " << planeDataConfig->spectral_modes[1] << std::endl;
		std::cout << "SPEC SIZE: " << planeDataConfig->spectral_data_size[0] << " x " << planeDataConfig->spectral_data_size[1] << std::endl;
		std::cout << std::endl;

#define PRINT_SPECTRUM	0

#if PRINT_SPECCTRUM
		int res_max = 32;
#endif


		PlaneData h(planeDataConfig);

		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////

		h.spectral_set_zero();

		h.spectral_space_data_valid = true;
		h.physical_space_data_valid = false;

		h = h.spectral_addScalarAll(1.0);

#if PRINT_SPECCTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "All one spectrum:" << std::endl;
			h.print_spectralData_zeroNumZero();
			std::cout << std::endl;
		}
#endif
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////

		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
			h.spectral_space_data[i] = {1.0,0.0};

		h.spectral_space_data_valid = true;
		h.physical_space_data_valid = false;

#if PRINT_SPECCTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "All one spectrum:" << std::endl;
			h.print_spectralData();
			std::cout << std::endl;
		}
#endif

		h.spectral_zeroAliasingModes();

#if PRINT_SPECCTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "Zero aliasing spectrum:" << std::endl;
			h.print_spectralData();
			std::cout << std::endl;
		}
#endif

		h = h.spectral_addScalarAll(1.0);

#if PRINT_SPECCTRUM
		if (res[0] < res_max && res[1] < res_max)
		{
			std::cout << "***************************************" << std::endl;
			std::cout << "Add scalar 1 to non-aliasing spectrum:" << std::endl;
			h.print_spectralData();
			std::cout << std::endl;
			std::cout << "(There should be only 2 and 0)" << std::endl;
		}
#endif

		for (std::size_t i = 0; i < planeDataConfig->spectral_array_data_number_of_elements; i++)
		{
			if (h.spectral_space_data[i] != 2.0 && h.spectral_space_data[i] != 0.0)
			{
				FatalError("INCONSISTENT ALIASING !!!");
			}
		}
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
