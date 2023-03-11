/*
 *  Created on: 14 Nov 2019
 *      Author: Pedro Peixoto <ppeixoto@usp.br>
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <sweet/core/defaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/ProgramArguments.hpp>

#include "pde_swePlane_normalModeConversion/SWEPlaneNormalModes.hpp"

#include <ostream>
#include <cmath>


typedef std::complex<double> complex;


int main(
		int i_argc,
		char *i_argv[]
)
{
	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackPlaneDataOps *shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	SWEPlaneNormalModes swePlaneNormalModes;
	swePlaneNormalModes.shackRegistration(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.printShackData();

	sweet::PlaneData_Config planeDataConfig;
	planeDataConfig.setupAuto(shackPlaneDataOps);

	sweet::PlaneOperators ops;
	ops.setup(planeDataConfig, shackPlaneDataOps);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(ops);



	double eps = 1e-9;

	std::cout << "*************************************************************" << std::endl;
	std::cout << "Testing normal mode conversion tool" << std::endl;
	std::cout << "*************************************************************" << std::endl;


	std::cout << "*************************************************************" << std::endl;
	planeDataConfig.printInformation();
	std::cout << "*************************************************************" << std::endl;
	std::cout << std::endl;

	sweet::PlaneData_Spectral h(planeDataConfig);
	sweet::PlaneData_Spectral u(planeDataConfig);
	sweet::PlaneData_Spectral v(planeDataConfig);

	//Normal modes to be generated
	double geo_mode=1.0;
	double igwest_mode=1.0;
	double igeast_mode=0.0;
	
	h.spectral_set_zero();
	u.spectral_set_zero();
	v.spectral_set_zero();

	std::cout << "Adding normal mode with coefficients:" << std::endl;
	std::cout << " Geost: " << geo_mode<<std::endl;
	std::cout << " IGWest: " << igwest_mode<<std::endl;
	std::cout << " IGEast: " << igeast_mode<<std::endl;
	std::cout << "To wavenumber  ";

	for (std::size_t ik1 = 0; ik1 < planeDataConfig.spectral_data_size[1]/4; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < planeDataConfig.spectral_data_size[0]/2; ik0++)
		{
			std::cout << " (" << ik0 << "," << ik1<< ") , ";

			swePlaneNormalModes.add_normal_mode(
									ik0, ik1,
									geo_mode,
									igwest_mode,
									igeast_mode,
									h,
									u,
									v
							);
		}
	}
	
	std::cout << "\n\n Extracting Normal Modes:" << std::endl;

	sweet::PlaneData_Spectral geo(planeDataConfig);
	sweet::PlaneData_Spectral igwest(planeDataConfig);
	sweet::PlaneData_Spectral igeast(planeDataConfig);
	complex geo_mode_c;
	complex igwest_mode_c;
	complex igeast_mode_c;

	swePlaneNormalModes.convert_allspectralmodes_to_normalmodes(
		h,	
		u,
		v,
		geo,
		igwest,
		igeast
	);

	for (std::size_t ik1 = 0; ik1 < planeDataConfig.spectral_data_size[1]/4; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < planeDataConfig.spectral_data_size[0]/2; ik0++)
		{
			std::cout << "From wavenumber (" << ik0 << "," << ik1<< "):";
			geo_mode_c = geo.spectral_get(ik1, ik0);
			igwest_mode_c = igwest.spectral_get(ik1, ik0);
			igeast_mode_c = igeast.spectral_get(ik1, ik0);

			std::cout << " Geost: " << geo_mode_c;
			std::cout << " IGWest: " << igwest_mode_c;
			std::cout << " IGEast: " << igeast_mode_c;

			double error = abs(geo_mode_c.real()-geo_mode)+abs(igwest_mode_c.real()-igwest_mode)+abs(igeast_mode_c.real()-igeast_mode);
			std::cout << " Error: " << error ;
			if (error<10e-10)
			{
				std::cout << " PASSED " << std::endl;
			}
			else //mode may have been split to mirror
			{
				error = abs(geo_mode_c.real()-geo_mode/2.0)+abs(igwest_mode_c.real()-igwest_mode/2.0)+abs(igeast_mode_c.real()-igeast_mode/2.0);
				if (error < eps)
				{
					std::cout << " PASSED " << std::endl;
				}
				else
				{
					std::cout << " FAIL " << std::endl;
					std::cerr << "NORMAL_MODES: Failed to convert back and forward to/from normal modes!" << error << std::endl;
					exit(-1);
				}
			}
		}
	}	

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
