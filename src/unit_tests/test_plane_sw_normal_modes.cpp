/*
 * test_plane_sw_normal_mode.cpp
 *
 *  Created on: 14 Nov 2019
 *      Author: Pedro Peixoto <ppeixoto@usp.br>
 */



#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "test_plane_sw_normal_modes/SWE_Plane_Normal_Modes.hpp"

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

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

	std::cout << "*************************************************************" << std::endl;
	std::cout << "Testing normal mode conversion tool" << std::endl;
	std::cout << "*************************************************************" << std::endl;
	//std::size_t res[2] = {res_x, res_y};

	//std::size_t test_max_freqx = planeDataConfig->spectral_real_modes[0];
	//std::size_t test_max_freqy = planeDataConfig->spectral_real_modes[1];

	std::cout << "*************************************************************" << std::endl;
	planeDataConfig->printInformation();
	std::cout << "*************************************************************" << std::endl;
	std::cout << std::endl;

	PlaneData h(planeDataConfig);
	PlaneData u(planeDataConfig);
	PlaneData v(planeDataConfig);

	// Spectral mode to be converted (x,y) respectively
	std::size_t k0, k1;

	double igwest_mode=1.0;
	double igeast_mode=0.0;
	double geo_mode=0.0;

	//Normal modes to be generated
	
	//std::size_t max_modes_x = planeDataConfig->spectral_real_modes[0];
	//std::size_t max_modes_y = planeDataConfig->spectral_real_modes[1];
	k0=5;
	k1=3;
	h.spectral_set_zero();
	u.spectral_set_zero();
	v.spectral_set_zero();
	//std::cout<< "out:" << h.planeDataConfig->spectral_data_size[1] << std::endl;

	SWE_Plane_Normal_Modes::add_normal_mode(
								k0, k1,
								igwest_mode,
								igeast_mode,
								geo_mode,
								h,
								u,
								v,
								simVars
						);
		

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
