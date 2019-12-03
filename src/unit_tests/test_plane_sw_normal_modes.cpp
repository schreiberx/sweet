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


#if SWEET_QUADMATH && 0
	typedef __float128 T;
#else
	typedef double T;
#endif
	typedef std::complex<T> complex;


int main(
		int i_argc,
		char *i_argv[]
)
{
	// override flag
	SimulationVariables simVars;

	double eps = 1e-9;

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

	//Normal modes to be generated
	double geo_mode=1.0;
	double igwest_mode=1.0;
	double igeast_mode=0.0;
	
	h.spectral_set_zero();
	u.spectral_set_zero();
	v.spectral_set_zero();
	std::cout<< "Adding normal mode with coefficients:"<< std::endl;
	std::cout<< " Geost: "<< geo_mode<<std::endl;
	std::cout<< " IGWest: "<< igwest_mode<<std::endl;
	std::cout<< " IGEast: "<< igeast_mode<<std::endl;
	std::cout<< "To wavenumber  ";
	for (std::size_t ik1 = 0; ik1 < planeDataConfig->spectral_data_size[1]/4; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < planeDataConfig->spectral_data_size[0]/2; ik0++)
		{
			std::cout<< " (" << ik0 << ","<< ik1<<") , ";

			SWE_Plane_Normal_Modes::add_normal_mode(
									ik0, ik1,
									geo_mode,
									igwest_mode,
									igeast_mode,
									h,
									u,
									v,
									simVars
							);
		}
	}
	//std::cout<<"spectral "<< std::endl;
	//std::cout<< std::endl;
	//v.print_spectralIndex();
	//std::cout<<"physical"<<std::endl;
	//v.print_physicalArrayData(3);
	

	//std::cout<< "Projecting fields on wavenumber (" << ik0 << ","<< ik1<<") to normal modes"<< std::endl;
	std::cout<< "\n\n Extracting Normal Modes:"<< std::endl;

	PlaneData geo(planeDataConfig);
	PlaneData igwest(planeDataConfig);
	PlaneData igeast(planeDataConfig);
	complex geo_mode_c;
	complex igwest_mode_c;
	complex igeast_mode_c;

	SWE_Plane_Normal_Modes::convert_allspectralmodes_to_normalmodes(
		h,	
		u,
		v,
		simVars,
		geo,
		igwest,
		igeast
	);

	for (std::size_t ik1 = 0; ik1 < planeDataConfig->spectral_data_size[1]/4; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < planeDataConfig->spectral_data_size[0]/2; ik0++)
		{
			std::cout<< "From wavenumber (" << ik0 << ","<< ik1<<"):";
			geo_mode_c = geo.p_spectral_get(ik1, ik0);
			igwest_mode_c = igwest.p_spectral_get(ik1, ik0);
			igeast_mode_c = igeast.p_spectral_get(ik1, ik0);

			std::cout<< " Geost: "<< geo_mode_c;
			std::cout<< " IGWest: "<< igwest_mode_c;
			std::cout<< " IGEast: "<< igeast_mode_c;
			double error = abs(geo_mode_c.real()-geo_mode)+abs(igwest_mode_c.real()-igwest_mode)+abs(igeast_mode_c.real()-igeast_mode);
			std::cout<< " Error: "<< error ;
			if(error<10e-10)
				std::cout<< " PASSED "<< std::endl;
			else //mode may have been split to mirror
			{
				error = abs(geo_mode_c.real()-geo_mode/2.0)+abs(igwest_mode_c.real()-igwest_mode/2.0)+abs(igeast_mode_c.real()-igeast_mode/2.0);
				if ( error < eps )
					std::cout<< " PASSED "<< std::endl;
				else
				{
					std::cout<< " FAIL "<< std::endl;
					std::cerr << "NORMAL_MODES: Failed to convert back and forward to/from normal modes!" << error << std::endl;
					exit(-1);
				}
			}
		}
	}	

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
