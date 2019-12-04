/*
 * SWE_bench_NormalModes.hpp
 *
 *  Created on: 03 Nov 2019
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef SWE_PLANE_NORMAL_MODES_HPP_
#define SWE_PLANE_NORMAL_MODES_HPP_

#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>

#include <sweet/sweetmath.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <libmath/GaussQuadrature.hpp>

#include "swe_plane/SWE_Plane_Normal_Modes.hpp"

/**
 * Implement normal mode initialization
 *
 *
 *
 **/
class SWE_bench_NormalModes
{

	SimulationVariables &simVars;

	PlaneOperators &op;

	double f = simVars.sim.plane_rotating_f0;
	double g = simVars.sim.gravitation;
	double sx = simVars.sim.plane_domain_size[0];
	double sy = simVars.sim.plane_domain_size[1];

	std::string bcasename;
	std::size_t k0, k1;
	double d0, dwest, deast;

public:
	SWE_bench_NormalModes(
		SimulationVariables &io_simVars,
		PlaneOperators &io_op
	)	:
		simVars(io_simVars),
		op(io_op)
	{
	}

	void setup(
			PlaneData &o_h,
			PlaneData &o_u,
			PlaneData &o_v
	)
	{
		std::cout<< "Generating Normal Modes Initial Conditions: ";

		extract_bench_info(simVars.benchmark.benchmark_normal_modes_case);
		
		//Naming convention for normal_mode_cases (for arbitraty use):
		//  name_k0_k1_d0_deast_dwest
		// (k0,k1) are wave numbers (set to -1 for all wavenumbers)
		// d0, dwest, deast are numbers (floats) that are coefficients for different normal wave types

		//zero initial conditions
		const PlaneDataConfig *planeDataConfig = o_h.planeDataConfig;

		o_h.spectral_set_zero();
		o_u.spectral_set_zero();
		o_v.spectral_set_zero();

		if(bcasename=="single")
		{
			//Set a single wavenumber with appropriate modes
			if(k0<0 && k1 <0 ){
				std::cout<<"Adding normal modes to all wavenumbers"<<std::endl;
				for (std::size_t ik1 = 0; ik1 < planeDataConfig->spectral_data_size[1]; ik1++)
				{
					for (std::size_t ik0 = 0; ik0 < planeDataConfig->spectral_data_size[0]; ik0++)
					{

						SWE_Plane_Normal_Modes::add_normal_mode(
											ik0, ik1,
											d0,
											dwest,
											deast,
											o_h,
											o_u,
											o_v,
											simVars
									);
					}
				}
			}
			else{
				if(k0>0 && k1 >0 ){
					SWE_Plane_Normal_Modes::add_normal_mode(
											k0, k1,
											d0,
											dwest,
											deast,
											o_h,
											o_u,
											o_v,
											simVars
									);
				}
				else
				{
					FatalError("SWE_bench_NormalModes: invalid wavenumber selection in --benchmark-normal-mode-case [string] (see SWE_bench_NormalModes.hpp file)");		
				}
			}
		};
		
		std::cout<< "   Done! " << std::endl;
	}

	//not in use
	private:
	void extract_bench_info(const std::string &bcase)
	{
		
		if(bcase==""){
			FatalError("SWE_bench_NormalModes: please choose the normal mode case with --benchmark-normal-mode-case [string] (see SWE_bench_NormalModes.hpp file)");
		};
		std::cout<< bcase <<std::endl;

		//Convert parameter to words
		std::string str = bcase;
		std::replace( str.begin(), str.end(), '_', ' ');
		//std::cout<< str<<std::endl;

		std::stringstream iss(str);
		iss >> bcasename;
		std::cout<< "Benchmark case: "<< bcasename << std::endl;
		iss >> k0;
		iss >> k1;
		std::cout<< "Wavenumbers: "<< k0 << "," <<k1 << std::endl;
		iss >> d0;
		iss >> dwest;
		iss >> deast;
		std::cout<< "Normal mode coefficients:" << d0 << " " << dwest << " " << deast << std::endl;
		FatalError("SWE_bench_NormalModes: please choose the normal mode case with --benchmark-normal-mode-case [string] (see SWE_bench_NormalModes.hpp file)");
		return;
		
	}

};


#endif
