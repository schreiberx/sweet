/*
 *  Created on: Feb 19, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRAL_TO_PLANEDATAPHYSICAL_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRAL_TO_PLANEDATAPHYSICAL_HPP_

#include "PlaneDataConfig.hpp"
#include "PlaneData_Physical.hpp"
#include "PlaneData_Spectral.hpp"

namespace sweet{

void convert(
		const PlaneData_Spectral& i_planeDataSpectral,
		PlaneData_Physical& o_planeDataPhysical
)
{
	const PlaneDataConfig* planeDataConfig = i_planeDataSpectral.planeDataConfig;
	PlaneData_Physical o_data(planeDataConfig);

	/**
	 * Warning: The fftw functions are in-situ operations.
	 * Therefore, the data in the source array will be destroyed.
	 * Hence, we create a copy
	 */
	PlaneData_Physical tmp(planeDataConfig);
	planeDataConfig->fft_spectral_to_physical(i_planeDataSpectral.spectral_space_data, o_planeDataPhysical.physical_space_data);
}

}

#endif
