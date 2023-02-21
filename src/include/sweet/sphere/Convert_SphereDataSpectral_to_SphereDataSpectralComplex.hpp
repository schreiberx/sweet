/*
 * SphereData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRAL_TO_SPHEREDATASPECTRALCOMPLEX_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRAL_TO_SPHEREDATASPECTRALCOMPLEX_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_SphereDataSpectral_To_SphereDataSpectralComplex
{
public:
	static
	SphereData_SpectralComplex physical_convert(
			const SphereData_Spectral &i_sphereData
	)
	{
		SphereData_Physical tmp = i_sphereData.toPhys();
		SphereData_PhysicalComplex tmpc(i_sphereData.sphereDataConfig);

#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < tmp.sphereDataConfig->physical_array_data_number_of_elements; i++)
			tmpc.physical_space_data[i] = tmp.physical_space_data[i];


		SphereData_SpectralComplex ret(tmpc);
		return ret;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
