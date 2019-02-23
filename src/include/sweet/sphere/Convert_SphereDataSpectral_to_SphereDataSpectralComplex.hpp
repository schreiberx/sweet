/*
 * SphereData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRAL_TO_SPHEREDATASPECTRALCOMPLEX_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRAL_TO_SPHEREDATASPECTRALCOMPLEX_HPP_

#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereDataSpectralComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_SphereDataSpectral_To_SphereDataSpectralComplex
{
public:
	static
	SphereDataSpectralComplex physical_convert(
			const SphereDataSpectral &i_sphereData
	)
	{
		SphereDataPhysical tmp = i_sphereData.getSphereDataPhysical();
		SphereDataPhysicalComplex tmpc(i_sphereData.sphereDataConfig);

#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (int i = 0; i < tmp.sphereDataConfig->physical_array_data_number_of_elements; i++)
			tmpc.physical_space_data[i] = tmp.physical_space_data[i];


		SphereDataSpectralComplex ret(tmpc);
		return ret;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
