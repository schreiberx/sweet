/*
 * SphereDataSpectral_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_SPHEREDATA_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_SPHEREDATA_TO_SCALARDATAARRAY_HPP_

#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/ScalarDataArray.hpp>

class Convert_ScalarDataArray_to_SphereDataPhysical
{
public:
	static
	SphereData_Physical convert(
			const ScalarDataArray &i_scalarDataArray,
			const SphereData_Config *i_sphereDataConfig
	)
	{
		SphereData_Physical out(i_sphereDataConfig);

		assert(out.sphereDataConfig->physical_array_data_number_of_elements == i_scalarDataArray.number_of_elements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < (std::size_t)out.sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_scalarDataArray.scalar_data[i];

		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
