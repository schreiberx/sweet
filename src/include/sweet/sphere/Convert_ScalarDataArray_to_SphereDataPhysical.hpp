/*
 * SphereDataSpectral_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_SPHEREDATA_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_SPHEREDATA_TO_SCALARDATAARRAY_HPP_

#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_ScalarDataArray_to_SphereDataPhysical
{
public:
	static
	SphereDataPhysical convert(
			const ScalarDataArray &i_scalarDataArray,
			const SphereDataConfig *i_sphereDataConfig
	)
	{
		SphereDataPhysical out(i_sphereDataConfig);

		for (std::size_t i = 0; i < (std::size_t)out.sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_scalarDataArray.scalar_data[i];

		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
