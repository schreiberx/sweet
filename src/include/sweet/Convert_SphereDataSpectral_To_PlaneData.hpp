/*
 * SphereDataSpectral_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_PLANEDATA_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_PLANEDATA_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>

class Convert_SphereDataSpectral_To_PlaneData
{
public:
	static
	PlaneData physical_convert(
			const SphereData_Spectral &i_sphereDataSpectral,
			PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(i_sphereDataSpectral.sphereDataConfig->physical_num_lon == (int)i_planeDataConfig->physical_res[0]);
		assert(i_sphereDataSpectral.sphereDataConfig->physical_num_lat == (int)i_planeDataConfig->physical_res[1]);
		assert(i_planeDataConfig->physical_array_data_number_of_elements == i_sphereDataSpectral.sphereDataConfig->physical_array_data_number_of_elements);

		SphereData_Physical i_sphereData = i_sphereDataSpectral.toPhys();

		PlaneData out(i_planeDataConfig);



#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int i = 0; i < i_sphereData.sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < i_sphereData.sphereDataConfig->physical_num_lat; j++)
				out.physical_space_data[(i_sphereData.sphereDataConfig->physical_num_lat-1-j)*i_sphereData.sphereDataConfig->physical_num_lon + i] = i_sphereData.physical_space_data[i*i_sphereData.sphereDataConfig->physical_num_lat + j];
#else
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int j = 0; j < i_sphereData.sphereDataConfig->physical_num_lat; j++)
			for (int i = 0; i < i_sphereData.sphereDataConfig->physical_num_lon; i++)
				out.physical_space_data[(i_sphereData.sphereDataConfig->physical_num_lat-1-j)*i_sphereData.sphereDataConfig->physical_num_lon + i] = i_sphereData.physical_space_data[j*i_sphereData.sphereDataConfig->physical_num_lon + i];
#endif

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;
#endif

		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
