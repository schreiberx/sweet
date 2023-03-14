/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_

#include <sweet/core/sphere/SphereData_Physical.hpp>
#include <sweet/core/ScalarDataArray.hpp>


namespace sweet
{

class Convert_SphereDataPhysical_to_ScalarDataArray
{
public:
	static
	ScalarDataArray physical_convert(
			const SphereData_Physical &i_sphereData
	)
	{
		ScalarDataArray out(i_sphereData.sphereDataConfig->physical_array_data_number_of_elements);

		for (std::size_t i = 0; i < out.number_of_elements; i++)
			out.scalar_data[i] = i_sphereData.physical_space_data[i];

		return out;
	}
};

}

#endif
