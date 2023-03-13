/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_

#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/ScalarDataArray.hpp>

namespace sweet
{

class Convert_PlaneDataPhysical_To_ScalarDataArray
{
public:
	static
	ScalarDataArray physical_convert(
			const PlaneData_Physical &i_planeData,
			bool i_raise_error_if_spectral = true
	)
	{
		ScalarDataArray out(i_planeData.planeDataConfig->physical_array_data_number_of_elements);

		for (std::size_t i = 0; i < out.number_of_elements; i++)
			out.scalar_data[i] = i_planeData.physical_space_data[i];

		return out;
	}
};

}

#endif
