/*
 * PlaneDataPhysical_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_PLANE_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_

#include <sweet/core/plane/PlaneData_Config.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/ScalarDataArray.hpp>

namespace sweet
{

class Convert_ScalarDataArray_to_PlaneDataPhysical
{
public:
	static
	PlaneData_Physical convert(
			const ScalarDataArray &i_scalarDataArray,
			const PlaneData_Config *i_planeDataConfig
	)
	{
		PlaneData_Physical out(i_planeDataConfig);

		for (std::size_t i = 0; i < out.planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_scalarDataArray.scalar_data[i];

		return out;
	}
};

}

#endif
