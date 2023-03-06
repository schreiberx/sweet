/*
 * PlaneDataPhysical_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_PLANE_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_

#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneDataConfig.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_ScalarDataArray_to_PlaneDataPhysical
{
public:
	static
	PlaneData_Physical convert(
			const ScalarDataArray &i_scalarDataArray,
			const PlaneDataConfig *i_planeDataConfig
	)
	{
		PlaneData_Physical out(i_planeDataConfig);

		for (std::size_t i = 0; i < out.planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_scalarDataArray.scalar_data[i];

////////#if SWEET_USE_PLANE_SPECTRAL_SPACE
////////		out.physical_space_data_valid = true;
////////		out.spectral_space_data_valid = false;
////////#endif
		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_ */
