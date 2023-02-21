/*
 * PlaneDataPhysical_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_

#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/ScalarDataArray.hpp>

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

////#if SWEET_USE_PLANE_SPECTRAL_SPACE
////		if (i_planeData.spectral_space_data_valid && i_raise_error_if_spectral)
////			SWEETError("This data should be typically never converted to spectral space");
////#endif

////		i_planeData.request_data_physical();

		for (std::size_t i = 0; i < out.number_of_elements; i++)
			out.scalar_data[i] = i_planeData.physical_space_data[i];

		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATAPHYSICAL_TO_SCALARDATAARRAY_HPP_ */
