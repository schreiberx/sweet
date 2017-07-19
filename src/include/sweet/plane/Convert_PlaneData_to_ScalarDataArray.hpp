/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_PlaneData_To_ScalarDataArray
{
public:
	static
	ScalarDataArray physical_convert(
			const PlaneData &i_planeData
	)
	{
		ScalarDataArray out(i_planeData.planeDataConfig->physical_array_data_number_of_elements);

		if (i_planeData.spectral_space_data_valid)
			FatalError("This data should be typically never converted to spectral space");

		i_planeData.request_data_physical();

		for (std::size_t i = 0; i < out.number_of_elements; i++)
			out.scalar_data[i] = i_planeData.physical_space_data[i];

		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_ */
