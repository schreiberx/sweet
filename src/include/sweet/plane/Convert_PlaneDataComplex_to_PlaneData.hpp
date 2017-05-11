/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATACOMPLEX_TO_PLANEDATA_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATACOMPLEX_TO_PLANEDATA_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_PlaneDataComplex_To_PlaneData
{
public:
	static
	PlaneData physical_convert(
			const PlaneDataComplex &i_planeData
	)
	{
		PlaneData out(i_planeData.planeDataConfig);

		i_planeData.request_data_physical();

		for (std::size_t i = 0; i < out.planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_planeData.physical_space_data[i].real();

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;
#endif

		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_ */
