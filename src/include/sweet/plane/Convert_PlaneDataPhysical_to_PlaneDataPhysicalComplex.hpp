/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATAPHYSICAL_TO_PLANEDATAPHYSICALCOMPLEX_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATAPHYSICAL_TO_PLANEDATAPHYSICALCOMPLEX_HPP_

#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex
{
public:
	static
	PlaneData_PhysicalComplex physical_convert(
			const PlaneData_Physical &i_planeData
	)
	{
		PlaneData_PhysicalComplex out(i_planeData.planeDataConfig);

#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < tmp.planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_planeData.physical_space_data[i];

		return out;
	}
};


#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_ */
