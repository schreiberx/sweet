/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATACOMPLEX_TO_PLANEDATA_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATACOMPLEX_TO_PLANEDATA_HPP_

#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical
{
public:
	static
	PlaneData_Physical physical_convert_real(
			const PlaneData_PhysicalComplex &i_planeData
	)
	{
		PlaneData_Physical out(i_planeData.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_planeData.physical_space_data[i].real();

		return out;
	}



public:
	static
	PlaneData_Physical physical_convert_imag(
			const PlaneData_PhysicalComplex &i_planeData
	)
	{
		PlaneData_Physical out(i_planeData.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_planeData.physical_space_data[i].imag();


		return out;
	}

};

#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_ */
