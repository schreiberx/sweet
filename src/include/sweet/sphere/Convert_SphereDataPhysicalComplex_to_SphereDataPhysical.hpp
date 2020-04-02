/*
 * SphereData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICALCOMPLEX_TO_SPHEREDATAPHYSICAL_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICALCOMPLEX_TO_SPHEREDATAPHYSICAL_HPP_

#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_SphereDataPhysicalComplex_To_SphereDataPhysical
{
public:
	static
	SphereData_Physical physical_convert_real(
			const SphereData_PhysicalComplex &i_sphereData
	)
	{
		SphereData_Physical out(i_sphereData.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_sphereData.physical_space_data[i].real();

		return out;
	}



public:
	static
	SphereData_Physical physical_convert_imag(
			const SphereData_PhysicalComplex &i_sphereData
	)
	{
		SphereData_Physical out(i_sphereData.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_sphereData.physical_space_data[i].imag();


		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
