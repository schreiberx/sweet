/*
 * SphereData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICALCOMPLEX_TO_SPHEREDATAPHYSICAL_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICALCOMPLEX_TO_SPHEREDATAPHYSICAL_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_SphereDataPhysicalComplex_To_SphereDataPhysical
{
public:
	static
	SphereDataPhysical physical_convert_real(
			const SphereDataPhysicalComplex &i_sphereData
	)
	{
		SphereDataPhysical out(i_sphereData.sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int i = 0; i < out.sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_sphereData.physical_space_data[i].real();

		return out;
	}



public:
	static
	SphereDataPhysical physical_convert_imag(
			const SphereDataPhysicalComplex &i_sphereData
	)
	{
		SphereDataPhysical out(i_sphereData.sphereDataConfig);

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (int i = 0; i < out.sphereDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_sphereData.physical_space_data[i].imag();


		return out;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
