/*
 * SphereDataSpectral_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRALCOMPLEX_TO_SPHEREDATASPECTRAL_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRALCOMPLEX_TO_SPHEREDATASPECTRAL_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/core/ScalarDataArray.hpp>

class Convert_SphereDataSpectralComplex_To_SphereDataSpectral
{
public:
	static
	SphereData_Spectral physical_convert_real(
			const SphereData_SpectralComplex &i_sphereData
	)
	{
		SphereData_PhysicalComplex tmp_cplx = i_sphereData.toPhys();
		SphereData_Physical tmp(i_sphereData.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.sphereDataConfig->physical_array_data_number_of_elements; i++)
			tmp.physical_space_data[i] = tmp_cplx.physical_space_data[i].real();

		return SphereData_Spectral(tmp);
	}



public:
	static
	SphereData_Spectral physical_convert_imag(
			const SphereData_SpectralComplex &i_sphereData
	)
	{
		SphereData_PhysicalComplex tmp_cplx = i_sphereData.toPhys();
		SphereData_Physical tmp(i_sphereData.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.sphereDataConfig->physical_array_data_number_of_elements; i++)
			tmp.physical_space_data[i] = tmp_cplx.physical_space_data[i].imag();

		return SphereData_Spectral(tmp);
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
