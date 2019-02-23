/*
 * SphereDataSpectral_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRALCOMPLEX_TO_SPHEREDATASPECTRAL_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATASPECTRALCOMPLEX_TO_SPHEREDATASPECTRAL_HPP_

#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereDataSpectralComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_SphereDataSpectralComplex_To_SphereDataSpectral
{
public:
	static
	SphereDataSpectral physical_convert_real(
			const SphereDataSpectralComplex &i_sphereData
	)
	{
		SphereDataPhysicalComplex tmp_cplx = i_sphereData.getSphereDataPhysicalComplex();
		SphereDataPhysical tmp(i_sphereData.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < tmp_cplx.sphereDataConfig->physical_array_data_number_of_elements; i++)
			tmp.physical_space_data[i] = tmp_cplx.physical_space_data[i].real();

		return SphereDataSpectral(tmp);
	}



public:
	static
	SphereDataSpectral physical_convert_imag(
			const SphereDataSpectralComplex &i_sphereData
	)
	{
		SphereDataPhysicalComplex tmp_cplx = i_sphereData.getSphereDataPhysicalComplex();
		SphereDataPhysical tmp(i_sphereData.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (int i = 0; i < tmp_cplx.sphereDataConfig->physical_array_data_number_of_elements; i++)
			tmp.physical_space_data[i] = tmp_cplx.physical_space_data[i].imag();

		return SphereDataSpectral(tmp);
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATA_TO_SCALARDATAARRAY_HPP_ */
