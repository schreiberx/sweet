/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRALCOMPLEX_TO_PLANEDATASPECTRAL_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRALCOMPLEX_TO_PLANEDATASPECTRAL_HPP_

#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral
{
public:
	static
	PlaneData_Spectral physical_convert_real(
			const PlaneData_SpectralComplex &i_planeData
	)
	{
		PlaneData_PhysicalComplex tmp_cplx = i_planeData.toPhys();
		PlaneData_Physical tmp(i_sphereData.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.planeDataConfig->physical_array_data_number_of_elements; i++)
			tmp.physical_space_data[i] = tmp_cplx.physical_space_data[i].real();

		return PlaneData_Spectral(tmp);
	}

public:
	static
	PlaneData_Spectral physical_convert_imag(
			const PlaneData_SpectralComplex &i_planeData
	)
	{
		PlaneData_PhysicalComplex tmp_cplx = i_planeData.toPhys();
		PlaneData_Physical tmp(i_sphereData.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.planeDataConfig->physical_array_data_number_of_elements; i++)
			tmp.physical_space_data[i] = tmp_cplx.physical_space_data[i].imag();

		return PlaneData_Spectral(tmp);
	}

};



#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRALCOMPLEX_TO_PLANEDATASPECTRAL_HPP_ */
