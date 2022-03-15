/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRAL_TO_PLANEDATASPECTRALCOMPLEX_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRAL_TO_PLANEDATASPECTRALCOMPLEX_HPP_

#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex
{
public:
	static
	PlaneData_SpectralComplex physical_convert(
			const PlaneData_Spectral &i_planeData
	)
	{
		PlaneData_Physical tmp = i_planeData.toPhys();
		PlaneData_PhysicalComplex tmpc(i_planeData.planeDataConfig);

#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < tmp.planeDataConfig->physical_array_data_number_of_elements; i++)
			tmpc.physical_space_data[i] = tmp.physical_space_data[i];


		PlaneData_SpectralComplex ret(tmpc);
		return ret;
	}
};


#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_ */
