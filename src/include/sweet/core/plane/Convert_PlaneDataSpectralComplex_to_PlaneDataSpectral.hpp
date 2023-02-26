/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRALCOMPLEX_TO_PLANEDATASPECTRAL_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRALCOMPLEX_TO_PLANEDATASPECTRAL_HPP_

#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/core/ScalarDataArray.hpp>

namespace sweet
{

class Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral
{
public:
	static
	PlaneData_Spectral physical_convert_real(
			const PlaneData_SpectralComplex &i_planeData
	)
	{
		PlaneData_PhysicalComplex tmp_cplx = i_planeData.toPhys();
		PlaneData_Physical tmp(i_planeData.planeDataConfig);

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
		PlaneData_Physical tmp(i_planeData.planeDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.planeDataConfig->physical_array_data_number_of_elements; i++)
			tmp.physical_space_data[i] = tmp_cplx.physical_space_data[i].imag();

		return PlaneData_Spectral(tmp);
	}

public:
	static
	PlaneData_Spectral spectral_convert_physical_real_only(
			const PlaneData_SpectralComplex &i_planeData
	)
	{
		PlaneData_Spectral out(i_planeData.planeDataConfig);
		out.spectral_set_zero();

		for (int r = 0; r < 2; r++)
		{
			for (	std::size_t j = out.planeDataConfig->spectral_data_iteration_ranges[r][1][0];
					j < out.planeDataConfig->spectral_data_iteration_ranges[r][1][1];
					j++
			) {
				for (	std::size_t i = out.planeDataConfig->spectral_data_iteration_ranges[r][0][0];
						i < out.planeDataConfig->spectral_data_iteration_ranges[r][0][1];
						i++
				) {
					const std::complex<double> &data = i_planeData.spectral_get(j, i);
					out.spectral_set(j, i, data);
				}
			}
		}

		return out;
	}

};

}

#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRALCOMPLEX_TO_PLANEDATASPECTRAL_HPP_ */
