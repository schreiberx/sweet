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
#if 1
//#error "Don't use physical_convert, since it looses the highest modes"
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
#endif


#if SWEET_USE_PLANE_SPECTRAL_SPACE
public:
	static
	PlaneData spectral_convert_physical_real_only(
			const PlaneDataComplex &i_planeData
	)
	{
		PlaneData out(i_planeData.planeDataConfig);
		out.spectral_set_zero();

		i_planeData.request_data_spectral();

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
					const std::complex<double> &data = i_planeData.p_spectral_get(j, i);
					out.p_spectral_set(j, i, data);
				}
			}
		}

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}
#endif
};



#endif /* SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_SCALARDATAARRAY_HPP_ */
