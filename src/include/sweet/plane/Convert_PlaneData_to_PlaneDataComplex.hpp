/*
 * PlaneData_To_ScalarDataArray.cpp
 *
 *  Created on: 20 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_PLANEDATACOMPLEX_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATA_TO_PLANEDATACOMPLEX_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/ScalarDataArray.hpp>

class Convert_PlaneData_To_PlaneDataComplex
{

#if 1
//#error "Don't use physical_convert, since it looses the highest modes"
public:
	static
	PlaneDataComplex physical_convert(
			const PlaneData &i_planeData
	)
	{
		PlaneDataComplex out(i_planeData.planeDataConfig);

		i_planeData.request_data_physical();

		for (std::size_t i = 0; i < out.planeDataConfig->physical_array_data_number_of_elements; i++)
			out.physical_space_data[i] = i_planeData.physical_space_data[i];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		out.physical_space_data_valid = true;
		out.spectral_space_data_valid = false;
#endif

#if SWEET_DEBUG
		out.test_realphysical();
#endif

		return out;
	}
#endif


#if SWEET_USE_PLANE_SPECTRAL_SPACE
public:
	static
	PlaneDataComplex spectral_convert(
			const PlaneData &i_planeData
	)
	{
		PlaneDataComplex out(i_planeData.planeDataConfig);
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

				for (	std::size_t i = out.planeDataConfig->spectral_data_iteration_ranges[r][0][0]+1;
						i < out.planeDataConfig->spectral_data_iteration_ranges[r][0][1];
						i++
				) {
					const std::complex<double> &data = i_planeData.p_spectral_get(j, i);
					std::complex<double> data2 = data;

					data2.imag(-data2.imag());
					if (j == 0)
						out.p_spectral_set(j, out.planeDataConfig->spectral_complex_data_size[0]-i, data2);
					else
						out.p_spectral_set(out.planeDataConfig->spectral_complex_data_size[1]-j, out.planeDataConfig->spectral_complex_data_size[0]-i, data2);
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
