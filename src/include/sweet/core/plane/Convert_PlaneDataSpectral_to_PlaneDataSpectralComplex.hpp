/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRAL_TO_PLANEDATASPECTRALCOMPLEX_HPP_
#define SRC_INCLUDE_SWEET_PLANE_CONVERT_PLANEDATASPECTRAL_TO_PLANEDATASPECTRALCOMPLEX_HPP_

#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/core/ScalarDataArray.hpp>

namespace sweet
{

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


public:
	static
	PlaneData_SpectralComplex spectral_convert(
			const PlaneData_Spectral &i_planeData
	)
	{
		PlaneData_SpectralComplex out(i_planeData.planeDataConfig);
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

				for (	std::size_t i = out.planeDataConfig->spectral_data_iteration_ranges[r][0][0]+1;
						i < out.planeDataConfig->spectral_data_iteration_ranges[r][0][1];
						i++
				) {
					const std::complex<double> &data = i_planeData.spectral_get(j, i);
					std::complex<double> data2 = data;

					data2.imag(-data2.imag());
					if (j == 0)
						out.spectral_set(j, out.planeDataConfig->spectral_complex_data_size[0]-i, data2);
					else
						out.spectral_set(out.planeDataConfig->spectral_complex_data_size[1]-j, out.planeDataConfig->spectral_complex_data_size[0]-i, data2);
				}
			}
		}

		return out;
	}


};

}

#endif
