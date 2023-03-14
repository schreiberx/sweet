/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICALCOMPLEX_TO_SPHEREDATAPHYSICAL_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_CONVERT_SPHEREDATAPHYSICALCOMPLEX_TO_SPHEREDATAPHYSICAL_HPP_

#include <sweet/core/sphere/SphereData_Physical.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/core/ScalarDataArray.hpp>


namespace sweet
{

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

}

#endif
