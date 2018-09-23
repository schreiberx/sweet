/*
 * SphereComplexErrorCheck.hpp
 *
 *  Created on: 10 Nov 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_SPHERECOMPLEXERRORCHECK_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_SPHERECOMPLEXERRORCHECK_HPP_

#include <sweet/sphere/SphereDataComplex.hpp>

class SphereComplexErrorCheck
{
	SphereComplexErrorCheck(
		SphereDataComplex &i_lhs,
		SphereDataComplex &i_rhs,
		const std::string &i_id,
		double i_error_threshold = 1.0,
		double i_ignore_error = false
	)
	{
		SphereDataComplex lhsr = i_lhs.spectral_returnWithDifferentModes(sphereDataConfig);
		SphereDataComplex rhsr = i_rhs.spectral_returnWithDifferentModes(sphereDataConfig);

		double normalize_fac = std::min(lhsr.physical_reduce_max_abs(), rhsr.physical_reduce_max_abs());

		SphereDataComplex diff = lhsr-rhsr;
		diff.physical_reduce_max_abs();

		if (normalize_fac == 0)
		{
			std::cout << "Error computation for '" << i_id << "' ignored since both fields are Zero" << std::endl;
			return;
		}
		double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.physical_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tThreshold: " << i_error_threshold << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(lhsr).physical_file_write("o_error_lhs.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(rhsr).physical_file_write("o_error_rhs.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(lhsr-rhsr).physical_file_write("o_error_diff.csv");

			if (i_ignore_error)
				std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
			else
				FatalError("Error too large");
		}
	}
};


#endif /* SRC_INCLUDE_SWEET_SPHERE_SPHERECOMPLEXERRORCHECK_HPP_ */
