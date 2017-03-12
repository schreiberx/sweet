/*
 * ErrorCheck.hpp
 *
 *  Created on: 3 Feb 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_ERRORCHECK_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_ERRORCHECK_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/Convert_SphereDataComplex_to_SphereData.hpp>
#include <sweet/sphere/Convert_SphereData_to_SphereDataComplex.hpp>


class ErrorCheck
{
public:
	static
	bool check(
			const SphereData &i_lhs,
			const SphereData &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		const SphereData lhs = i_lhs;
		const SphereData rhs = i_rhs;

		SphereData diff = i_lhs.physical_diff_realconst(rhs);

		double lhs_maxabs = lhs.physical_reduce_max_abs();
		double rhs_maxabs = rhs.physical_reduce_max_abs();

		double normalize_fac = 1.0;

		if (i_normalization)
		{
			normalize_fac = std::max(lhs_maxabs, rhs_maxabs);

			if (normalize_fac < i_error_threshold)
			{
				std::cout << "Normalization for '" << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				normalize_fac = 1.0;
			}
		}

		double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.physical_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\terror threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			lhs.physical_file_write("o_error_lhs_values.csv");
			rhs.physical_file_write("o_error_rhs_values.csv");
			(lhs-rhs).physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			diff.physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			if (i_ignore_error)
				std::cerr << "Error ignored" << std::endl;
			else
				FatalError("Error too large");

			return true;
		}
		return false;
	}



public:
	static
	bool check(
			const SphereDataComplex &i_lhs,
			const SphereDataComplex &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
//		SphereData diff = i_lhs - i_rhs;
		SphereDataComplex diff = i_lhs.physical_diff_realconst(i_rhs);

		double lhs_maxabs = SphereDataComplex(i_lhs).physical_reduce_max_abs();
		double rhs_maxabs = SphereDataComplex(i_rhs).physical_reduce_max_abs();

		double normalize_fac = 1.0;

		if (i_normalization)
		{
			normalize_fac = std::max(lhs_maxabs, rhs_maxabs);

			if (normalize_fac < i_error_threshold)
			{
				std::cout << "Normalization for '" << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				normalize_fac = 1.0;
			}
		}

		double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.physical_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\terror threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(i_lhs).physical_file_write("o_error_lhs_values.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(i_rhs).physical_file_write("o_error_rhs_values.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(i_lhs-i_rhs).physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(diff).physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			if (i_ignore_error)
				std::cerr << "Error ignored" << std::endl;
			else
				FatalError("Error too large");

			return true;
		}
		return false;
	}



public:
	static
	bool checkTruncated(
			const SphereData &i_lhs,
			const SphereData &i_rhs,
			const SphereDataConfig *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold,	// = 1.0,
			double i_ignore_error,		// = false,
			bool i_normalization		// = true
	)
	{
		SphereData lhsr = SphereData(i_lhs).spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereData rhsr = SphereData(i_rhs).spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereData diff = SphereData(lhsr).physical_diff_realconst(SphereData(rhsr));

		double lhs_maxabs = SphereData(lhsr).physical_reduce_max_abs();
		double rhs_maxabs = SphereData(rhsr).physical_reduce_max_abs();

		double normalize_fac = std::min(lhs_maxabs, rhs_maxabs);

		if (std::max(lhs_maxabs, rhs_maxabs) < i_error_threshold)
		{
			std::cout << "Error computation for '" << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
			return false;
		}


		double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.physical_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tError threshold: " << i_error_threshold << "\tNormalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
				return true;
			}

			lhsr.physical_file_write("o_error_lhs.csv");
			rhsr.physical_file_write("o_error_rhs.csv");
			(lhsr-rhsr).physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			diff.physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			FatalError("Error too large");
			return true;
		}
		return false;
	}


public:
	static
	bool checkTruncated(
			const SphereDataComplex &i_lhs,
			const SphereDataComplex &i_rhs,
			const SphereDataConfig *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		SphereDataComplex lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereDataComplex rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereDataComplex diff = lhsr.physical_diff_realconst(rhsr);

		double normalize_fac;

		if (i_normalization)
		{
			double lhs_maxabs = lhsr.physical_reduce_max_abs();
			double rhs_maxabs = rhsr.physical_reduce_max_abs();

			if (std::max(lhs_maxabs, rhs_maxabs) < i_error_threshold)
			{
				std::cout << "Error for " << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				return false;
			}

			normalize_fac = std::min(lhsr.physical_reduce_max_abs(), rhsr.physical_reduce_max_abs());

			if (normalize_fac == 0)
			{
				std::cout << "Error for " << i_id << "' ignored since at least one field is Zero" << std::endl;
				return false;
			}
		}
		else
		{
			normalize_fac = 1.0;
		}

		double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.physical_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tError threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
				return false;
			}

			Convert_SphereDataComplex_To_SphereData::physical_convert_real(lhsr).physical_file_write("o_error_lhs.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(rhsr).physical_file_write("o_error_rhs.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(lhsr-rhsr).physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(diff).physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			FatalError("Error too large");
			return true;
		}
		return false;
	}



public:
	static
	bool checkTruncatedSpectral(
			const SphereDataComplex &i_lhs,
			const SphereDataComplex &i_rhs,
			const SphereDataConfig *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		SphereDataComplex lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereDataComplex rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereDataComplex diff = lhsr-rhsr;

		double normalize_fac = 1.0;

		if (i_normalization)
		{
			double lhs_maxabs = lhsr.physical_reduce_max_abs();
			double rhs_maxabs = rhsr.physical_reduce_max_abs();

			if (std::max(lhs_maxabs, rhs_maxabs) < i_error_threshold)
			{
				std::cout << "Error for " << i_id << "' ignored since both fields are below threshold tolerance" << std::endl;
				return false;
			}

			normalize_fac = std::min(lhsr.physical_reduce_max_abs(), rhsr.physical_reduce_max_abs());

			if (normalize_fac == 0)
			{
				std::cout << "Error for " << i_id << "' ignored since at least one field is Zero" << std::endl;
				return false;
			}
		}

		double rel_max_abs = diff.physical_reduce_max_abs() / normalize_fac;
		double rel_rms = diff.physical_reduce_rms() / normalize_fac;

		std::cout << "Error for " << i_id << ": \t" << rel_max_abs << "\t" << rel_rms << "\t\tError threshold: " << i_error_threshold << " with normalization factor " << normalize_fac << std::endl;

		if (rel_max_abs > i_error_threshold)
		{
			if (i_ignore_error)
			{
				std::cerr << "Error ignored (probably because extended modes not >= 2)" << std::endl;
				return false;
			}

			Convert_SphereDataComplex_To_SphereData::physical_convert_real(lhsr).physical_file_write("o_error_lhs.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(rhsr).physical_file_write("o_error_rhs.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(lhsr-rhsr).physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			Convert_SphereDataComplex_To_SphereData::physical_convert_real(diff).physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			FatalError("Error too large");
			return true;
		}
		return false;
	}

};


#endif /* SRC_INCLUDE_SWEET_SPHERE_ERRORCHECK_HPP_ */
