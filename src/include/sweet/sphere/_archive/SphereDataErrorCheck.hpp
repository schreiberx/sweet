/*
 * SphereDataErrorCheck.hpp
 *
 *  Created on: 3 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE__ARCHIVE_SPHEREDATAERRORCHECK_HPP_
#define SRC_INCLUDE_SWEET_SPHERE__ARCHIVE_SPHEREDATAERRORCHECK_HPP_

#include <sweet/sphere/Convert_SphereDataSpectral_to_SphereDataSpectralComplex.hpp>
#include <sweet/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>


class SphereDataErrorCheck
{
public:
	static
	bool check(
			const SphereData_Spectral &i_lhs,
			const SphereData_Spectral &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		const SphereData_Spectral lhs = i_lhs;
		const SphereData_Spectral rhs = i_rhs;

		SphereData_Physical diff = i_lhs.getSphereDataPhysical()-rhs.getSphereDataPhysical();

		double lhs_maxabs = lhs.getSphereDataPhysical().physical_reduce_max_abs();
		double rhs_maxabs = rhs.getSphereDataPhysical().physical_reduce_max_abs();

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
			if (i_ignore_error)
			{
				std::cerr << "Error ignored" << std::endl;
			}
			else
			{
				lhs.getSphereDataPhysical().physical_file_write("o_error_lhs_values.csv");
				rhs.getSphereDataPhysical().physical_file_write("o_error_rhs_values.csv");
				(lhs-rhs).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
				diff.physical_file_write("o_error_lhs_rhs_diff_physical.csv");

				FatalError("Error too large");
			}

			return true;
		}
		return false;
	}



public:
	static
	bool check(
			const SphereData_SpectralComplex &i_lhs,
			const SphereData_SpectralComplex &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
//		SphereDataSpectral diff = i_lhs - i_rhs;
		SphereData_SpectralComplex diff = i_lhs.physical_diff_realconst(i_rhs);

		double lhs_maxabs = SphereData_SpectralComplex(i_lhs).physical_reduce_max_abs();
		double rhs_maxabs = SphereData_SpectralComplex(i_rhs).physical_reduce_max_abs();

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
			if (i_ignore_error)
			{
				std::cerr << "Error ignored" << std::endl;
			}
			else
			{
				Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(i_lhs).getSphereDataPhysical().physical_file_write("o_error_lhs_values.csv");
				Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(i_rhs).getSphereDataPhysical().physical_file_write("o_error_rhs_values.csv");
				Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(i_lhs-i_rhs).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
				Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(diff).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_physical.csv");

				FatalError("Error too large");
			}

			return true;
		}
		return false;
	}



public:
	static
	bool checkTruncated(
			const SphereData_Spectral &i_lhs,
			const SphereData_Spectral &i_rhs,
			const SphereData_Config *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold,	// = 1.0,
			double i_ignore_error,		// = false,
			bool i_normalization		// = true
	)
	{
		SphereData_Spectral lhsr = SphereData_Spectral(i_lhs).spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereData_Spectral rhsr = SphereData_Spectral(i_rhs).spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereData_Physical diff = lhsr.getSphereDataPhysical()-rhsr.getSphereDataPhysical();

		double lhs_maxabs = SphereData_Spectral(lhsr).getSphereDataPhysical().physical_reduce_max_abs();
		double rhs_maxabs = SphereData_Spectral(rhsr).getSphereDataPhysical().physical_reduce_max_abs();

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

			lhsr.getSphereDataPhysical().physical_file_write("o_error_lhs.csv");
			rhsr.getSphereDataPhysical().physical_file_write("o_error_rhs.csv");
			(lhsr-rhsr).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			diff.physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			FatalError("Error too large");
			return true;
		}
		return false;
	}


public:
	static
	bool checkTruncated(
			const SphereData_SpectralComplex &i_lhs,
			const SphereData_SpectralComplex &i_rhs,
			const SphereData_Config *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		SphereData_SpectralComplex lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereData_SpectralComplex rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereData_SpectralComplex diff = lhsr.physical_diff_realconst(rhsr);

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

			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(lhsr).getSphereDataPhysical().physical_file_write("o_error_lhs.csv");
			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(rhsr).getSphereDataPhysical().physical_file_write("o_error_rhs.csv");
			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(lhsr-rhsr).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(diff).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			FatalError("Error too large");
			return true;
		}
		return false;
	}



public:
	static
	bool checkTruncatedSpectral(
			const SphereData_SpectralComplex &i_lhs,
			const SphereData_SpectralComplex &i_rhs,
			const SphereData_Config *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		SphereData_SpectralComplex lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereData_SpectralComplex rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereData_SpectralComplex diff = lhsr-rhsr;

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

			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(lhsr).getSphereDataPhysical().physical_file_write("o_error_lhs.csv");
			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(rhsr).getSphereDataPhysical().physical_file_write("o_error_rhs.csv");
			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(lhsr-rhsr).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_spectral.csv");
			Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(diff).getSphereDataPhysical().physical_file_write("o_error_lhs_rhs_diff_physical.csv");

			FatalError("Error too large");
			return true;
		}
		return false;
	}

};


#endif /* SRC_INCLUDE_SWEET_SPHERE__ARCHIVE_SPHEREDATAERRORCHECK_HPP_ */
