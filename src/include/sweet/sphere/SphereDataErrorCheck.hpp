/*
 * SphereDataErrorCheck.hpp
 *
 *  Created on: 3 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_SPHEREDATAERRORCHECK_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_SPHEREDATAERRORCHECK_HPP_

#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/Convert_SphereDataSpectral_to_SphereDataSpectralComplex.hpp>
#include <sweet/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/sphere/SphereDataSpectralComplex.hpp>


class SphereDataErrorCheck
{
public:
	static
	bool check(
			const SphereDataSpectral &i_lhs,
			const SphereDataSpectral &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		const SphereDataSpectral lhs = i_lhs;
		const SphereDataSpectral rhs = i_rhs;

		SphereDataPhysical diff = i_lhs.getSphereDataPhysical()-rhs.getSphereDataPhysical();

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
			const SphereDataSpectralComplex &i_lhs,
			const SphereDataSpectralComplex &i_rhs,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
//		SphereDataSpectral diff = i_lhs - i_rhs;
		SphereDataSpectralComplex diff = i_lhs.physical_diff_realconst(i_rhs);

		double lhs_maxabs = SphereDataSpectralComplex(i_lhs).physical_reduce_max_abs();
		double rhs_maxabs = SphereDataSpectralComplex(i_rhs).physical_reduce_max_abs();

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
			const SphereDataSpectral &i_lhs,
			const SphereDataSpectral &i_rhs,
			const SphereDataConfig *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold,	// = 1.0,
			double i_ignore_error,		// = false,
			bool i_normalization		// = true
	)
	{
		SphereDataSpectral lhsr = SphereDataSpectral(i_lhs).spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereDataSpectral rhsr = SphereDataSpectral(i_rhs).spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereDataPhysical diff = lhsr.getSphereDataPhysical()-rhsr.getSphereDataPhysical();

		double lhs_maxabs = SphereDataSpectral(lhsr).getSphereDataPhysical().physical_reduce_max_abs();
		double rhs_maxabs = SphereDataSpectral(rhsr).getSphereDataPhysical().physical_reduce_max_abs();

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
			const SphereDataSpectralComplex &i_lhs,
			const SphereDataSpectralComplex &i_rhs,
			const SphereDataConfig *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		SphereDataSpectralComplex lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereDataSpectralComplex rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereDataSpectralComplex diff = lhsr.physical_diff_realconst(rhsr);

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
			const SphereDataSpectralComplex &i_lhs,
			const SphereDataSpectralComplex &i_rhs,
			const SphereDataConfig *i_sphereDataConfig,
			const std::string &i_id,
			double i_error_threshold = 1.0,
			bool i_ignore_error = false,
			bool i_normalization = true
	)
	{
		SphereDataSpectralComplex lhsr = i_lhs.spectral_returnWithDifferentModes(i_sphereDataConfig);
		SphereDataSpectralComplex rhsr = i_rhs.spectral_returnWithDifferentModes(i_sphereDataConfig);

		SphereDataSpectralComplex diff = lhsr-rhsr;

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


#endif /* SRC_INCLUDE_SWEET_SPHERE_SPHEREDATAERRORCHECK_HPP_ */
