/*
 * test_sph_quadrature_nodes.hpp
 *
 *  Created on: 3 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_TESTSPH_QUADRATURE_NODES_HPP_
#define SRC_TESTSPH_QUADRATURE_NODES_HPP_


#include <sweet/SimulationVariables.hpp>
#include <libmath/BandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>
//#include <sweet/sphere/SphereOperators_SphereData.hpp>





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



SimulationVariables simVars;

void run_tests(
		SphereData_Config *sphereDataConfig
)
{
	double epsilon = 1e-12;
	epsilon *= (sphereDataConfig->spectral_modes_n_max);
	std::cout << "Using max allowed error of " << epsilon << std::endl;

	std::cout << std::setprecision(10);

	{
		SphereData_Physical physical(sphereDataConfig);
		physical.physical_update_lambda_cogaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(2.0)*0.5;
			}
		);

		SphereData_Spectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 0 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereDataErrorCheck::check(SphereData_Spectral(physical), spectral, "n=0, m=0", epsilon, false, true);
	}

	{
		SphereData_Physical physical(sphereDataConfig);
		physical.physical_update_lambda_gaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(6.0)*mu*0.5;
			}
		);

		SphereData_Spectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 1 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereData_Spectral(physical).spectral_print();
		spectral.spectral_print();

		SphereDataErrorCheck::check(SphereData_Spectral(physical), spectral, "n=1, m=0", epsilon, false, true);
	}

	{
		SphereData_Physical physical(sphereDataConfig);
		physical.physical_update_lambda_gaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(10.0)/4.0 * (3.0*mu*mu - 1.0);
			}
		);

		SphereData_Spectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 2 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereDataErrorCheck::check(SphereData_Spectral(physical), spectral, "n=2, m=0", epsilon, false, true);
	}


	{
		SphereData_Physical physical(sphereDataConfig);
		physical.physical_update_lambda_gaussian_grid(
			[&](double lat, double mu, double &io_data)
			{
				io_data = std::sqrt(14.0)/4.0*mu * (5.0*mu*mu - 3.0);
			}
		);

		SphereData_Spectral spectral(sphereDataConfig);
		spectral.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &io_data)
			{
				io_data = (n == 3 && m == 0);

				// rescale because of 2pi FT along longitude
				io_data *= sqrt(2.0*M_PI);
			}
		);

		SphereDataErrorCheck::check(SphereData_Spectral(physical), spectral, "n=3, m=0", epsilon, false, true);
	}
}




int main(
		int i_argc,
		char *const i_argv[]
)
{
	MemBlockAlloc numaBlockAlloc;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.disc.space_res_spectral[0] == 0)
		FatalError("Set number of spectral modes to use SPH!");

	SphereData_Config sphereDataConfig;
	sphereDataConfig.setupAutoPhysicalSpace(
					simVars.disc.space_res_spectral[0],
					simVars.disc.space_res_spectral[1],
					&simVars.disc.space_res_physical[0],
					&simVars.disc.space_res_physical[1],
					simVars.misc.reuse_spectral_transformation_plans
			);

	run_tests(&sphereDataConfig);
}



#endif
