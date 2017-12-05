/*
 * RexiNG.hpp
 *
 *  Created on: 3 July 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_INCLUDE_REXI_FILE_HPP_
#define SRC_INCLUDE_REXI_FILE_HPP_

#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <sweet/FatalError.hpp>
#include <libmath/DQStuff.hpp>
#include <rexi/REXI_File_Coefficients.hpp>
#include <rexi/REXI_File.hpp>




/**
 * Next generation REXI implementation.
 *
 * This class provides the abstraction layer for a
 * Rational approximation of EXponential Integrators.
 *
 * In contrast to the REXI.hpp implementation which is
 * based on [Terry et al. paper], this work is not based on
 * an rational approximation of a Gaussian basis function
 * which is used to approximate one phiN function.
 *
 * 	phi0(x) = e^{ix}
 * 	phi1(x) = (e^{ix} - 1)/(ix)
 * 	phi2(x) = (e^{ix} - ix - 1)/(ix)^2
 *
 * This implementation uses a database of coefficients which
 * were computed as a preprocessing step. These coefficients
 *  - result in higher accuracy (at least in the test cases)
 *    and
 *  - require less poles for the approximation
 *
 * This is based on the "FAF" Python code
 * (Function approximation of a function)
 *
 */
template <typename T = double>
class REXI_File
{
public:
	typedef std::complex<T> complexT;

	bool reduce_to_half;
	std::string faf_data_directory;

	std::string filename;

public:
	std::vector<complexT> alpha;
	std::vector<complexT> beta;
//	std::vector<complexT> beta_im;

	REXI_File_Coefficients<T> fafcoeffs;


	static
	constexpr
	T None()
	{
		return std::numeric_limits<T>::quiet_NaN();
	}

	static
	constexpr
	bool isNone(T i_value)
	{
#if __GNUC__ < 6
		return isnan(i_value);
#else
		return std::isnan(i_value);
#endif
	}


#if 0
public:
	bool auto_load(
			const std::string &i_function_name,		///< "phi0", "phi1", "phi2", etc.

			int i_N = 0,						///< Number of approximation poles

			T i_max_error = None(),
			T i_max_error_double_precision = None(),
			T i_test_min = None(),
			T i_test_max = None(),
			T i_basis_function_scaling = None(),
			T i_basis_function_spacing = None(),
			T i_basis_function_rat_shift = None(),

			bool i_reduce_to_half = false,
			const std::string &i_faf_data_directory = "data/faf_data/"
	)
	{
		/// fix for default program parameter value of h
		if (i_basis_function_spacing <= 0)
			i_basis_function_spacing = None();

		std::string faf_data_dir = i_faf_data_directory + "/faf_data_rationalcplx_"+i_function_name;
		REXI_File_Coefficients<T> target_fafcoeffs;

		target_fafcoeffs.N = i_N;
		target_fafcoeffs.max_error = i_max_error;
		target_fafcoeffs.max_error_double_precision = i_max_error_double_precision;

		target_fafcoeffs.test_min = i_test_min;
		target_fafcoeffs.test_max = i_test_max;

		target_fafcoeffs.basis_function_scaling = i_basis_function_scaling;
		target_fafcoeffs.basis_function_spacing = i_basis_function_spacing;
		target_fafcoeffs.basis_function_rat_shift = i_basis_function_rat_shift;

		reduce_to_half = i_reduce_to_half;
		faf_data_directory = i_faf_data_directory;

		DIR *dirp = opendir(faf_data_dir.c_str());

		if (dirp == nullptr)
			FatalError(std::string("Unable to open directory ") + i_faf_data_directory);

		bool best_found = false;
		REXI_File_Coefficients<T> &best = fafcoeffs;

		struct dirent *dp;
		while (dirp)
		{
			errno = 0;	/// required! set errno to zero!
			if ((dp = readdir(dirp)) == nullptr)
			{
				if (errno == 0)
					break;

				FatalError(std::string("Error while reading faf coefficients from directory '")+faf_data_dir+"'");
			}

			std::string filename = dp->d_name;

			if (filename == ".")
				continue;

			if (filename == "..")
				continue;

			// skip hidden files
			if (filename.data()[0] == '.')
				continue;

			std::string filepath = faf_data_dir + "/"+ filename;

			REXI_File_Coefficients<T> test_fafcoeffs;
			test_fafcoeffs.load_from_file(filepath);

			T eps = 1e-10;

			if (target_fafcoeffs.N != 0)
				if (test_fafcoeffs.N != target_fafcoeffs.N)
					continue;


			if (!isNone(target_fafcoeffs.max_error))
				if (test_fafcoeffs.max_error > target_fafcoeffs.max_error)
					continue;

			if (!isNone(target_fafcoeffs.max_error_double_precision))
				if (test_fafcoeffs.max_error_double_precision > target_fafcoeffs.max_error_double_precision)
					continue;

			if (!isNone(target_fafcoeffs.test_min))
				if (test_fafcoeffs.test_min > target_fafcoeffs.test_min)
					continue;

			if (!isNone(target_fafcoeffs.test_max))
				if (test_fafcoeffs.test_max < target_fafcoeffs.test_max)
					continue;

			if (target_fafcoeffs.function_name != "")
				if (test_fafcoeffs.function_name != target_fafcoeffs.function_name)
					continue;
#if 0
			if (!isNone(target_fafcoeffs.function_scaling))
				if (std::abs(test_fafcoeffs.function_scaling-target_fafcoeffs.function_scaling) > eps)
					continue;

			if (!isNone(target_fafcoeffs.basis_function_name))
				if (test_fafcoeffs.basis_function_name != target_fafcoeffs.basis_function_name)
					continue;
#endif

			if (!isNone(target_fafcoeffs.basis_function_spacing))
				if (std::abs(test_fafcoeffs.basis_function_spacing-target_fafcoeffs.basis_function_spacing) > eps)
					continue;


			if (!isNone(target_fafcoeffs.basis_function_scaling))
				if (std::abs(test_fafcoeffs.basis_function_scaling-target_fafcoeffs.basis_function_scaling) > eps)
					continue;

			if (!isNone(target_fafcoeffs.basis_function_rat_shift))
				if (std::abs(test_fafcoeffs.basis_function_rat_shift-target_fafcoeffs.basis_function_rat_shift) > eps)
					continue;


			/*
			 * Soft constraints next
			 */

			if (	isNone(target_fafcoeffs.max_error) &&
					isNone(target_fafcoeffs.max_error_double_precision) &&
					target_fafcoeffs.N == 0
			)
			{
				if (!best_found)
				{
					best = test_fafcoeffs;
					best_found = true;
					continue;
				}

				// If nothing is specified, search for highest accuracy!
				if (best.max_error > test_fafcoeffs.max_error)
				{
					best = test_fafcoeffs;
					best_found = true;
					continue;
				}

				continue;
			}

			if (target_fafcoeffs.N == 0)
			{
				if (!best_found)
				{
					best = test_fafcoeffs;
					best_found = true;
					continue;
				}

				// Minimize number of poles
				if (best.N < test_fafcoeffs.N)
					continue;

				// If poles are the same, try to minimize error
				if (best.N == test_fafcoeffs.N)
				{
					if (best.max_error < test_fafcoeffs.max_error)
					{
						if (!isNone(best.max_error_double_precision) && !isNone(test_fafcoeffs.max_error_double_precision))
						{
							if (best.max_error_double_precision < test_fafcoeffs.max_error_double_precision)
								continue;
						}
						else
							continue;
					}
				}

				best = test_fafcoeffs;
				best_found = true;
				continue;
			}


			if (true)	// max_error == None
			{
				// Hard constraint of test_fafcoeffs.N == N was already triggered before
				if (!best_found)
				{
					best = test_fafcoeffs;
					best_found = true;
					continue;
				}

				// Minimize error
				if (best.max_error > test_fafcoeffs.max_error)
					continue;

				// Minimize error
				if (!isNone(best.max_error_double_precision) && !isNone(test_fafcoeffs.max_error_double_precision))
					if (best.max_error_double_precision > test_fafcoeffs.max_error_double_precision)
						continue;

				best = test_fafcoeffs;
				best_found = true;
				continue;
			}
		}

		closedir(dirp);

		if (!best_found)
			return false;
			//FatalError("RexiNG: No coefficients found for given constraints");


		return true;
	}
#endif

	bool load_from_file(
			std::string& i_filename_alpha,
			std::string& i_filename_beta
	)
	{
		bool retval;

		retval = fafcoeffs.load_from_file(i_filename_alpha, alpha);
		if (!retval)
			return false;

		retval = fafcoeffs.load_from_file(i_filename_beta, beta);
		if (!retval)
			return false;

		return true;
	}
#if 0
	void generate_alpha_and_beta(bool i_reduce_to_half)
	{
		int N = fafcoeffs.N;

		alpha.resize(N);
		beta.resize(N);

		for (int i = 0; i < N; i++)
		{
			// generate alphas
			int K = i - (N/2);
			alpha[i] = std::complex<T>(fafcoeffs.basis_function_rat_shift/fafcoeffs.basis_function_scaling, -(T)K*fafcoeffs.basis_function_spacing);

			// generate betas
			beta[i] = fafcoeffs.weights_cplx[i];

			beta[i] /= fafcoeffs.basis_function_scaling;
		}

		if ((N & 1) != 1)
			FatalError("N must be odd!");

		if (i_reduce_to_half)
		{
			int newN = N/2+1;
			for (int i = 0; i < newN-1; i++)
				beta[i] += conj(beta[N-1-i]);

			alpha.resize(newN);
			beta.resize(newN);
		}
	}
#endif

#if 0
	void output()
	{
		int N = alpha.size();
		std::cout << "N: " << N << std::endl;

		for (int i = 0; i < N; i++)
			std::cout << "alpha[" << i << "] = " << alpha[i] << std::endl;

		for (int i = 0; i < N; i++)
			std::cout << "beta[" << i << "] = " << beta[i] << std::endl;
	}



	/**
	 * \return \f$ cos(x) + i*sin(x) \f$
	 */
	T eval_real(
			T i_x	///< sampling position
	)
	{
		double eps = 1e-8;

		if (fafcoeffs.function_name == "phi0")
			return DQStuff::Re(DQStuff::expIm(i_x));

		if (fafcoeffs.function_name == "phi1")
		{
			if (std::abs(i_x) < eps)
				return DQStuff::Re(DQStuff::expIm(i_x));
			else
				return DQStuff::Re((DQStuff::expIm(i_x)-1.0)/DQStuff::I(i_x));
		}

		if (fafcoeffs.function_name == "phi2")
		{
			if (std::abs(i_x) < eps)
				return (T)0.5;
			else
				return DQStuff::Re((DQStuff::expIm(i_x) - 1.0 - DQStuff::I(i_x))/(-i_x*i_x));
		}

		FatalError("Function "+fafcoeffs.function_name+" not supported");
		return 0;
	}



	/**
	 * \return \f$ f(x) \approx Re(cos(x) + i*sin(x)) = cos(x) \f$
	 */
	T approx_real(
			T i_x
	)
	{
		T sum = 0;

		std::size_t S = alpha.size();

		for (std::size_t n = 0; n < S; n++)
			sum += (beta[n] / (complexT(0.0, i_x) + alpha[n])).real();

		return sum;
	}
#endif
};


#endif
