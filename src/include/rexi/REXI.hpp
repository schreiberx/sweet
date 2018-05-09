/*
 * REXI_OLD_and_NG.hpp
 *
 *  Created on: 4 Aug 2017
 *      Author: martin
 */

#ifndef SRC_INCLUDE_REXI_HPP_
#define SRC_INCLUDE_REXI_HPP_


#include <rexi/REXI_SimulationVariables.hpp>
#include <rexi/REXI_Terry.hpp>
#include <rexi/REXI_File.hpp>
#include <rexi/REXI_CI.hpp>
#include <vector>
#include <complex>


#if SWEET_MPI
#include <mpi.h>
#endif

/**
 * Interface to load either REXI via Terrys method or via File
 */
class REXI
{
#if SWEET_MPI
	static
	int& getMPIRank()
	{
		static int mpi_rank;
		return mpi_rank;
	}
#endif



public:
	REXI()
	{
#if SWEET_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &getMPIRank());
#endif
	}

public:
	static
	void load(
			REXI_SimulationVariables *i_rexiSimVars,
			const std::string &i_function_name,

			std::vector<std::complex<double>> &alpha,
			std::vector<std::complex<double>> &beta,

			double i_dt,	///< Timestep size, this gets important for the Gaussian filter

			int i_verbosity
	)
	{
		alpha.clear();
		beta.clear();

		if (i_rexiSimVars->rexi_method == "file")
		{
			/// REXI next generation stuff
			REXI_File<> rexi_file;

			bool retval = rexi_file.load_from_file(i_rexiSimVars->file_filename_alpha, i_rexiSimVars->file_filename_beta);
#if 0
			if (i_rexiSimVars->file_filename_alpha == "")
			{
				retval = rexi_file.auto_load(
						i_function_name,
						i_rexiSimVars->file_N,	/// N
						rexi_file.None(),			/// max_error
						i_rexiSimVars->file_max_error_double_precision,			/// max_error_double_precision
						i_rexiSimVars->file_test_min,
						i_rexiSimVars->file_test_max,
						rexi_file.None(),			/// basis_function_scaling
						i_rexiSimVars->file_h, 	/// basis_function_spacing
						rexi_file.None(),			/// basis_function rat shift

						i_rexiSimVars->use_half_poles,
						i_rexiSimVars->file_faf_dir
					);
			}
			else
			{
				retval = rexi_file.load_from_file(i_rexiSimVars->file_filename, i_rexiSimVars->use_half_poles);
			}
#endif
			if (!retval)
				FatalError(std::string("Not able to find coefficients for given constraints for function "+i_function_name));

//			if (i_verbosity > 0)
//				std::cout << "Loaded REXI coefficients from file '" << rexi_file.fafcoeffs.filename << "'" << std::endl;

			alpha = rexi_file.alpha;
			beta = rexi_file.beta;
		}
		else if (i_rexiSimVars->rexi_method == "terry")
		{
			/// REXI stuff
			REXI_Terry<> rexi_terry;
			rexi_terry.setup(i_function_name, i_rexiSimVars->terry_h, i_rexiSimVars->terry_M, i_rexiSimVars->terry_L, i_rexiSimVars->use_half_poles, i_rexiSimVars->normalization);

			alpha = rexi_terry.alpha;
			beta = rexi_terry.beta;
		}
		else if (i_rexiSimVars->rexi_method == "ci")
		{
			/// REXI stuff
			REXI_CI<> rexi_ci;

			if (i_rexiSimVars->ci_max_real >= 0)
			{
				rexi_ci.setup_shifted_circle(
						i_function_name,
						i_rexiSimVars->ci_n, i_rexiSimVars->ci_max_real, i_rexiSimVars->ci_max_imag,
						i_rexiSimVars->ci_gaussian_filter_scale_a, i_rexiSimVars->ci_gaussian_filter_dt_norm, i_rexiSimVars->ci_gaussian_filter_exp_N, i_dt
					);
			}
			else
			{
				rexi_ci.setup(
						i_function_name,
						i_rexiSimVars->ci_n, i_rexiSimVars->ci_primitive, i_rexiSimVars->ci_s_real, i_rexiSimVars->ci_s_imag, i_rexiSimVars->ci_mu,
						i_rexiSimVars->ci_gaussian_filter_scale_a, i_rexiSimVars->ci_gaussian_filter_dt_norm, i_rexiSimVars->ci_gaussian_filter_exp_N, i_dt
					);
			}

			alpha = rexi_ci.alpha;
			beta = rexi_ci.beta;
		}
		else if (i_rexiSimVars->rexi_method == "direct")
		{
			// no REXI
		}
		else
		{
			FatalError("REXI Mode not supported");
		}


		if (i_rexiSimVars->beta_cutoff != 0)
		{
			std::vector<std::complex<double>>::iterator ia = alpha.begin();
			std::vector<std::complex<double>>::iterator ib = beta.begin();
			while (ia != alpha.end())
			{
				if (std::abs(*ib) < i_rexiSimVars->beta_cutoff)
				{
					ia = alpha.erase(ia);
					ib = beta.erase(ib);
				}
				else
				{
					ia++;
					ib++;
				}
			}
		}

#if SWEET_MPI
		if (getMPIRank() == 0)
#endif
		if (i_verbosity > 2)
		{
			int N = alpha.size();
			std::cout << "N: " << N << std::endl;

	//		std::cout << "Alpha:" << std::endl;
			for (int i = 0; i < N; i++)
				std::cout << "alpha[" << i << "] = " << alpha[i] << std::endl;

	//		std::cout << "Beta:" << std::endl;
			for (int i = 0; i < N; i++)
				std::cout << "beta[" << i << "] = " << beta[i] << std::endl;
		}
	}


	static
	void testREXIphi0(
			std::vector<std::complex<double>> &alpha,
			std::vector<std::complex<double>> &beta_re,
			double i_test_abs,
			double i_max_error = 1e-10
	)
	{
		/*
		 * WARNING
		 * WARNING
		 * WARNING
		 *
		 * This only works without halving the poles!!!
		 *
		 * WARNING
		 * WARNING
		 * WARNING
		 */
		std::cout << "************************************************************" << std::endl;
		std::cout << "* RUNNING TESTS " << std::endl;
		std::cout << "************************************************************" << std::endl;

		typedef double T;

#if 1

		double start = -M_PI;//*0.5;
		double end = M_PI;//*0.5;
		//double step_size = 1e-2;
		double step_size = M_PI/80.0;//1e-1;

#else

		double start = -1.0;
		double end = 1.0;
		double step_size = 1e-4;
#endif

		for (double x = start; x < end; x += step_size)
		{
			T correct_real = std::exp(std::complex<double>(0, 1)*x).real();

			std::complex<T> approx = 0;
			for (std::size_t n = 0; n < alpha.size(); n++)
				approx += (beta_re[n] / (std::complex<double>(0.0, x) + alpha[n]));
			T approx_real = approx.real();

			if (DQStuff::abs(approx_real) > 1.0+1e-12)
			{
				std::cerr << "approx value " << approx_real << " not bounded by unity (just a warning and not a problem) at x=" << x << std::endl;
				std::cerr << "correct: " << correct_real << std::endl;
				std::cerr << "approx: " << approx_real << std::endl;
			}
			T error_real = DQStuff::abs(correct_real - approx_real);

			if (error_real > i_max_error && 0)
			{
				i_max_error = error_real;
				std::cout << "ERROR " << error_real << " too large at " << x << std::endl;
			}
		}
	}

};


#endif /* SRC_INCLUDE_REXI_HPP_ */
