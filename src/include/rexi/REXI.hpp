/*
 * REXI.hpp
 *
 *  Created on: 4 Aug 2017
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_REXI_HPP_
#define SRC_INCLUDE_REXI_HPP_


#include <rexi/REXI_SimulationVariables.hpp>
#include <rexi/REXI_Terry.hpp>
#include <rexi/REXI_File.hpp>
#include <rexi/REXI_CI.hpp>
#include <rexi/REXICoefficients.hpp>
#include <vector>
#include <complex>


/**
 * Interface to load REXI coefficients via various methods
 */
template <typename T=double>
class REXI
{

public:
	REXI()
	{
	}



public:
	static
	void load(
			REXI_SimulationVariables *i_rexiSimVars,
			const std::string &i_function_name,

			std::vector<std::complex<T>> &o_alpha,
			std::vector<std::complex<T>> &o_beta,
			std::complex<T> &o_gamma,

			int i_verbosity
	)
	{
		o_alpha.clear();
		o_beta.clear();

		if (i_rexiSimVars->rexi_method == "terry")
		{
			/// REXI stuff
			REXI_Terry<T, T> rexi_terry;
			rexi_terry.setup(i_function_name, i_rexiSimVars->terry_h, i_rexiSimVars->terry_M, i_rexiSimVars->terry_L, i_rexiSimVars->terry_reduce_to_half, i_rexiSimVars->terry_normalization);

			o_alpha = rexi_terry.alpha;
			o_beta = rexi_terry.beta;
		}
		else if (i_rexiSimVars->rexi_method == "ci")
		{
			/// REXI stuff
			REXI_CI<T, T> rexi_ci;

			if (i_rexiSimVars->ci_max_real >= 0)
			{
				rexi_ci.setup_shifted_circle(
						i_function_name,
						i_rexiSimVars->ci_n, i_rexiSimVars->ci_max_real, i_rexiSimVars->ci_max_imag
					);
			}
			else
			{
				rexi_ci.setup(
						i_function_name,
						i_rexiSimVars->ci_n, i_rexiSimVars->ci_primitive, i_rexiSimVars->ci_s_real, i_rexiSimVars->ci_s_imag, i_rexiSimVars->ci_mu
					);
			}

			o_alpha = rexi_ci.alpha;
			o_beta = rexi_ci.beta;
		}
		else if (i_rexiSimVars->rexi_method == "direct")
		{
			// no REXI
		}
		else
		{
			FatalError("REXI Mode not supported");
		}


		if (i_verbosity > 2)
		{
			int N = o_alpha.size();
			std::cout << "N: " << N << std::endl;

			for (int i = 0; i < N; i++)
				std::cout << "alpha[" << i << "] = " << (double)o_alpha[i].real() << ", " << (double)o_alpha[i].imag() << std::endl;

			for (int i = 0; i < N; i++)
				std::cout << "beta[" << i << "] = " << (double)o_beta[i].real() << ", " << (double)o_beta[i].imag() << std::endl;
		}
	}


public:
	static
	bool load(
			REXI_SimulationVariables *i_rexiSimVars,
			const std::string &i_function_name,

			REXICoefficients<T> &o_rexiCoefficients,

			int i_verbosity
	)
	{
		if (i_rexiSimVars->rexi_method == "file")
		{
			for (
				std::vector<REXI_SimulationVariables::REXIFile>::iterator iter = i_rexiSimVars->p_rexi_files_processed.begin();
				iter != i_rexiSimVars->p_rexi_files_processed.end();
				iter++
			)
			{
				if (iter->function_name == i_function_name)
					return o_rexiCoefficients.load_from_file(iter->filename);
			}

			return false;
		}
		else
		{
			load(
					i_rexiSimVars,
					i_function_name,

					o_rexiCoefficients.alphas,
					o_rexiCoefficients.betas,
					o_rexiCoefficients.gamma,

					i_verbosity
				);

			return true;
		}
	}

#if 0
	static
	void testREXIphi0(
			std::vector<std::complex<T>> &alpha,
			std::vector<std::complex<T>> &beta_re,
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
			T correct_real = std::exp(std::complex<T>(0, 1)*x).real();

			std::complex<T> approx = 0;
			for (std::size_t n = 0; n < alpha.size(); n++)
				approx += (beta_re[n] / (std::complex<T>(0.0, x) + alpha[n]));
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
#endif

};


#endif /* SRC_INCLUDE_REXI_HPP_ */
