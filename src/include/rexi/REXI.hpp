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
	bool load(
		EXP_SimulationVariables *i_rexiSimVars,
		const std::string &i_function_name,

		std::vector<std::complex<T>> &o_alpha,
		std::vector<std::complex<T>> &o_beta,
		std::complex<T> &o_gamma,

		int i_verbosity
	)
	{
		o_alpha.clear();
		o_beta.clear();

		if (i_rexiSimVars->exp_method == "terry")
		{
			std::cout << "WARNING: This way of using REXI is deprecated" << std::endl;
			/// REXI stuff
			REXI_Terry<T, T> rexi_terry;
			rexi_terry.setup(i_function_name, i_rexiSimVars->terry_h, i_rexiSimVars->terry_M, i_rexiSimVars->terry_L, i_rexiSimVars->terry_reduce_to_half, i_rexiSimVars->terry_normalization);

			o_alpha = rexi_terry.alpha;
			o_beta = rexi_terry.beta;
		}
		else if (i_rexiSimVars->exp_method == "ci")
		{
			std::cout << "WARNING: This way of using REXI is deprecated" << std::endl;
			std::cout << "WARNING: Compile SWEET with quad precision if using an order > 2." << std::endl;

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
		else if (i_rexiSimVars->exp_method == "direct")
		{
			// no REXI, but direct exponential time integration
		}
		else
		{
			if (i_rexiSimVars->exp_method == "")
				SWEETError("Please specify rexi method via --rexi-method=[str]");
			else
				SWEETError(std::string("REXI method '")+ i_rexiSimVars->exp_method+"' is not supported!");
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

		return true;
	}

public:
	static
	bool is_rexi_method_supported(
			std::string &i_rexi_method
	)
	{
		return (
				i_rexi_method == "file" ||
				i_rexi_method == "terry" ||
				i_rexi_method == "ci" ||
				0
			);
	}


public:
	static
	void get_available_rexi_methods(
			std::stringstream &i_ss
	)
	{
		i_ss << "        'file': File-based REXI" << std::endl;
		i_ss << "        'terry': T-REXI (DEPRECATED, use 'file' interface)" << std::endl;
		i_ss << "        'ci': Cauchy Contour-based REXI method  (DEPRECATED, use 'file' interface)" << std::endl;
	}


public:
	static
	bool load(
			EXP_SimulationVariables *i_rexiSimVars,
			const std::string &i_function_name,

			REXICoefficients<T> &o_rexiCoefficients,

			int i_verbosity
	)
	{
		if (i_rexiSimVars->exp_method == "file")
		{
			for (
				std::vector<EXP_SimulationVariables::REXIFile>::iterator iter = i_rexiSimVars->p_rexi_files_processed.begin();
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

};


#endif /* SRC_INCLUDE_REXI_HPP_ */
