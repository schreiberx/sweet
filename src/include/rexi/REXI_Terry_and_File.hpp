/*
 * REXI_OLD_and_NG.hpp
 *
 *  Created on: 4 Aug 2017
 *      Author: martin
 */

#ifndef SRC_INCLUDE_REXI_REXI_TERRY_AND_FILE_HPP_
#define SRC_INCLUDE_REXI_REXI_TERRY_AND_FILE_HPP_


#include <rexi/RexiFile.hpp>
#include <rexi/REXI_SimulationVariables.hpp>
#include <rexi/REXI.hpp>
#include <vector>
#include <complex>


/**
 * Interface to load either REXI via Terrys method or via File
 */
class REXI_Terry_or_File
{
public:
	static
	void load(
			REXI_SimulationVariables *i_rexiSimVars,
			const std::string &i_function_name,

			std::vector<std::complex<double>> &alpha,
			std::vector<std::complex<double>> &beta_re,

			int i_verbosity = 0
	)
	{
		if (i_rexiSimVars->use_next_generation)
		{
			/// REXI next generation stuff
			RexiFile<double> rexiNG;

			bool retval = rexiNG.auto_load(
					i_function_name,
					i_rexiSimVars->ng_N,	/// N
					rexiNG.None(),			/// max_error
					i_rexiSimVars->ng_max_error_double_precision,			/// max_error_double_precision
					i_rexiSimVars->ng_test_min,
					i_rexiSimVars->ng_test_max,
					rexiNG.None(),			/// basis_function_scaling
					i_rexiSimVars->ng_h, 	/// basis_function_spacing
					rexiNG.None(),			/// basis_function rat shift

					i_rexiSimVars->use_half_poles,
					i_rexiSimVars->ng_faf_dir
				);

			if (!retval)
				FatalError(std::string("Not able to find coefficients for given constraints for function "+i_function_name));

			if (i_verbosity > 0)
				std::cout << "Loaded REXI coefficients from file '" << rexiNG.fafcoeffs.filename << "'" << std::endl;

			if (i_verbosity > 3)
			{
				rexiNG.fafcoeffs.output();
				rexiNG.fafcoeffs.outputWeights();
				rexiNG.output();
			}

			alpha = rexiNG.alpha;
			beta_re = rexiNG.beta_re;

#if 0
			std::size_t N = alpha.size();
			alpha.resize(N*2);
			beta_re.resize(N*2);

			for (std::size_t i = 0; i < N; i++)
			{
				alpha[i+N] = alpha[i] - std::complex<double>(0, M_PI/2.0);
				beta_re[i+N] = beta_re[i];
//				alpha[i] = -std::conj(alpha[i]);
//				beta_re[i] = -std::conj(beta_re[i]);
			}
#endif
		}
		else
		{
			/// REXI stuff
			REXI<double, double> rexi;

			rexi.setup(i_function_name, i_rexiSimVars->h, i_rexiSimVars->M, i_rexiSimVars->L, i_rexiSimVars->use_half_poles, i_rexiSimVars->normalization);

			alpha = rexi.alpha;
			beta_re = rexi.beta_re;

			if (i_verbosity > 3)
				rexi.output();
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

//			std::cout << "exp(I*" << x << ") ~~ " << approx_real << "\t" << correct_real << std::endl;
		}
	}

};


#endif /* SRC_INCLUDE_REXI_REXI_TERRY_AND_FILE_HPP_ */
