/*
 * test_rexi.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <rexi/REXI_Terry.hpp>
#include <rexi/REXI_Terry_ExponentialApproximation.hpp>
#include <rexi/REXI_Terry_GaussianApproximation.hpp>
#include <iostream>
#include <rexi/REXIFunctions.hpp>
#include <sweet/SimulationVariables.hpp>
#include "../include/sweet/plane/PlaneDataComplex.hpp"


#if 0
	typedef __float128 TGeneration;
#else
	typedef double TGeneration;
#endif
typedef std::complex<TGeneration> TComplexGeneration;



typedef double T;
typedef std::complex<T> TComplex;



int main(
		int i_argc,
		char *const i_argv[]
)
{
	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr, false))
	{
		return -1;
	}

	double max_error_threshold = 1e-9;


#if 1
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING GAUSSIAN APPROXIMATION (real)" << std::endl;
		std::cout << "******************************************************" << std::endl;
		REXI_Terry_GaussianApproximation<TGeneration> ga(simVars.rexi.L);

		for (double h = 0.2; h > 0.001; h *= 0.5)
		{
			double start = -100.0;
			double end = 100.0;
			double step_size = 0.001;

			double max_error_real = 0.0;
			double max_error_imag = 0.0;

#if 1
			start = -5.0;
			end = 5.0;
			step_size = 0.05;
#endif

			for (double x = start; x < end; x += step_size)
			{
				TGeneration pi = DQStuff::fromString<TGeneration>("3.14159265358979323846264338327950288");
				TGeneration pi4 = pi*DQStuff::fromString<TGeneration>("4.0");
				TGeneration sqrtpi4 = DQStuff::sqrt(pi4);

				std::complex<double> approx = ga.approxGaussian_Complex(x, h);
				std::complex<double> analyt = ga.evalGaussian(x, h);

				double t = 2*x/std::pow(x*x+1.0, 2.0);
				analyt.imag(t);

				std::complex<double> error = analyt - approx;

				max_error_real = DQStuff::max(max_error_real, DQStuff::abs(error.real()));
				max_error_imag = DQStuff::max(max_error_imag, DQStuff::abs(error.imag()));

				std::cout << x << "\t";
				std::cout << approx.real() << "\t" << approx.imag() << "\t";
				std::cout << analyt.real() << "\t" << analyt.imag();
				std::cout << std::endl;
			}
//			exit(-1);

			std::cout << "max_error_real: " << max_error_real << " for h " << h << std::endl;
//				std::cout << "max_error_imag: " << max_error_imag << " for h " << h << std::endl;

			if (DQStuff::abs(max_error_real) > max_error_threshold)
			{
				std::cerr << "MAX ERROR THRESHOLD EXCEEDED!" << std::endl;
				exit(-1);
			}
		}
	}
#endif

// TODO: REACTIVATE!!!!!!!!!!!!!!
// TODO: REACTIVATE!!!!!!!!!!!!!!
// TODO: REACTIVATE!!!!!!!!!!!!!!
#if 1
	if (1)
	{
		std::cout << "******************************************************" << std::endl;
		std::cout << "EVALUATING EXPONENTIAL (e^(ix)) APPROXIMATION with approx Gaussian" << std::endl;
		std::cout << "******************************************************" << std::endl;

		for (T h = 0.2; h > 0.01; h *= 0.5)
		{
			int M = 32/h;
			REXI_Terry_ExponentialApproximation<T> ea(h, M);

			T start = -M*h*0.95;
			T end = -start;
			T step_size = 0.1;

			T max_error = 0;

			for (T x = start; x < end; x += step_size)
			{
				std::complex<T> diff = ea.eval(x) - ea.approx(x);
				T error = DQStuff::max(DQStuff::abs(diff.real()), DQStuff::abs(diff.imag()));
				max_error = DQStuff::max(max_error, error);
			}

			std::cout << "max_error: " << (double)max_error << " for h " << h << " and M " << M << std::endl;

			if (DQStuff::abs(max_error) > max_error_threshold)
			{
				std::cerr << "MAX ERROR THRESHOLD EXCEEDED!" << std::endl;
				exit(-1);
			}
		}
	}
#endif



	for (int fun_id = 1; fun_id <= 1; fun_id++)
	{
		std::string function_name;
		switch(fun_id)
		{
		case 0:
			function_name = "phi0";
			break;

		case 1:
			function_name = "phi1";
			break;

		case 2:
			function_name = "phi2";
			break;

		default:
			FatalError("this phi function is not implemented");
		}


		if (simVars.rexi.use_half_poles)
		{
			std::cout << "Skipping testing for valid real parts since halving of poles is used" << std::endl;
		}
		else
		{
			std::cout << "******************************************************" << std::endl;
			std::cout << "PHI " << fun_id << " - REXI real: Test for partition of unity and accuracy" << std::endl;
			std::cout << "******************************************************" << std::endl;

			for (T h = 0.2; h >= 0.05; h *= 0.5)
			{
				int Mstart = 128;
				if (fun_id > 0)
					Mstart = 1024;

				for (int M = Mstart; M <= 1024; M *= 2)
				{
					std::cout << "******************************************************" << std::endl;
					std::cout << "M: " << M << std::endl;
					std::cout << "h: " << (double)h << std::endl;
					std::cout << "******************************************************" << std::endl;
					std::cout << "Setup coefficients... " << std::flush;
					REXI_Terry<TGeneration, T> rexi(function_name, h, M, simVars.rexi.L, simVars.rexi.use_half_poles, simVars.rexi.normalization);
					std::cout << "OK" << std::endl;


					REXIFunctions<__float128> rexiFunctions(function_name);

					// REXI approximates the interval [-M*h;M*h] but gets inaccurate close to the interval boundaries
					T start = -M*h*0.9;
					if (fun_id > 0)
						start = start*0.5;

					T end = -start;
					T step_size = 0.011;

					T max_error_real = 0.0;
					T max_error_imag = 0.0;

					double max_error_threshold_imag = max_error_threshold;

					/*
					 * TODO: Analyze this required increase error threshold for the imaginary axis
					 */
//					if (fun_id > 0)
//						max_error_threshold_imag *= 1e+5;

					for (T x = start; x < end; x += step_size)
					{
						std::complex<__float128> correct_ = rexiFunctions.eval(std::complex<__float128>(0, x));
						TComplex correct(correct_.real(), correct_.imag());
						TComplex approx = rexi.approx_returnComplex(x);

//						std::cout << (double)x << ": " << correct << "\t" << approx << "\t" << (correct-approx) << std::endl;
//						continue;

						if (DQStuff::abs(approx.real()) > 1.0)
							std::cerr << "approx value_real " << approx.real() << " not bounded by unity (just a warning and not a problem) at x=" << (double)x << std::endl;

//						if (DQStuff::abs(approx.imag()) > 1.0)
//							std::cerr << "approx value_imag " << approx.imag() << " not bounded by unity (just a warning and not a problem) at x=" << (double)x << std::endl;

						T error_real = DQStuff::abs(correct.real() - approx.real());
						T error_imag = DQStuff::abs(correct.imag() - approx.imag());

						max_error_real = DQStuff::max(max_error_real, error_real);
						max_error_imag = DQStuff::max(max_error_imag, error_imag);
					}
//					exit(1);

					std::cout << "max_error_real: " << (double)max_error_real << " for h " << (double)h << " and M " << M << std::endl;

					if (max_error_real > max_error_threshold)
					{
						std::cout << "max_error_threshold: " << max_error_threshold << std::endl;
						FatalError("MAX ERROR THRESHOLD EXCEEDED for real values!");
					}

					std::cout << "max_error_imag: " << (double)max_error_imag << " for h " << (double)h << " and M " << M << std::endl;

#if 0
					if (max_error_imag > max_error_threshold_imag)
					{
						std::cout << "max_error_threshold_imag: " << max_error_threshold_imag << std::endl;
						FatalError("MAX ERROR THRESHOLD EXCEEDED for imaginary value!");
					}
#endif
				}
			}
		}
	}


	return 0;
}
