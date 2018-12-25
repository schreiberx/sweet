/*
 * test_rexi_ng.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <sweet/SimulationVariables.hpp>
#include <quadmath.h>
#include <rexi/REXI.hpp>
#include <rexi/REXI_File.hpp>
#include <rexi/REXIFunctions.hpp>


#define TEST_REXI_PDE_QUADPRECISION SWEET_QUADMATH

#if TEST_REXI_PDE_QUADPRECISION
	typedef __float128 T;
#else
	typedef double T;
#endif


REXIFunctions<T> rexiFunctions;



std::complex<double> fromTtoComplexDouble(
		const std::complex<T> &i_value
)
{
	std::complex<double> value;
	value.real(i_value.real());
	value.imag(i_value.imag());

	return value;
}


std::complex<T> fromComplexDoubleToT(
		const std::complex<double> &i_value
)
{
	std::complex<T> value;
	value.real(i_value.real());
	value.imag(i_value.imag());

	return value;
}



typedef std::complex<double> cplx;


void solveLalpha(
		const cplx i_lambda,	///< stiffness
		const double i_dt,		///< timestep size
		const cplx &i_alpha,	///< REXI shift
		const cplx &i_beta,		///< REXI shift
		const cplx i_u[2],		///< state variable
		cplx o_u[2]		///< output after REXI computation
)
{
	cplx i(0, 1);

	const cplx alpha = i_alpha/i_dt;
	const cplx beta = i_beta/i_dt;

	cplx val = cplx(1.0)/(i_lambda*i_lambda - alpha*alpha);
	cplx ia = i*i_lambda;

	o_u[0] = beta*(val*(-alpha*i_u[0] + ia*i_u[1]));
	o_u[1] = beta*(val*(-ia*i_u[0] - alpha*i_u[1]));
}



void computeLU(
		const cplx &i_lambda,
		const cplx i_u[2],		///< state variable
		cplx o_LU[2]		///< output after REXI computation
)
{
	o_LU[0] = (i_u[1]*i_lambda);
	o_LU[1] = (-i_u[0]*i_lambda);
}



void analyticalIntegration(
		const cplx &i_lambda,	///< stiffness
		double i_dt,		///< timestep size
		const cplx i_u[2],		///< state variable
		cplx o_u[2]		///< output after REXI computation
)
{
	cplx I(0, 1);

	cplx tmp[2];
	tmp[0] = cplx(0.5)*(-I*i_u[0] + i_u[1]);
	tmp[1] = cplx(0.5)*(I*i_u[0] + i_u[1]);

	cplx K = i_dt*i_lambda;
	tmp[0] = fromTtoComplexDouble(rexiFunctions.eval(fromComplexDoubleToT( K)))*tmp[0];
	tmp[1] = fromTtoComplexDouble(rexiFunctions.eval(fromComplexDoubleToT(-K)))*tmp[1];

	o_u[0] = I*tmp[0] - I*tmp[1];
	o_u[1] = tmp[0] + tmp[1];
}



void rexiIntegration(
		const cplx &i_lambda,	///< stiffness
		double i_dt,		///< timestep size
		std::vector<cplx> &i_alpha,
		std::vector<cplx> &i_beta,
		cplx io_u[2]		///< state variable
)
{
	cplx o_u[2] = {0.0, 0.0};

	for (std::size_t i = 0; i < i_alpha.size(); i++)
	{
		cplx ru[2];

		solveLalpha(
				i_lambda,	///< stiffness
				i_dt,
				i_alpha[i],	///< REXI shift
				i_beta[i],
				io_u,		///< state variable
				ru
		);

		o_u[0] += ru[0];
		o_u[1] += ru[1];
	}

	io_u[0] = o_u[0];
	io_u[1] = o_u[1];
}


void rexiIntegrationEValues(
		const cplx &i_lambda,	///< stiffness
		double i_dt,		///< timestep size
		std::vector<cplx> &i_alpha,
		std::vector<cplx> &i_beta,
		cplx io_u[2]		///< state variable
)
{
	cplx I(0, 1);

	cplx tmp[2];
	tmp[0] = cplx(0.5)*(-I*io_u[0] + io_u[1]);
	tmp[1] = cplx(0.5)*(I*io_u[0] + io_u[1]);

	cplx tmp_exp[2] = {0.0, 0.0};
	for (std::size_t i = 0; i < i_alpha.size(); i++)
	{
		const cplx alpha = i_alpha[i]/i_dt;
		const cplx beta = i_beta[i]/i_dt;

		tmp_exp[0] += beta/(i_lambda + alpha);
		tmp_exp[1] += beta/(-i_lambda + alpha);
	}

	tmp_exp[0] = tmp_exp[0]*tmp[0];
	tmp_exp[1] = tmp_exp[1]*tmp[1];

	io_u[0] = I*tmp_exp[0] - I*tmp_exp[1];
	io_u[1] = tmp_exp[0] + tmp_exp[1];
}



#if 0
T lenreal(
		const cplx &a,
		const cplx &b
)
{
	T av = a.real();
	T bv = b.real();
	return rexiFunctions.l_sqrt(av*av+bv*bv);
}
#endif



/**
 * \return \f$ Re(cos(x) + i*sin(x)) = cos(x) \f$
 */
std::complex<double> approx_returnComplex(
		std::vector< std::complex<double> > &alpha,
		std::vector< std::complex<double> > &beta_re,
		double i_x
)
{
	std::complex<double> sum = 0;

	std::size_t S = alpha.size();

	for (std::size_t n = 0; n < S; n++)
		sum += beta_re[n] / (std::complex<double>(0, i_x) + alpha[n]);

	return sum;
}



int main(
		int i_argc,
		char *const i_argv[]
)
{
	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr, false))
		return -1;

	if (simVars.rexi.rexi_method != "ci")
		FatalError("This test is for rexi_method=='ci' only");

	std::string function_names[9] =
	{
			"phi0",
			"phi1",
			"phi2",
			"phi3",
			"phi4",
			"phi5",

			"ups1",
			"ups2",
			"ups3",
	};

	double max_error_threshold = 1e-8;

	for (int fun_id = 0; fun_id < 9; fun_id++)
	{
		std::string &function_name = function_names[fun_id];

		rexiFunctions.setup(function_name);

		std::cout << "******************************************************" << std::endl;
		std::cout << function_name << " - Testing time stepping with PDE" << std::endl;
		std::cout << "******************************************************" << std::endl;

		/*
		 * Stiffness is specified here
		 *
		 * Oscillatory stiffness: imaginary-only
		 */
		cplx lambda = {0.0, 1.0};
		std::vector<cplx> alpha, beta;

		REXI<T>::load(
				&simVars.rexi,
				function_name,
				alpha,
				beta,
				1,
				simVars.misc.verbosity
			);

#if 0
		std::cout << std::endl;
		std::cout << "ALPHA:" << std::endl;
		for (std::size_t i = 0; i < alpha.size(); i++)
			std::cout << i << ": " << alpha[i] << std::endl;

		std::cout << std::endl;
		std::cout << "BETA:" << std::endl;
		for (std::size_t i = 0; i < beta.size(); i++)
			std::cout << i << ": " << beta[i] << std::endl;
#endif


		cplx U0[2];
		U0[0] = {1.0, 1.0};
		U0[1] = {2.0, 3.0};

		cplx U[2];
		U[0] = U0[0];
		U[1] = U0[1];

		analyticalIntegration(
				lambda,
				0,
				U0,
				U
			);

		double timestep_size = 0.3;
		double simtime = 0;

		cplx evalue = lambda*timestep_size;
		std::cout << "Eigenvalue 1: " << -fromTtoComplexDouble(evalue) << std::endl;
		std::cout << "Eigenvalue 2: " << fromTtoComplexDouble(evalue) << std::endl;

		int tnr_max = (int)(M_PI*2.0/(timestep_size*lambda.imag()));
		int tnr = 0;
		while (true)
		{
			cplx Uanal[2];	// analytical solution

			analyticalIntegration(
					lambda,
					simtime,
					U0,
					Uanal
				);

			double error = std::max(
					std::abs(	(double)Uanal[0].real() - (double)U[0].real()	),
					std::abs(	(double)Uanal[1].real() - (double)U[1].real()	)
				);

			std::cout << "t = " << (double)simtime;
			std::cout << "\tU=(" << fromTtoComplexDouble(U[0]) << ", " << fromTtoComplexDouble(U[1]) << ")";
			std::cout << "\tUanal=(" << fromTtoComplexDouble(Uanal[0]) << ", " << fromTtoComplexDouble(Uanal[1]) << ")";
//			std::cout << "\tUreallen=(" << (double)lenreal(U[0], U[1]) << ")";
//			std::cout << "\tUanalreallen=(" << (double)lenreal(Uanal[0], Uanal[1]) << ")";
			std::cout << "\tUmaxerr=(" << error << ")";
			std::cout << std::endl;

#if 1
			if (error > max_error_threshold)
				FatalError("Error too large");
#endif

			if (tnr >= tnr_max)
				break;

			if (function_name == "phi0")
			{
	#if 1
				rexiIntegration(
						lambda,
						timestep_size,
						alpha,
						beta,
						U
					);
	#else
				rexiIntegrationEValues(
						lambda,
						timestep_size,
						alpha,
						beta,
						U
					);
	#endif
			}
			else
			{
				/*
				 * Integrate over entire interval since there's a damping,
				 * hence a memorization of the state
				 */
				U[0] = U0[0];
				U[1] = U0[1];
				rexiIntegration(
						lambda,
						simtime+timestep_size,
						alpha,
						beta,
						U
					);
			}

#if 0
			if (simVars.rexi.use_half_poles)
			{
				// eliminate imaginary poles
				// this is required for solving a PDE
				U[0].imag(0);
				U[1].imag(0);
			}
#endif

			tnr++;
			simtime += timestep_size;
		}
	}


	return 0;
}
