/*
 * test_rexi_ng.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <sweet/SimulationVariables.hpp>
#include <quadmath.h>
#include <rexi/RexiFile.hpp>
#include <rexi/REXI_Terry_and_File.hpp>


#define TEST_REXI_PDE_QUADPRECISION 1

#if TEST_REXI_PDE_QUADPRECISION

	typedef __float128 T;
	typedef std::complex<T> cplx;

	cplx l_expcplx(const cplx &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = cexpq(z);

		return std::complex<double>(crealq(val), cimagq(val));
	}


	inline
	T l_sqrt(const T &i_value)
	{
		return sqrtq(i_value);
	}

	const cplx l_sqrtcplx(const cplx &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = csqrtq(z);

		return std::complex<double>(crealq(val), cimagq(val));
	}

	inline
	T eps_phi()
	{
		return 1e-10;
	}

	inline
	T eps_ups()
	{
		return 1e-10;
	}

	inline
	T pi2()
	{
		char *sp;
		T retval = (T)2.0*strtoflt128("3.1415926535897932384626433832795029", &sp);
		return retval;
	}

#else

	/*
	 * Double precision
	 * Might suffer of numerical double precision limited effects
	 */
	typedef double T;
	typedef std::complex<T> cplx;

	cplx l_expcplx(const cplx &i_value)
	{
		return std::exp(i_value);
	};

	T l_sqrt(T &i_value)
	{
		return std::sqrt(i_value);
	};

	cplx l_sqrtcplx(cplx &i_value)
	{
		return std::exp(i_value);
	};


	T eps_phi()
	{
		return 1e-10;
	}

	T eps_ups()
	{
		return 1e-10;
	}

	T pi2()
	{
		return (T)2.0*(T)M_PI;
	}

#endif



std::complex<double> toComplexDouble(
		const std::complex<__float128> &i_value
)
{
	std::complex<double> value;
	value.real(i_value.real());
	value.imag(i_value.imag());

	return value;
}







void solveLalpha(
		const cplx i_lambda,	///< stiffness
		const T i_dt,		///< timestep size
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

//	std::cout << toComplexDouble(i_lambda*i_lambda - alpha*alpha) << std::endl;

#if 1

	o_u[0] = beta*(val*(-alpha*i_u[0] + ia*i_u[1]));
	o_u[1] = beta*(val*(-ia*i_u[0] - alpha*i_u[1]));

#else

	cplx utmp[2] = {i_u[0], i_u[1]};
	utmp[0] *= beta;
	utmp[1] *= beta;

	o_u[0] = val*(-alpha*utmp[0] + ia*utmp[1]);
	o_u[1] = val*(-ia*utmp[0] - alpha*utmp[1]);

#endif

//	std::cout << toComplexDouble(o_u[0]) << "\t" << toComplexDouble(o_u[1]) << std::endl;
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
		T i_dt,		///< timestep size
		const cplx i_u[2],		///< state variable
		cplx o_u[2]		///< output after REXI computation
)
{
	cplx I(0, 1);

	cplx tmp[2];
	tmp[0] = cplx(0.5)*(-I*i_u[0] + i_u[1]);
	tmp[1] = cplx(0.5)*(I*i_u[0] + i_u[1]);

	tmp[0] = l_expcplx(i_dt*i_lambda)*tmp[0];
	tmp[1] = l_expcplx(-i_dt*i_lambda)*tmp[1];

	o_u[0] = I*tmp[0] - I*tmp[1];
	o_u[1] = tmp[0] + tmp[1];
}



void rexiIntegration(
		const cplx &i_lambda,	///< stiffness
		T i_dt,		///< timestep size
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
		T i_dt,		///< timestep size
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



T lenreal(
		const cplx &a,
		const cplx &b
)
{
	T av = a.real();
	T bv = b.real();
	return l_sqrt(av*av+bv*bv);
}



int main(
		int i_argc,
		char *const i_argv[]
)
{
	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr, false))
		return -1;

//	T max_error_threshold = 1e-9;

	for (int fun_id = 0; fun_id < 1; fun_id++)
	{
		std::ostringstream os;
		os << "phi" << fun_id;
		std::string function_name = os.str();

		std::cout << "******************************************************" << std::endl;
		std::cout << function_name << " - Testing time stepping" << std::endl;
		std::cout << "******************************************************" << std::endl;

		/*
		 * Stiffness is specified here
		 *
		 * Oscillatory stiffness: imaginary-only
		 */
		cplx lambda = {0.0, 1.0};

		std::vector<std::complex<double>> alpha_tmp, beta_tmp;

		REXI_Terry_or_File::load(
				&simVars.rexi,
				function_name,
				alpha_tmp,
				beta_tmp,
				simVars.misc.verbosity
			);

		if (!simVars.rexi.use_half_poles && function_name == "phi0")
		{
			REXI_Terry_or_File::testREXIphi0(
					alpha_tmp,
					beta_tmp,
					2.0
			);
		}

		std::vector<cplx> alpha, beta;
		alpha.resize(alpha_tmp.size());
		beta.resize(beta_tmp.size());

		for (std::size_t i = 0; i < alpha.size(); i++)
		{
			alpha[i].real(alpha_tmp[i].real());
			alpha[i].imag(alpha_tmp[i].imag());

			beta[i].real(beta_tmp[i].real());
			beta[i].imag(beta_tmp[i].imag());
		}


		cplx U0[2] = {1.0, 0.0};

		cplx U[2] = {U0[0], U0[1]};

		T timestep_size = 0.1;
		T simtime = 0;

		cplx evalue = lambda*timestep_size;
		std::cout << "Eigenvalue 1: " << -toComplexDouble(evalue) << std::endl;
		std::cout << "Eigenvalue 2: " << toComplexDouble(evalue) << std::endl;

		int tnr_max = 10;
		int tnr = 0;
		while (true)
		{
			cplx Ubutt[2];	// analytical solution

			analyticalIntegration(
					lambda,
					simtime,
					U0,
					Ubutt
				);

			std::cout << "t = " << (double)simtime;
			std::cout << "\tU=(" << toComplexDouble(U[0]) << ", " << toComplexDouble(U[1]) << ")";
			std::cout << "\tUbutt=(" << toComplexDouble(Ubutt[0]) << ", " << toComplexDouble(Ubutt[1]) << ")";
			std::cout << "\tUreallen=(" << (double)lenreal(U[0], U[1]) << ")";
			std::cout << "\tUbuttreallen=(" << (double)lenreal(Ubutt[0], Ubutt[1]) << ")";
			std::cout << "\tUmaxerr=(" << std::max(
									std::abs(	(double)Ubutt[0].real() - (double)U[0].real()	),
									std::abs(	(double)Ubutt[1].real() - (double)U[1].real()	)
								) << ")";
			std::cout << std::endl;

			if (tnr >= tnr_max)
				break;

#if 0
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
			if (simVars.rexi.use_half_poles)
			{
				// eliminate imaginary poles
				U[0].imag(0);
				U[1].imag(0);
			}

			tnr++;
			simtime += timestep_size;

		}
	}


	return 0;
}
