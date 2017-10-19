/*
 * PhiApproximation.hpp
 *
 *  Created on: 10 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_INCLUDE_REXI_REXI_TERRY_FUNAPPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_REXI_TERRY_FUNAPPROXIMATION_HPP_

#include <sweet/sweetmath.hpp>
#include <libmath/DQStuff.hpp>
#include <libmath/GaussQuadrature.hpp>
#include <rexi/REXI_Terry_GaussianApproximation.hpp>



/**
 * This class computes an approximation of functions suitable for the
 * computation of exponential integrator formulations including non-linearities
 *
 * \f$
 * 	phi1(x) := (e^{ix}-1)/ix
 * \f$
 *
 * with a combination of Gaussians.
 *
 * For more information, see Terry Haut et al. paper
 * "A high-order time-parallel scheme for solving wave propagation problems via
 * the direct construction of an approximate time-evolution operator"
 *
 * We use a numerical quadrature to compute the following coefficients (see Terry's paper):
 *
 * \f$
 * 		c_m := \int_{-1/(2h)}^{1/2h} h * \exp(-2 \pi i m h \xi) \tilde\phi_N(\xi) / \tilde\psi_h(\xi) d \xi
 * \f$
 * with the Gaussian
 * \f$
 * 		\psi_h(\xi) = 1.0/\sqrt(4.0*\pi) \exp(-x*x/(4*h*h))
 * \f$
 * which is in Fourier space
 * [see https://www.wolframalpha.com/input/?i=Fourier+transform+sqrt(2*pi)*(1%2Fsqrt(4*pi)+*+exp(-x*x%2F(4*h*h))))) ]
 * \f$
 * 		\tilde\psi_h(\xi) = h*\exp(-h*h*\xi*\xi)
 * \f$
 *
 * Combining both equations leads to
 *
 * \f$
 * 		c_m := \int_{-1/(2h)}^{1/2h}  exp(-2 \pi i m h \xi) / (h* \exp(-h*h*\xi*\xi)) * phiN(\xi) d \xi
 * \f$
 */
template <
#if SWEET_QUADMATH
	typename T = __float128
#else
	typename T = double
#endif
>
class REXI_Terry_FunApproximation
{
	typedef std::complex<T> TComplex;

	REXI_Terry_GaussianApproximation<T> ga;


	T h;
	int M;

	const T pi = DQStuff::fromString<T>("3.14159265358979323846264338327950288");
	const T pi2 = pi*DQStuff::fromString<T>("2.0");
	const T pi4 = pi*DQStuff::fromString<T>("4.0");
	const T sqrtpi4 = DQStuff::sqrt(pi4);
	const std::complex<T> I = std::complex<T>(0, 1);
	T int_threshold;

public:
	std::vector<TComplex> b;



	T heaviside(T i_x)
	{
		/*
		 * Maple defines 0 to be undefined
		 * https://www.maplesoft.com/support/help/maple/view.aspx?path=Heaviside
		 */
		if (i_x < 0)
			return 0;

		if (i_x > 0)
			return 1;

		FatalError("Heaviside not defined at 0");
		return 0;
	}



	T signum(T i_x)
	{
		/*
		 * Maple defines 0 to be undefined
		 * https://www.maplesoft.com/support/help/maple/view.aspx?path=Heaviside
		 */
		if (i_x < 0)
			return -1;

		if (i_x > 0)
			return 1;

		return 0;
	}



	/**
	 * Fourier transformed approximation rational functions in spectral space
	 *
	 * This sum represents the Gaussian basis function, but includes the imaginary components
	 */
	std::complex<T> gaussian_rat_approx_spec(
			T i_xi
	)
	{
		std::complex<T> sum = {0,0};

		for (int l = -ga.L; l <= ga.L; l++)
		{
			std::complex<T> term = -pi*h*DQStuff::exp(
					std::complex<T>(pi2*h*(I*(T)l + ga.mu)*i_xi)
				)*(
					(T)2.0*heaviside(i_xi)-(T)1-signum(ga.mu)
				);

			term *= ga.a[l+ga.L];
			term *= 0.5;		/// TODO: For some reason, there's a factor 0.5 required

			sum += term;
		}

		return sum;
	}

	std::complex<T> rexi_int_fun(
			int m,
			T i_xi,
			std::function<std::complex<T>(T)> fun_phi_tilde	///< spectral representation of function to be approximated,
	)
	{
		/*
		 * Note that we had to change the equation from
		 *
		 * DQStuff::expIm(-pi2*(T)m*h*i_xi)/(h*DQStuff::exp(-h*h*i_xi*i_xi)) * fun_phi_tilde(i_xi);
		 *
		 * to
		 *
		 * DQStuff::expIm(pi2*(T)m*h*i_xi)/(h*DQStuff::exp(-h*h*i_xi*i_xi*pi2*pi2)) * fun_phi_tilde(i_xi*pi2) * h;
		 *
		 * to get correct coefficients. This is since FT's are differently defined.
		 */
#if 1
		/*
		 * Use analytical representation of Gaussian function in spectral space
		 */
		T psi_h = h*DQStuff::exp(-h*h*i_xi*i_xi*pi2*pi2);
		return h * DQStuff::expIm(pi2*(T)m*h*i_xi)/psi_h * fun_phi_tilde(i_xi*pi2);

#elif 1

		/*
		 * Use rational approximation (of Gaussian function) in spectral space
		 */
		std::complex<T> psi_h = gaussian_rat_approx_spec(1.0/pi2);

		return h * DQStuff::expIm(pi2*(T)m*h*i_xi)/psi_h * fun_phi_tilde(i_xi*pi2);

#endif

	};

	std::complex<T> computeREXICoefficient(
			int m,
			std::function<std::complex<T>(T)> fun_phi_tilde,	///< spectral representation of function to be approximated,
			T i_int_start,
			T i_int_end
	)
	{
		T real = GaussQuadrature::integrate5_intervals_adaptive_recursive<T>(
						i_int_start,	// start of quadrature
						i_int_end,	// end of quadrature
						[&](T xi) -> T
						{
							return rexi_int_fun(m, xi, fun_phi_tilde).real();
						},
						int_threshold
					);

		T imag = GaussQuadrature::integrate5_intervals_adaptive_recursive<T>(
						(T)i_int_start,	// start of quadrature
						(T)i_int_end,	// end of quadrature
						[&](T xi) -> T
						{
							return rexi_int_fun(m, xi, fun_phi_tilde).imag();
						},
						int_threshold
					);

		return std::complex<T>(real, imag);
	}



	REXI_Terry_FunApproximation()
	{
	}



	REXI_Terry_FunApproximation(
			const std::string &i_function_name,
			T i_h,
			int i_M
	)
	{
		setup(i_function_name, i_h, i_M);
	}



	void setup(
			const std::string &i_function_name,
			T i_h,
			int i_M
	)
	{
		if (typeid(T) == typeid(double))
		{
			int_threshold = 1e-12;
		}
#if SWEET_QUADMATH
		else if (typeid(T) == typeid(__float128))
		{
			int_threshold = 1e-14;
		}
#endif
		else
		{
			FatalError("Type not supported!");
		}

		h = i_h;
		M = (i_M == -1 ? (T)i_M/h : i_M);
		b.resize(i_M*2+1);

#if SWEET_QUADMATH
#if 1
		__float128 asdf = DQStuff::fromString<T>("3.14159265358979323846264338327950288");
		if (asdf - pi != 0)
			FatalError("Compiled constant not equal to string-induced constant!");
#endif
#endif

		if (i_function_name == "phi0")
		{
			/*
			 * This is only for testing purpose, since phi0 can be directly computed
			 *
			 * Since it's a dirac function, we have to compute the quadrature over a tiny interval at [1-eps;1+eps]
			 */

			auto fun = [&](T i_xi)	-> std::complex<T>
			{
				return pi2;
			};
			for (int m = -i_M; m <= i_M; m++)
			{
				b[m+M] = rexi_int_fun(m, 1.0/pi2, fun);
			}
		}
		else if (i_function_name == "phi1")
		{
			/*
			 * \f$
			 * 		\phi_1(\xi) = 2\pi for -1/(2*\pi) <= \xi <= 0
			 * \f$
			 *
			 * https://www.wolframalpha.com/input/?i=Fourier+transform+(exp(i*x%2F(2*pi))-1)%2F(i*x)*sqrt(2*pi)
			 */
			auto fun = [&](T i_xi)	-> std::complex<T>
			{
				return pi2;
			};

			for (int m = -i_M; m <= i_M; m++)
			{
//				std::cout << "Computing REXI coefficient " << m << std::endl;
				b[m+M] = computeREXICoefficient(m, fun, -1.0/pi2, 0);
//				b[m+M] = computeREXICoefficient(m, fun, -1.0/(2.0*h), 0);
//				b[m+M] = computeREXICoefficient(m, fun, -1.0/pi2, -1e-12);

//				std::cout << (double)b[m+M].real() << ", " << (double)b[m+M].imag() << std::endl;
			}
		}
		else if (i_function_name == "phi2")
		{
			/*
			 * \f$
			 * 		\phi_2(\xi) = 2\pi for -1/(2*\pi) <= \xi <= 0
			 * \f$
			 *
			 * https://www.wolframalpha.com/input/?i=Fourier+transform+(exp(i*x%2F(2*pi))-1)%2F(i*x)*sqrt(2*pi)
			 */
			auto fun = [&](T i_xi)	-> std::complex<T>
			{
				return pi2;
			};

			for (int m = -i_M; m <= i_M; m++)
				b[m+M] = computeREXICoefficient(m, fun, -1.0/pi2, 0);
//				b[m+M] = computeREXICoefficient(m, fun, -1.0/pi2, -1e-12);
		}
		else
		{
			FatalError(std::string("Function ")+i_function_name+std::string(" is not supported"));
		}
	}



	void print()
	{
		for (std::size_t i = 0; i < b.size(); i++)
		{
			std::cout << "b[" << i << "]: " << b[i] << std::endl;
		}
	}


	TComplex approx(
			T i_x
	)
	{
		/// \f$ \sum_{m=-M}^{M}{b_m \psi_h(x+m*h)} \f$

		TComplex sum = 0;
		for (int m = -M; m < M+1; m++)
		{
			sum += b[m+M] * ga.approxGaussian(i_x+((T)m)*h, h);
		}
		return sum;
	}
};



#endif
