/*
 * Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_SPHERE_HELPERS_SPH_IDENTITIES_HPP_
#define SRC_SPHERE_HELPERS_SPH_IDENTITIES_HPP_

#include <cassert>

class SphereHelpers_SPHIdentities
{
public:
	typedef std::complex<double> Tcomplex;
	typedef double Treal;


	/*
	 * Normalization factor
	 */
	static inline
	double implicit_eps(double n, double m)
	{
		assert(n > 0);
		/*
		 * If you read this, check for n=0 and set the result to 0 since the nominator equals zero for n=0
		 */
		return std::sqrt(Treal(n*n-m*m)/Treal(4.0*n*n-1.0));
	}


	static inline
	Tcomplex implicit_J_scalar(Treal nr, Treal mr, Treal i_dt_two_omega)
	{
		assert(nr >= 0);

		if (nr == 0)
			return 1.0;		// Set integration constant (the imaginary part of equation below) to 0, hence return 1.0

		return Tcomplex(1.0, -i_dt_two_omega*mr/(nr*(nr+1.0)));
	}

	static inline
	Tcomplex implicit_J_scalar(Treal nr, Treal mr, Tcomplex i_dt_two_omega)
	{
		assert(nr >= 0);

		if (nr == 0)
			return 1.0;		// Set integration constant (the imaginary part of equation below) to 0, hence return 1.0

		return Tcomplex(1.0, 0) - Tcomplex(0, 1)*i_dt_two_omega*mr/(nr*(nr+1.0));
	}

	static inline
	Treal implicit_f_plus(Treal nr, Treal mr)
	{
		return nr/(nr+1.0) * implicit_eps(nr+1, mr);
	}


	static inline
	Treal implicit_f_minus(Treal nr, Treal mr)
	{
		return (nr+1.0)/nr * implicit_eps(nr, mr);
	}



public:

	inline
	static double D(double k, double m)
	{
		double n=k+1;
		assert(n >= 0);
//		if (n < 0)
//			n = -n-1;
		return ((2.0*n+1.0)*std::sqrt((n*n-m*m)/(4.0*n*n-1.0)));
	}

	inline
	static double E(double n, double m)
	{
		assert(n >= 0);
		return -n;
	}


	inline
	static double G(double n, double m)
	{
		n = n+1;
		return D(n-1,m) + E(n,m)*R(n-1,m);
	}


	inline
	static double H(double n, double m)
	{
		n = n-1;
		return E(n,m)*S(n+1,m);
	}

	inline
	static double R(double k, double m)
	{
		double n=k+1;

		if (n < 0)
			n = -n-1;	// if you remove me, the results will be really sad!

		if (std::abs(m) > n)
			return 0;

		assert (n*n-m*m >= 0);

		assert(n >= 0);
		return std::sqrt((n*n-m*m)/(4.0*n*n-1.0));
	}


	inline
	static double R_real(double k, double m)
	{
		double n=k+1;

		assert(std::abs(m) <= n);
		assert(n*n-m*m >= 0);
		assert(n >= 0);

		return std::sqrt((n*n-m*m)/(4.0*n*n-1.0));
	}


	inline
	static double S(double k, double m)
	{
		double n=k-1;

		if (n < 0)
			n = -n-1;

		if (std::abs(m) > n)
			return 0;

		assert (n*n-m*m >= 0);

		assert(n >= 0);
		return std::sqrt(((n+1.0)*(n+1.0)-m*m)/((2.0*n+1.0)*(2.0*n+3.0)));
	}


	inline
	static double S_real(double k, double m)
	{
		double n=k-1;

		assert(n >= 0);
		assert(std::abs(m) <= n);
		assert (n*n-m*m >= 0);

		return std::sqrt(((n+1.0)*(n+1.0)-m*m)/((2.0*n+1.0)*(2.0*n+3.0)));
	}

	inline
	static double A(double k, double m)
	{
		double n = k+2.0;
		return R(n-1,m)*R(n-2,m);
	}

	inline
	static double B(double n, double m)
	{
		return R(n-1,m)*S(n,m) + S(n+1,m)*R(n,m);
	}

	inline
	static double C(double k, double m)
	{
		double n = k-2.0;
		return S(n+1,m)*S(n+2,m);
	}
};


#endif
