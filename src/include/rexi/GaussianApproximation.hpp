/*
 * GaussianApproximation.hpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_

#include <complex>
#include <vector>


/**
 * This class provides the weights and coefficients for the
 * approximation of a Gaussian
 *
 * \f$
 * 	  exp(-(x*x)/(4*h*h))/sqrt(4 \pi)
 * \f$
 *
 * with a sum over complex rational functions.
 *
 * See e.g. Near optimal rational approximations of large data sets, Damle et. al.
 */
class GaussianApproximation
{
	typedef std::complex<double> complex;

public:
	complex mu;				///< average
	std::vector<complex> a;	///< weights for approximation
	int L;					///< 2*L+1 = number of weights

	GaussianApproximation()
	{
		/*
		 * mu and a coefficients from
		 * "A high-order time-parallel scheme for solving wave propagation problems via the direct construction of an approximate time-evolution operator", Haut et.al.
		 */
		mu = {	-4.315321510875024, 0};

		L = 11;
		a.resize(L*2+1);
		a = {
				{	-1.0845749544592896e-7,		2.77075431662228e-8		},
				{	1.858753344202957e-8,		-9.105375434750162e-7	},
				{	3.6743713227243024e-6,		7.073284346322969e-7	},
				{	-2.7990058083347696e-6,		0.0000112564827639346	},
				{	0.000014918577548849352,	-0.0000316278486761932	},
				{	-0.0010751767283285608,		-0.00047282220513073084	},
				{	0.003816465653840016,		0.017839810396560574	},
				{	0.12124105653274578,		-0.12327042473830248	},
				{	-0.9774980792734348,		-0.1877130220537587		},
				{	1.3432866123333178,			3.2034715228495942		},
				{	4.072408546157305,			-6.123755543580666		},
				{	-9.442699917778205,			0.						},
				{	4.072408620272648,			6.123755841848161		},
				{	1.3432860877712938,			-3.2034712658530275		},
				{	-0.9774985292598916,		0.18771238018072134		},
				{	0.1212417070363373,			0.12326987628935386		},
				{	0.0038169724770333343,		-0.017839242222443888	},
				{	-0.0010756025812659208,		0.0004731874917343858	},
				{	0.000014713754789095218,	0.000031358475831136815	},
				{	-2.659323898804944e-6,		-0.000011341571201752273	},
				{	3.6970377676364553e-6,		-6.517457477594937e-7	},
				{	3.883933649142257e-9,		9.128496023863376e-7	},
				{	-1.0816457995911385e-7,		-2.954309729192276e-8	}
		};
	}



	/**
	 * directly evaluate basis function which is to be approximated
	 */
	double evalGaussian(
			double x,	///< x-coefficient for Gaussian basis function
			double h	///< h-coefficient for Gaussian basis function
	)
	{
		return std::exp(-(x*x)/(4*h*h))/std::sqrt(4.0*M_PIl);
	}


	/**
	 * evaluate approximation of Gaussian basis function
	 *
	 * with sum of complex rational functions
	 */
	double approxGaussian(
			double x,	///< x-coefficient for Gaussian basis function
			double h	///< h-coefficient for Gaussian basis function
	)
	{
		// scale x, since it depends linearly on h:
		// x^2 ~ h^2
		x /= h;

		double sum = 0;

		for (int l = 0; l < 2*L+1; l++)
		{
			int j = l-L;

			// WORKS with max error 7.15344e-13
			sum += (a[l]/(complex(0, x) + mu + complex(0, j))).real();
		}

		return sum;
	}

	void print()
	{
		std::cout << "mu: " << mu << std::endl;
		for (std::size_t i = 0; i < a.size(); i++)
		{
			std::cout << "a[" << i << "]: " << a[i] << std::endl;
		}
	}
};


#endif /* SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_ */
