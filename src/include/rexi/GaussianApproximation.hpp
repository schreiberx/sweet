/*
 * GaussianApproximation.hpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_

#include <sweetmath.hpp>
#include <iostream>
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

	GaussianApproximation(
			int i_L = 11	///< L
	)
	{
		L = i_L;

		if (L == 0)
			L = 11;

		/*
		 * mu and a coefficients from
		 * "A high-order time-parallel scheme for solving wave propagation problems via the direct construction of an approximate time-evolution operator", Haut et.al.
		 */
		if (L == 11)
		{
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
		else if (L == -11)
		{
			/**
			 * Recomputed coefficients with higher accuracy
			 */
			mu = {	-4.315321510875024024755930440733209252357, 0};

			L = 11;
			a.resize(L*2+1);
			a = {
					{-1.08463547512905881040e-7,	+2.9363558093024306997e-8},
					{5.77935685299374656e-9,		-9.1481151228521658961e-7},
					{3.7038358208237924102e-6,		+6.600044740446944982e-7},
					{-2.6814569755414357050e-6,		+0.0000113596973319111859617},
					{0.000014671580056339240006,	-0.000031400644103252651400},
					{-0.00107554419399313741859,	-0.00047327406930292107703},
					{0.003817126636360349944,		+0.017839292571877041899},
					{0.12124171418932352950,		-0.12326964235131572523},
					{-0.97749882474018490175,		-0.18771224971482819432},
					{1.3432857754271598782,			+3.2034709703859047814},
					{4.0724088186063580608,			-6.1237563578921685448},
					{-9.4426992314956841312,		+0.},
					{4.0724088069324508865,			+6.1237563189376063022},
					{1.3432858438021879621,			-3.2034710111648174724},
					{-0.97749875271985153802,		+0.18771233382364543862},
					{0.12124162750886421924,		+0.12326973063490963278},
					{0.003817044967316809891,		-0.017839369918773727991},
					{-0.00107548532437296570988,	+0.00047321522033559951237},
					{0.000014704798043942080216,	+0.000031438120209716725927},
					{-2.7009271345856269452e-6,		-0.0000113454977527733507622},
					{3.6996919184365413322e-6,		-6.678127118632243357e-7},
					{7.90833808793635723e-9,		+9.1419674399960741488e-7},
					{-1.08466768909193194543e-7,	-2.9079438997683162566e-8}
			};
		}
		else
		{
			std::cerr << "invalid value of L: " << i_L << std::endl;
			exit(1);
		}
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
