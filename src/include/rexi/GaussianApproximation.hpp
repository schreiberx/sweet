/*
 * GaussianApproximation.hpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_
#define SRC_INCLUDE_REXI_GAUSSIANAPPROXIMATION_HPP_

#include <sweet/sweetmath.hpp>
#include <libmath/DQStuff.hpp>
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
template <
	typename Tevaluation,	///< evaluation accuracy of coefficients
	typename TStorage		///< storage precision of coefficients - use quad precision per default
>
class GaussianApproximation
{
	typedef std::complex<Tevaluation> complex;
	typedef std::complex<TStorage> complexStorage;

	const TStorage pi = DQStuff::fromString<TStorage>("3.14159265358979323846264338327950288");
	const TStorage pi2 = pi*DQStuff::fromString<TStorage>("2.0");
	const TStorage pi4 = pi*DQStuff::fromString<TStorage>("4.0");
	const TStorage sqrtpi4 = DQStuff::sqrt(pi4);

public:
	TStorage mu;					///< average
	std::vector<complexStorage> a;	///< weights for approximation
	int L;							///< 2*L+1 = number of weights

	GaussianApproximation(
			int i_L = 0	///< L
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
			L = 11;
			a.resize(L*2+1);

#if 0
			mu = DQStuff::fromString<TStorage>("-4.315321510875024");

			a[0].real(DQStuff::fromString<TStorage>("-1.0845749544592896e-7"));		a[0].imag(DQStuff::fromString<TStorage>("2.77075431662228e-8"));
			a[1].real(DQStuff::fromString<TStorage>("1.858753344202957e-8"));		a[1].imag(DQStuff::fromString<TStorage>("-9.105375434750162e-7"));
			a[2].real(DQStuff::fromString<TStorage>("3.6743713227243024e-6"));		a[2].imag(DQStuff::fromString<TStorage>("7.073284346322969e-7"));
			a[3].real(DQStuff::fromString<TStorage>("-2.7990058083347696e-6"));		a[3].imag(DQStuff::fromString<TStorage>("0.0000112564827639346"));
			a[4].real(DQStuff::fromString<TStorage>("0.000014918577548849352"));	a[4].imag(DQStuff::fromString<TStorage>("-0.0000316278486761932"));
			a[5].real(DQStuff::fromString<TStorage>("-0.0010751767283285608"));		a[5].imag(DQStuff::fromString<TStorage>("-0.00047282220513073084"));
			a[6].real(DQStuff::fromString<TStorage>("0.003816465653840016"));		a[6].imag(DQStuff::fromString<TStorage>("0.017839810396560574"));
			a[7].real(DQStuff::fromString<TStorage>("0.12124105653274578"));		a[7].imag(DQStuff::fromString<TStorage>("-0.12327042473830248"));
			a[8].real(DQStuff::fromString<TStorage>("-0.9774980792734348"));		a[8].imag(DQStuff::fromString<TStorage>("-0.1877130220537587"));
			a[9].real(DQStuff::fromString<TStorage>("1.3432866123333178"));			a[9].imag(DQStuff::fromString<TStorage>("3.2034715228495942"));
			a[10].real(DQStuff::fromString<TStorage>("4.072408546157305"));			a[10].imag(DQStuff::fromString<TStorage>("-6.123755543580666"));
			a[11].real(DQStuff::fromString<TStorage>("-9.442699917778205"));		a[11].imag(DQStuff::fromString<TStorage>("0.0"));
			a[12].real(DQStuff::fromString<TStorage>("4.072408620272648"));			a[12].imag(DQStuff::fromString<TStorage>("6.123755841848161"));
			a[13].real(DQStuff::fromString<TStorage>("1.3432860877712938"));		a[13].imag(DQStuff::fromString<TStorage>("-3.2034712658530275"));
			a[14].real(DQStuff::fromString<TStorage>("-0.9774985292598916"));		a[14].imag(DQStuff::fromString<TStorage>("0.18771238018072134"));
			a[15].real(DQStuff::fromString<TStorage>("0.1212417070363373"));		a[15].imag(DQStuff::fromString<TStorage>("0.12326987628935386"));
			a[16].real(DQStuff::fromString<TStorage>("0.0038169724770333343"));		a[16].imag(DQStuff::fromString<TStorage>("-0.017839242222443888"));
			a[17].real(DQStuff::fromString<TStorage>("-0.0010756025812659208"));	a[17].imag(DQStuff::fromString<TStorage>("0.0004731874917343858"));
			a[18].real(DQStuff::fromString<TStorage>("0.000014713754789095218"));	a[18].imag(DQStuff::fromString<TStorage>("0.000031358475831136815"));
			a[19].real(DQStuff::fromString<TStorage>("-2.659323898804944e-6"));		a[19].imag(DQStuff::fromString<TStorage>("-0.000011341571201752273"));
			a[20].real(DQStuff::fromString<TStorage>("3.6970377676364553e-6"));		a[20].imag(DQStuff::fromString<TStorage>("-6.517457477594937e-7"));
			a[21].real(DQStuff::fromString<TStorage>("3.883933649142257e-9"));		a[21].imag(DQStuff::fromString<TStorage>("9.128496023863376e-7"));
			a[22].real(DQStuff::fromString<TStorage>("-1.0816457995911385e-7"));	a[22].imag(DQStuff::fromString<TStorage>("-2.954309729192276e-8"));

#elif 1
			mu = -4.315321510875024;

			/*
			 * Coefficients from Terry
			 */
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
#else
			/*
			 * Poles which are forced to be symmetric and computed with optimizer
			 * See output_rexi_gaussian_approx_output_dpsfloat128_N23_eps1e-15_minimizeeps1e-15_MCG_gs1.0_2017-06-12_09_35_50.444617.log
			 */

			mu = -4.315321510875044452859583543614e+00;

			a = {
					{	-1.083110282911149456307994863938e-07,	2.862532517802907476007183299252e-08	},
					{	1.123573978342406675797923146145e-08,	-9.116935703342052795387939236049e-07	},
					{	3.685704547609694109520754590203e-06,	6.795370931272857571990879052559e-07	},
					{	-2.729164855027209818832999513316e-06,	1.129902698668245529479634003955e-05	},
					{	1.481616616464353772058496622188e-05,	-3.149316224546558226843676053797e-05	},
					{	-1.075389654802298916927427718804e-03,	-4.730048484188074443303195781851e-04	},
					{	3.816719065433554509275682065095e-03,	1.783952630952095538829915710721e-02	},
					{	1.212413817845429536701473693938e-01,	-1.232701505138063119426661273792e-01	},
					{	-9.774983042666550714372419861320e-01,	-1.877127011172174952946534176590e-01	},
					{	1.343286350052321553860679159698e+00,	3.203471394351331102967606057064e+00	},
					{	4.072408583214998323285271908389e+00,	-6.123755692714398790599261701573e+00	},
					{	-9.442699917778186957662001077551e+00,	0.000000000000000000000000000000e+00	},
			};

			for (int i = 0; i < L; i++)
			{
				a[2*L-i].real(a[i].real());
				a[2*L-i].imag(-a[i].imag());
			}
#endif
		}
		else if (L == -11)
		{
			/**
			 * Recomputed coefficients with higher accuracy
			 */
			mu = -4.315321510875024024755930440733209252357;

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
	Tevaluation evalGaussian(
			Tevaluation x,	///< x-coefficient for Gaussian basis function
			Tevaluation h	///< h-coefficient for Gaussian basis function
	)
	{
		return DQStuff::exp(-(x*x)/(4.0*h*h))/sqrtpi4;
	}


	/**
	 * evaluate approximation of Gaussian basis function
	 *
	 * with sum of complex rational functions
	 */
	Tevaluation approxGaussian(
			Tevaluation x,	///< x-coefficient for Gaussian basis function
			Tevaluation h	///< h-coefficient for Gaussian basis function
	)
	{
		// scale x, since it depends linearly on h:
		// x^2 ~ h^2
		x /= h;

		Tevaluation sum = 0;

		for (int l = 0; l < 2*L+1; l++)
		{
			int j = l-L;

			std::complex<Tevaluation> aEval(a[l].real(), a[l].imag());

			// WORKS with max error 7.15344e-13
			sum += (aEval/(complex(0, x) + (Tevaluation)mu + complex(0, j))).real();
		}

		return sum;
	}


	/**
	 * evaluate approximation of Gaussian basis function
	 *
	 * with sum of complex rational functions
	 */
	std::complex<Tevaluation> approxGaussian_Complex(
			Tevaluation x,	///< x-coefficient for Gaussian basis function
			Tevaluation h	///< h-coefficient for Gaussian basis function
	)
	{
		// scale x, since it depends linearly on h:
		// x^2 ~ h^2
		x /= h;

		std::complex<Tevaluation> sum = 0;

		for (int l = 0; l < 2*L+1; l++)
		{
			int j = l-L;

			sum += (a[l]/(complex(0, x) + (Tevaluation)mu + complex(0, j)));
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
