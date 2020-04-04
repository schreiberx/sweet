/*
 * SPHTestSolutions.hpp
 *
 *  Created on: 24 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SPHTESTSOLUTIONSGAUSSIAN_HPP_
#define SRC_INCLUDE_SPHTESTSOLUTIONSGAUSSIAN_HPP_


class SphereTestSolutions_Gaussian
{
	double exp_fac = 10.0;
#if 0
	double center_lon = M_PI/3;
	double center_lat = M_PI/2;
#else
	double center_lon = M_PI/5;
	double center_lat = M_PI/4;
//	double center_lat = 0;
#endif


public:
	SphereTestSolutions_Gaussian(
			double i_lon = M_PI/5.0,
			double i_lat = M_PI/4.0
	)
	{
		center_lon = i_lon;
		center_lat = i_lat;
	}


public:
	// EXPONENTIAL
	void test_function__grid_gaussian(double lon, double mu, double &o_data)
	{
		// https://en.wikipedia.org/wiki/Great-circle_distance
		// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
		// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

		double phi1 = asin(mu);
		double phi2 = center_lat;
		double lambda1 = lon;
		double lambda2 = center_lon;

		double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		o_data = exp(-d*d*exp_fac);
	};

	void test_function_phi__grid_phi(double lon, double phi1, double &o_data)
	{
		// https://en.wikipedia.org/wiki/Great-circle_distance
		// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
		// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

		double phi2 = center_lat;
		double lambda1 = lon;
		double lambda2 = center_lon;

		double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		o_data = exp(-d*d*exp_fac);
	};

	void correct_result_diff_lambda__grid_gaussian(double lon, double mu, double &o_data)
	{
		double phi1 = asin(mu);
		double phi2 = center_lat;
		double lambda1 = lon;
		double lambda2 = center_lon;

		double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		// http://www.wolframalpha.com/input/?i=diff(exp(-pow(arccos(sin(a)*sin(b)+%2B+cos(a)*cos(b)*cos(c-d)),2)*A),+c)
		o_data = -
				(exp_fac*2.0*cos(phi1)*cos(phi2)*sin(lambda1-lambda2))
				*d
				*exp(-d*d*exp_fac)
				/sqrt(1-pow(cos(phi1)*cos(phi2)*cos(lambda1-lambda2) + sin(phi1)*sin(phi2), 2.0))
			;
	};

	void correct_result_diff_phi__grid_gaussian(double lon, double mu, double &o_data)
	{
		// OK
		double phi1 = asin(mu);
		double phi2 = center_lat;
		double lambda1 = lon;
		double lambda2 = center_lon;

		// http://www.wolframalpha.com/input/?i=diff(exp(-pow(arccos(sin(a)*sin(b)+%2B+cos(a)*cos(b)*cos(c-d)),2)*A),+a)
		double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		o_data = exp_fac*2.0
				*d
				*(cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda1-lambda2))
				*exp(-d*d*exp_fac)
				/sqrt(1.0-pow(cos(phi1)*cos(phi2)*cos(lambda1-lambda2) + sin(phi1)*sin(phi2), 2.0))
			;
	};

	void correct_result_diff_mu__grid_gaussian(double lon, double mu, double &o_data)
	{
		correct_result_diff_phi__grid_gaussian(lon, mu, o_data);
		o_data /= sqrt(1.0-mu*mu);
	};

	void correct_result_grad_lambda__grid_gaussian(double lon, double mu, double &o_data)
	{
		correct_result_diff_lambda__grid_gaussian(lon, mu, o_data);
		o_data /= sqrt(1.0-mu*mu);
	};

	void correct_result_mu__grid_gaussian(double lon, double mu, double &o_data)
	{
		test_function__grid_gaussian(lon, mu, o_data);
		o_data *= mu;
	};
	void correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(double lon, double mu, double &o_data)
	{
		correct_result_diff_mu__grid_gaussian(lon, mu, o_data);
		o_data *= (1.0-mu*mu);
	};
	void correct_result_grad_phi__grid_gaussian(double lon, double mu, double &o_data)
	{
		correct_result_diff_phi__grid_gaussian(lon, mu, o_data);
	};

	void correct_result_div_lambda__grid_gaussian(double lon, double mu, double &o_data)
	{
		correct_result_grad_lambda__grid_gaussian(lon, mu, o_data);
	};

	void correct_result_div_mu__grid_gaussian(double lon, double mu, double &o_data)
	{
		double phi1 = asin(mu);
		double phi2 = center_lat;
		double lambda1 = lon;
		double lambda2 = center_lon;

		// http://www.wolframalpha.com/input/?i=diff(cos(a)*exp(-pow(arccos(sin(a)*sin(b)+%2B+cos(a)*cos(b)*cos(c-d)),2)*A),+a)
		double da = sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2);
		double d = acos(da);

		o_data = exp(-exp_fac*d*d)*
				(
					2.0*exp_fac*cos(phi1)
						*d
						*(cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda1-lambda2))
						/sqrt(1.0-da*da)

					-sin(phi1)
				)
			;
		o_data /= cos(phi1);
	};
};


#endif /* SRC_INCLUDE_SPHTESTSOLUTIONSGAUSSIAN_HPP_ */
