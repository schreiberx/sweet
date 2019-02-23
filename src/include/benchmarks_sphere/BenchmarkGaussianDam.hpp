/*
 * BenchmarkGaussianDam.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SPHERE_BENCHMARKGAUSSIANDAM_HPP_
#define SRC_INCLUDE_BENCHMARKS_SPHERE_BENCHMARKGAUSSIANDAM_HPP_


class BenchmarkGaussianDam
{

public:
	static
	void setup_initial_conditions_gaussian(
			SphereData_Physical &o_h,
			SimulationVariables &i_simVars,
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{

		double center_lat = i_center_lat;
		double center_lon = i_center_lon;

		auto initial_condition_h = [&](double lon, double mu, double &o_data)
		{
			// https://en.wikipedia.org/wiki/Great-circle_distance
			// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
			// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			o_data = std::exp(-d*d*i_exp_fac)*0.1*i_simVars.sim.h0 + i_simVars.sim.h0;
		};

		o_h.physical_update_lambda_gaussian_grid(initial_condition_h);
	}



public:
	static
	void setup_initial_conditions_gaussian_normalized(
			SphereData_Physical &o_data,
			const SimulationVariables &i_simVars,
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{
		o_data.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					// https://en.wikipedia.org/wiki/Great-circle_distance
					// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
					// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

					double phi1 = asin(mu);
					double phi2 = i_center_lat;
					double lambda1 = lon;
					double lambda2 = i_center_lon;

					double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

					o_data = std::exp(-d*d*i_exp_fac);
				}
		);
	}



public:
	static
	void setup_initial_conditions_gaussian(
			SphereData_Physical &o_h,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v,
			SimulationVariables &i_simVars,
			double i_center_lon = M_PI/3,
			double i_center_lat = M_PI/3,
			double i_exp_fac = 10.0
	)
	{
		setup_initial_conditions_gaussian(o_h, i_simVars, i_center_lon, i_center_lat, i_exp_fac);

		o_u.physical_set_zero();
		o_v.physical_set_zero();
	}


};


#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_BENCHMARKGAUSSIANDAM_HPP_ */
