/*
 * BenchmarkFlowOverMountain.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SPHERE_BENCHMARKFLOWOVERMOUNTAIN_HPP_
#define SRC_INCLUDE_BENCHMARKS_SPHERE_BENCHMARKFLOWOVERMOUNTAIN_HPP_


class BenchmarkFlowOverMountain
{
public:
	static
	void setup_topography(
			SphereDataPhysical          &o_h_topo,
			SimulationVariables &i_simVars,
			double i_R            = M_PI/9.,
			double i_h_topo_0     = 2000.,  
			double i_center_lon   = 3.*M_PI/2.,
			double i_center_lat   = M_PI/6.
	)
	{
	        const double center_lat = i_center_lat;
		const double center_lon = i_center_lon;

		auto topography = [&](double lon, double mu, double &o_data)
		{

			const double phi1    = asin(mu);
			const double phi2    = center_lat;
			const double lambda1 = lon;
			const double lambda2 = center_lon;

			const double r_squared = std::min( i_R*i_R, (phi1-phi2)*(phi1-phi2) + (lambda1-lambda2)*(lambda1-lambda2) );
			
			o_data = i_h_topo_0 * ( 1. - sqrt(r_squared) / i_R );
		};

		o_h_topo.physical_update_lambda_gaussian_grid(topography);
	}

};


#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_BENCHMARKFLOWOVERMOUTAIN_HPP_ */
