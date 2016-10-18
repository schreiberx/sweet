/*
 * sh_helpers.hpp
 *
 *  Created on: 10 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHHELPER_HPP_
#define SPHHELPER_HPP_

#include <array>

class SPHHelper
{
public:
	static void angles_to_cart(
			double lon,
			double lat,
			std::array<double,3> &o_vec
	)
	{
		o_vec[0] = std::cos(lon)*cos(lat);
		o_vec[1] = std::sin(lon)*cos(lat);
		o_vec[2] = std::sin(lat);
	};

public:
	static void angles_gauss_to_cart(
			double lon,
			double mu,
			std::array<double,3> &o_vec
	)
	{
		double cos_lat = std::sqrt(1-mu*mu);

		o_vec[0] = std::cos(lon)*cos_lat;
		o_vec[1] = std::sin(lon)*cos_lat;
		o_vec[2] = mu;
	};

public:
	static double cart_dist(
			const std::array<double,3> &i_vec1,
			const std::array<double,3> &i_vec2
	)
	{
		double d0 = i_vec1[0] - i_vec2[0];
		double d1 = i_vec1[1] - i_vec2[1];
		double d2 = i_vec1[2] - i_vec2[2];

		return std::sqrt(d0*d0 + d1*d1 + d2*d2);
	}
};





#endif /* SPHHELPER_HPP_ */
