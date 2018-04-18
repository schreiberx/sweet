/*
 * test_sphere_coordinates.cpp
 *
 *  Created on: 17 Apr 2018
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <sweet/sphere/SphereDataSemiLagrangian.hpp>
#include <sweet/FatalError.hpp>




int main(int i_argc, char *i_argv[])
{
	//double N = 16;
	double N = 1024;
	double dlon = 1.0/(N+1);
	double dlat = 1.0/(N+1);
	double eps = 1e-10;

	for (double lon = dlon*0.5; lon < 2.0*M_PI; lon += dlon)
	{
//		for (double lat = -M_PI*0.5+dlat*0.5; lat <= M_PI*0.5; lat += dlat)
		double lat = 0;
		{
			//std::cout << "input: " << lon << "\t" << lat << std::endl;
			double cart_coord[3];

			SphereDataSemiLagrangian::angleToCartCoord(lon, lat, cart_coord);

			double o_lon, o_lat;
			SphereDataSemiLagrangian::cartToAngleCoord(cart_coord, &o_lon, &o_lat);

			double error = std::max(std::abs(o_lon-lon), std::abs(o_lat-lat));

			//std::cout << "cart: " << cart_coord[0] << "\t" << cart_coord[1] << "\t" << cart_coord[2] << std::endl;
			//std::cout << "output: " << o_lon << "\t" << o_lat << std::endl;

			if (error > eps)
			{
				std::cout << std::endl;
				std::cout << "input: " << lon << "\t" << lat << std::endl;
				std::cout << "cart: " << cart_coord[0] << "\t" << cart_coord[1] << "\t" << cart_coord[2] << std::endl;
				std::cout << "output: " << o_lon << "\t" << o_lat << std::endl;
				FatalError("Error too large!");
			}
		}
	}

	std::cout << "Tests passed" << std::endl;

	return 0;
}
