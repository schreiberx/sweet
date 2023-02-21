/*
 * test_sphere_coordinates.cpp
 *
 *  Created on: 17 Apr 2018
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <iostream>
#include <sweet/ScalarDataArray.hpp>
#include <sweet/SWEETError.hpp>
#include <sweet/SWEETVectorMath.hpp>


void check_error(
		double i_lon_input1,
		double i_lat_input1,
		double i_lon_input2,
		double i_lat_input2
)
{
	double error = std::max(std::abs(i_lon_input1-i_lon_input2), std::abs(i_lat_input1-i_lat_input2));
	double eps = 1e-10;

	if (error > eps)
	{
		std::cout << std::endl;
		std::cout << "input 1: (lon/lat): " << i_lon_input1 << "\t" << i_lon_input1 << std::endl;
		std::cout << "input 2: (lon/lat): " << i_lat_input2 << "\t" << i_lat_input2 << std::endl;

		SWEETError("Error too large!");
	}
}

void check_error(
		double i_x_input1,
		double i_y_input1,
		double i_z_input1,
		double i_x_input2,
		double i_y_input2,
		double i_z_input2
)
{
	double error = std::max(
			std::max(
				std::abs(i_x_input1-i_x_input2),
				std::abs(i_y_input1-i_y_input2)
			),
				std::abs(i_z_input1-i_z_input2)
		);

	double eps = 1e-10;

	if (error > eps)
	{
		std::cout << std::endl;
		std::cout << "input 1: (x,y,z): " << i_x_input1 << "\t" << i_y_input1 << "\t" << i_z_input1 << std::endl;
		std::cout << "input 2: (x,y,z): " << i_x_input2 << "\t" << i_y_input2 << "\t" << i_z_input2 << std::endl;

		SWEETError("Error too large!");
	}
}


int main(int i_argc, char *i_argv[])
{
	double N = 32;
	//double N = 128;

	double dlon = 2.0*M_PI/(N+1);
	double dlat = 2.0*M_PI/(N+1);


	{
		std::cout << "************************************************************" << std::endl;
		std::cout << "* Some information on longitude <=> Cartesian transformations" << std::endl;
		std::cout << "* Used functions: latlon_to_cartesian, cartesian_to_latlon" << std::endl;
		std::cout << "************************************************************" << std::endl;
		for (double a_lon = dlon*0.5; a_lon < 2.0*M_PI; a_lon += 0.3)
		{
			double a_lat = 0;

			double a_x, a_y, a_z;
			SWEETVectorMath::point_latlon_to_cartesian__scalar(a_lon, a_lat, a_x, a_y, a_z);

			double o_lon, o_lat;
			SWEETVectorMath::point_cartesian_to_latlon__scalar(a_x, a_y, a_z, o_lon, o_lat);

			std::cout << "lat/lon: " << a_lon << ", " << a_lat;
			std::cout << "\t=> Cartesian: " << a_x << ", " << a_y << ", " << a_z;
			std::cout << "\t=> lat/lon: " << o_lon << ", " << o_lat;
			std::cout << std::endl;
		}
	}


	{
		std::cout << "************************************************************" << std::endl;
		std::cout << "* Some information on uv velocity <=> Cartesian transformations" << std::endl;
		std::cout << "* Used functions: latlon_to_cartesian, cartesian_to_latlon" << std::endl;
		std::cout << "************************************************************" << std::endl;
		{
			double vel_u = 2.0;
			double vel_v = 0.0;

			std::cout << "* u,v velocity: " << vel_u << ", " << vel_v << std::endl;
			for (double a_lon = 0; a_lon < 2.0*M_PI; a_lon += 0.3)
			{
				double a_lat = 0;

				double vel_x, vel_y, vel_z;
				SWEETVectorMath::velocity_latlon_to_cartesian__scalar(a_lon, a_lat, vel_u, vel_v, vel_x, vel_y, vel_z);

				std::cout << "lat/lon: " << a_lon << ", " << a_lat;
				std::cout << "\t=> u,v velocity: " << vel_u << ", " << vel_v;
				std::cout << "\t=> x,y,z velocity: " << vel_x << ", " << vel_y << ", " << vel_z;
				std::cout << "\t=> len(x,y,z velocity): " << std::sqrt(vel_x*vel_x + vel_y*vel_y + vel_z*vel_z);
				std::cout << std::endl;
			}
		}
	}

	{
		std::cout << "************************************************************" << std::endl;
		std::cout << "* Some information on uv velocity <=> Cartesian transformations" << std::endl;
		std::cout << "* Used functions: latlon_to_cartesian, cartesian_to_latlon" << std::endl;
		std::cout << "************************************************************" << std::endl;
		{
			double vel_u = 0.0;
			double vel_v = 1.0;

			std::cout << "* u,v velocity: " << vel_u << ", " << vel_v << std::endl;
			for (double a_lon = 0; a_lon < 2.0*M_PI; a_lon += 0.3)
			{
				double a_lat = 0.5*M_PI*0.5;

				double vel_x, vel_y, vel_z;
				SWEETVectorMath::velocity_latlon_to_cartesian__scalar(a_lon, a_lat, vel_u, vel_v, vel_x, vel_y, vel_z);

				std::cout << "lat/lon: " << a_lon << ", " << a_lat;
				std::cout << "\t=> u,v velocity: " << vel_u << ", " << vel_v;
				std::cout << "\t=> x,y,z velocity: " << vel_x << ", " << vel_y << ", " << vel_z;
				std::cout << "\t=> len(x,y,z velocity): " << std::sqrt(vel_x*vel_x + vel_y*vel_y + vel_z*vel_z);
				std::cout << std::endl;
			}
		}
	}


	{
		std::cout << "************************************************************" << std::endl;
		std::cout << "* Tested functions: latlon_normalize" << std::endl;
		std::cout << "************************************************************" << std::endl;
		int k = 0;
		{
			std::cout << "Test " << k++ << std::endl;
			double a_lon = -2.0;
			double a_lat = 0;

			SWEETVectorMath::point_latlon_normalize__scalar(a_lon, a_lat);

			double o_lon = -2.0+M_PI*2.0;
			double o_lat = 0;
			check_error(a_lon, a_lat, o_lon, o_lat);
		}
		{
			std::cout << "Test " << k++ << std::endl;
			double a_lon = 8.0;
			double a_lat = 0;

			SWEETVectorMath::point_latlon_normalize__scalar(a_lon, a_lat);

			double o_lon = 8.0-M_PI*2.0;
			double o_lat = 0;
			check_error(a_lon, a_lat, o_lon, o_lat);
		}
		{
			std::cout << "Test " << k++ << std::endl;
			double a_lon = 8.0;
			double a_lat = 0;

			SWEETVectorMath::point_latlon_normalize__scalar(a_lon, a_lat);

			double o_lon = 8.0-M_PI*2.0;
			double o_lat = 0;
			check_error(a_lon, a_lat, o_lon, o_lat);
		}
		{
			std::cout << "Test " << k++ << std::endl;
			double a_lon = 2.0;
			double a_lat = M_PI*0.5+0.3;

			SWEETVectorMath::point_latlon_normalize__scalar(a_lon, a_lat);

			double o_lon = 2.0+M_PI;
			double o_lat = M_PI*0.5-0.3;
			check_error(a_lon, a_lat, o_lon, o_lat);
		}
		{
			std::cout << "Test " << k++ << std::endl;
			double a_lon = -2.0;
			double a_lat = -M_PI*0.5-0.3;

			SWEETVectorMath::point_latlon_normalize__scalar(a_lon, a_lat);

			double o_lon = -2.0+M_PI;
			double o_lat = -M_PI*0.5+0.3;
			check_error(a_lon, a_lat, o_lon, o_lat);
		}
	}


	{
		std::cout << "************************************************************" << std::endl;
		std::cout << "* Tests for various longitudes and latitudes" << std::endl;
		std::cout << "* Tested functions: latlon_to_cartesian, cartesian_to_latlon" << std::endl;
		std::cout << "************************************************************" << std::endl;

		for (double a_lon = dlon*0.5; a_lon < 2.0*M_PI; a_lon += dlon)
		{
			std::cout << "Testing for longitude " << a_lon << std::endl;

			for (double a_lat = -M_PI*0.5+dlat*0.5; a_lat <= M_PI*0.5; a_lat += dlat)
			{
				double o_lon, o_lat;
				double a_x, a_y, a_z;
				SWEETVectorMath::point_latlon_to_cartesian__scalar(a_lon, a_lat, a_x, a_y, a_z);
				SWEETVectorMath::point_cartesian_to_latlon__scalar(a_x, a_y, a_z, o_lon, o_lat);

				check_error(o_lon, o_lat, a_lon, a_lat);
			}
		}
	}


	{
		std::cout << "************************************************************" << std::endl;
		std::cout << "* Tests for various lon/lat velocity conversions" << std::endl;
		std::cout << "* Tested functions: velocity_latlon_to_cartesian_scalar, velocity_cartesian_to_latlon_scalar" << std::endl;
		std::cout << "************************************************************" << std::endl;

		for (double a_lon = dlon*0.5; a_lon < 2.0*M_PI; a_lon += dlon)
		{
			std::cout << "Testing for longitude " << a_lon << std::endl;

			for (double a_lat = -M_PI*0.5+dlat*0.5; a_lat <= M_PI*0.5; a_lat += dlat)
			{
				for (double a_u = -2.0; a_u <= 2.0; a_u += 0.5)
				{
					for (double a_v = -2.0; a_v <= 2.0; a_v += 0.25)
					{
						double a_V_x, a_V_y, a_V_z;
						SWEETVectorMath::velocity_latlon_to_cartesian__scalar(a_lon, a_lat, a_u, a_v, a_V_x, a_V_y, a_V_z);

						double o_u, o_v;
						SWEETVectorMath::velocity_cartesian_to_latlon__scalar(a_lon, a_lat, a_V_x, a_V_y, a_V_z, o_u, o_v);

						check_error(a_u, a_v, o_u, o_v);
					}
				}
			}
		}
	}


	{
		std::cout << "************************************************************" << std::endl;
		std::cout << "* Tests for Cartesian point rotations" << std::endl;
		std::cout << "* Tested function: point_rotate_3d__scalar" << std::endl;
		std::cout << "************************************************************" << std::endl;

		double testcases[][10] =
			{
				/*
				 * starting point, angle, axis, end point
				 *
				 * Point with your index finger along the rotation axis and then rotate your hand to the left.
				 */

				/*
				 * angle = 0
				 */
				// Starting point placed towards the reader of this code! :-)
				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ 0.0,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ 0.0,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 0.0, 1.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ 0.0,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 0.0, 0.0, 1.0},

				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ 0.0,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ 0.0,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 0.0, 1.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ 0.0,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 0.0, 0.0, 1.0},

				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ 0.0,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ 1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ 0.0,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ 0.0, 1.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ 0.0,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ 0.0, 0.0, 1.0},

				/*
				 * angle = pi
				 */
				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ M_PI,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ M_PI,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 0.0, -1.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ M_PI,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 0.0, 0.0, -1.0},

				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ M_PI,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ -1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ M_PI,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 0.0, 1.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ M_PI,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 0.0, 0.0, -1.0},

				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ M_PI,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ -1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ M_PI,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ 0.0, -1.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ M_PI,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ 0.0, 0.0, 1.0},


				/*
				 * angle = pi*0.5,
				 *
				 * Right-hand coordinate system:
				 * Use your right hand, point with the thumb along the rotation axis.
				 * Your fingers point out the rotation direction.
				 */
				/*
				 * angle = pi*0.5, rotation axis (1,0,0)
				 */
				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ M_PI*0.5,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ M_PI*0.5,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 0.0, 0.0, 1.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ M_PI*0.5,	/* axis */ 1.0, 0.0, 0.0,	/* end point */ 0.0, -1.0, 0.0},

				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ M_PI*0.5,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 0.0, 0.0, -1.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ M_PI*0.5,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 0.0, 1.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ M_PI*0.5,	/* axis */ 0.0, 1.0, 0.0,	/* end point */ 1.0, 0.0, 0.0},

				{/* starting point */ 1.0, 0.0, 0.0,	/* angle */ M_PI*0.5,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ 0.0, 1.0, 0.0},
				{/* starting point */ 0.0, 1.0, 0.0,	/* angle */ M_PI*0.5,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ -1.0, 0.0, 0.0},
				{/* starting point */ 0.0, 0.0, 1.0,	/* angle */ M_PI*0.5,	/* axis */ 0.0, 0.0, 1.0,	/* end point */ 0.0, 0.0, 1.0},

			};

		int num_test_cases = sizeof(testcases)/sizeof(testcases[1]);

		for (int i = 0; i < num_test_cases; i++)
		{
			std::cout << "Testcase " << i << std::endl;

			double o_P_x, o_P_y, o_P_z;
			SWEETVectorMath::point_rotate_3d__scalar(
					testcases[i][0],	// start
					testcases[i][1],
					testcases[i][2],
					testcases[i][3],	// angle
					testcases[i][4],	// rotation vector
					testcases[i][5],
					testcases[i][6],
					o_P_x,
					o_P_y,
					o_P_z
				);

			std::cout << " +       starting point: " << testcases[i][0] << ", " << testcases[i][1] << ", " << testcases[i][2] << std::endl;
			std::cout << " +       rotation angle: " << testcases[i][3] << std::endl;
			std::cout << " +      rotation vector: " << testcases[i][4] << ", " << testcases[i][5] << ", " << testcases[i][6] << std::endl;
			std::cout << " +     (test) end point: " << testcases[i][7] << ", " << testcases[i][8] << ", " << testcases[i][9] << std::endl;
			std::cout << " + (computed) end point: " << o_P_x << ", " << o_P_y << ", " << o_P_z << std::endl;

			check_error(o_P_x, o_P_y, o_P_z, testcases[i][7], testcases[i][8], testcases[i][9]);
			std::cout << std::endl;
		}
	}




	std::cout << "************************************************************" << std::endl;
	std::cout << "* Tests passed" << std::endl;
	std::cout << "************************************************************" << std::endl;

	return 0;
}
