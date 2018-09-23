/*
 * SphereStaggering.hpp
 *
 *  Created on: 26 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_SPHERESTAGGERING_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_SPHERESTAGGERING_HPP_



/**
 * Class which defines staggering for all variables
 */
class SphereStaggering
{
public:
	double h[2] = {-0.5,-0.5};
	double u[2] = {-0.5,-0.5};
	double v[2] = {-0.5,-0.5};

	char staggering_type = 'a';


public:
	SphereStaggering()
	{
		setup_a_staggering();
	}




public:
	void setup_a_staggering()
	{
		staggering_type = 'a';

		h[0] = -0.5;
		h[1] = -0.5;

		u[0] = -0.5;
		u[1] = -0.5;

		v[0] = -0.5;
		v[1] = -0.5;
	}



	/*
	 *              ^
	 *              |
	 *       ______v0,1_____
	 *       |             |
	 *       |			   |
	 *       |   (0.5,0.5) |
	 *  u0,0 |->  H/P0,0   |u1,0 ->
	 *(0,0.5)|	           |
	 *       |      ^      |
	 *   q0,0|______|______|
	 * (0,0)      v0,0
	 *           (0.5,0)
	 *
	 * These staggering should be used when interpolating from a staggered variable
	 * If interpolating from A grid to C staggered, use negative of displacements.
	 */
public:
	void setup_c_staggering()
	{
		staggering_type = 'c';

		h[0] = -0.5;
		h[1] = -0.5;

		u[0] = -0.0;
		u[1] = -0.5;

		v[0] = -0.5;
		v[1] = -0.0;
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_SPHERESTAGGERING_HPP_ */
