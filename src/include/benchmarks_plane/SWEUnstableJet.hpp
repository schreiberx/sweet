/*
 * SWEUnstableJet.hpp
 *
 *  Created on: 03 Nov 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef SWE_PLANE_UNSTABLEJET_HPP_
#define SWE_PLANE_UNSTABLEJET_HPP_


#include <stdlib.h>
#include <sweet/sweetmath.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>


/**
 * Implement unstable jet initial conditions
 *
 * Mimics Spherical Mountain Wave test case 5
 *
 *
 **/
class SWEUnstableJet
{
	SimulationVariables &simVars;

	PlaneOperators &op;

	double f = simVars.sim.f0;
	double g = simVars.sim.gravitation;
	double sx = simVars.sim.domain_size[0];
	double sy = simVars.sim.domain_size[1];


	/*
	 * The depth function is numerically integrated to ensure
	 * balanced initial conditions for the jet
	 */
	double depth(
			double x,
			double y
			)
	{
		//TODO

				return 1.0;

	}

	double u(
			double x,
			double y
			)
	{
		//TODO
				return std::pow(std::sin(2.0*M_PI*y/sy),20);

	}

	void setup_depth(
			PlaneData &o_depth
	)
	{

		for (int j = 0; j < simVars.disc.res_physical[1]; j++)
		{
			for (int i = 0; i < simVars.disc.res_physical[0]; i++)
			{

				// h - lives in the center of the cell
				double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

				o_depth.p_physical_set(j, i, depth(x, y));
			}
		}

	}

	void setup_velocity(
			PlaneData &o_u,
			PlaneData &o_v
	)
	{

		o_v.physical_set_zero();

		for (int j = 0; j < simVars.disc.res_physical[1]; j++)
		{
			for (int i = 0; i < simVars.disc.res_physical[0]; i++)
			{

				// h - lives in the center of the cell
				double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

				o_u.p_physical_set(j, i, u(x, y));
			}
		}

	}

public:
	SWEUnstableJet(
		SimulationVariables &io_simVars,
		PlaneOperators &io_op
	)	:
		simVars(io_simVars),
		op(io_op)
	{
	}

	void setup(
			PlaneData &o_h,
			PlaneData &o_u,
			PlaneData &o_v
	)
	{
		/*
		 * Setup velocities
		 */
		setup_velocity(o_u,o_v);

		/*
		 * Setup depth function
		 * based on velocities
		 */
		setup_depth(o_h);
		//o_h.file_physical_saveData_ascii("ouput_depth");

		/*
		 * Add perturbation to depth
		 */

	}


};


#endif
