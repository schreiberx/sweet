/*
 * SWEUnstableJetFast.hpp
 *
 *  Created on: 05 Mar 2018
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef SWE_PLANE_UNSTABLEJETFAST_HPP_
#define SWE_PLANE_UNSTABLEJETFAST_HPP_


#include <stdlib.h>

#include <sweet/sweetmath.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <libmath/GaussQuadrature.hpp>


/**
 * Implement unstable jet initial conditions
 *
 * Mimics Spherical Mountain Wave test case 5
 *
 *
 **/
class SWE_bench_UnstableJetFast
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
	 *
	 */
	double depth(
			double x,
			double y
	)
	{

		return -(f/g)*sy*GaussQuadrature::integrate5_intervals_adaptive_linear<double>(0, y, u_fun, 10e-13);
		//return -(f/g)*GaussQuadrature::integrate5_intervals<double>(0, y, u_fun, 200);
	}

	/*
	 * Velocity
	 * On (x,y) \in [0,1]x[0,1]
	 */
	static double u_fun(
			double y
	)
	{
		 //power has to be odd to ensure periodicity
		// the larger the thiner the jet
		// Max speed is 50m/s
		return 300.0*std::pow(std::sin(2.0*M_PI*y), 81);
	}

	double u(
			double x,
			double y
	)
	{
		return u_fun(y);

	}

	/*
	 * Depth perturbation (gaussian bumps)
	 * On (x,y) \in [0,1]x[0,1]
	 */
	double bump(
			double x,
			double y
	)
	{
		//double radius = simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
		double radius = 1.0; //simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
		double factor = 1000.0;


		// Gaussian Bump top
		double dx = x-0.85;
		double dy = y-0.75;

		dx /= radius;
		dy /= radius;

		double exp1 = std::exp(-factor*(dx*dx + dy*dy));

		// Gaussian Bump bottom
		dx = x-0.15;
		dy = y-0.25;

		dx /= radius;
		dy /= radius;

		double exp2 = std::exp(-factor*(dx*dx + dy*dy));

		double pert = 0.01*simVars.sim.h0;
		//double pert = 0.000;

		return (pert)*(exp1+exp2);

	}

	void setup_depth(
			PlaneData &o_depth
	)
	{
		// First set for the first column (one vertical slice)

		for (int j = 0; j < simVars.disc.res_physical[1]; j++)
		{
			int i = 0;
			double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0]); //*simVars.sim.domain_size[0];
			double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1]); //*simVars.sim.domain_size[1];

			o_depth.p_physical_set(j, i, depth(x, y));
		}

		//Now set for other "x" and add bump
		for (int j = 0; j < simVars.disc.res_physical[1]; j++)
		{
			for (int i = 1; i < simVars.disc.res_physical[0]; i++)
			{

				// h - lives in the center of the cell
				// (x,y) \in [0,1]x[0,1]
				double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0]); //*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1]); //*simVars.sim.domain_size[1];

				o_depth.p_physical_set(j, i, o_depth.p_physical_get(j, 0) + bump(x,y));

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

				// (u,v) - lives in the center of the cell
				double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0]); //*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1]); //*simVars.sim.domain_size[1];
				// (x,y) \in [0,1]x[0,1]
				o_u.p_physical_set(j, i, u(x, y));
			}
		}
	}

public:
	SWE_bench_UnstableJetFast(
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
		std::cout<< "Generating Unstable Jet initial conditions.";
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
		std::cout<< "   Done! " << std::endl;
	}


};


#endif
