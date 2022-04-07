/*
 * SWEUnstableJetFast.hpp
 *
 *  Created on: 05 Mar 2018
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef SWE_PLANE_UNSTABLEJETADV_HPP_
#define SWE_PLANE_UNSTABLEJETADV_HPP_


#include <stdlib.h>

#include <cmath>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>
#include <libmath/GaussQuadrature.hpp>


/**
 * Implement unstable jet initial conditions
 *
 * Mimics Spherical Mountain Wave test case 5
 *
 *
 **/
class SWE_bench_UnstableJetAdv
{
	SimulationVariables &simVars;

	PlaneOperators &op;

	double f = simVars.sim.plane_rotating_f0;
	double g = simVars.sim.gravitation;
	double sx = simVars.sim.plane_domain_size[0];
	double sy = simVars.sim.plane_domain_size[1];


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
	  //return simVars.sim.h0;
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
		return 100.0*std::pow(std::sin(2.0*M_PI*y), 81);
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

		//default is 1%
		double pert = 0.01*simVars.sim.h0;
		//double pert = 0.000;

		return (pert)*(exp1+exp2);

	}

	void setup_depth(
			PlaneData_Spectral &o_depth
	)
	{
		// First set for the first column (one vertical slice)

		PlaneData_Physical depth_phys(o_depth.planeDataConfig);

		for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
		{
			int i = 0;
			double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0]); //*simVars.sim.domain_size[0];
			double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1]); //*simVars.sim.domain_size[1];

			depth_phys.physical_set_value(j, i, depth(x, y));
		}

		//Now set for other "x" and add bump
		for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
		{
			for (int i = 1; i < simVars.disc.space_res_physical[0]; i++)
			{

				// h - lives in the center of the cell
				// (x,y) \in [0,1]x[0,1]
				double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0]); //*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1]); //*simVars.sim.domain_size[1];

				depth_phys.physical_set_value(j, i, depth_phys.physical_get(j, 0) + bump(x,y));

			}
		}


		o_depth.loadPlaneDataPhysical(depth_phys);

	}

	void setup_velocity(
			PlaneData_Spectral &o_u,
			PlaneData_Spectral &o_v
	)
	{

		o_v.spectral_set_zero();

		PlaneData_Physical u_phys(o_u.planeDataConfig);


		for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
		{
			for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
			{

				// (u,v) - lives in the center of the cell
				double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0]); //*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1]); //*simVars.sim.domain_size[1];
				// (x,y) \in [0,1]x[0,1]
				u_phys.physical_set_value(j, i, u(x, y));
			}
		}

		o_u.loadPlaneDataPhysical(u_phys);
	}

public:
	SWE_bench_UnstableJetAdv(
		SimulationVariables &io_simVars,
		PlaneOperators &io_op
	)	:
		simVars(io_simVars),
		op(io_op)
	{
	}

	void setup(
			PlaneData_Spectral &o_h,
			PlaneData_Spectral &o_u,
			PlaneData_Spectral &o_v
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
