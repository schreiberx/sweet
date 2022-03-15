/*
 * SWEMergeVortex.hpp
 *
 *  Created on: 01 Nov 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef SWE_PLANE_MERGEVORTEX_HPP_
#define SWE_PLANE_MERGEVORTEX_HPP_


#include <stdlib.h>
#include <cmath>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>


/**
 * Implement merging vortex initial conditions
 * See Energy- and enstrophy-conserving schemes for the shallow-water
 * equations, based on mimetic finite elements
 * Andrew T. T. McRae and Colin J. Cotter
 *
 * IMPORTANT: TO BE USED IN [0,1]x[0,1] domain
 *
 *  To match paper use:
 * f = 1
 * g = 1
 * h0 = 1
 * [0,1]x[0,1
 **/
class SWE_bench_MergeVortex
{
	SimulationVariables &simVars;

	PlaneOperators &op;

	double f = simVars.sim.plane_rotating_f0;
	double g = simVars.sim.gravitation;
	double sx = simVars.sim.plane_domain_size[0];
	double sy = simVars.sim.plane_domain_size[1];


	double stream(
			double x,
			double y
			)
	{

		//double radius = simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
		double radius = 1; //simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
		double factor = 500.0;

		//Andrew's parameters
		double k = 0.03;
		double a = 10.0;
		//my parameters
		//a=20.0;
		factor=(a*a);

		// Gaussian Vortice 1
		//double dx = x-0.45*sx;
		//double dy = y-0.5*sy;
		double dx = x-0.425*sx;
		double dy = y-0.5*sy;

		dx /= radius;
		dy /= radius;

		double exp1 = std::exp(-factor*(dx*dx + dy*dy));

		// Gaussian Vortice 2
		//dx = x-0.55*sx;
		//dy = y-0.5*sy;
		dx = x-0.575*sx;
		dy = y-0.5*sy;

		dx /= radius;
		dy /= radius;

		double exp2 = std::exp(-factor*(dx*dx + dy*dy));

		return (k/a)*(exp1+exp2);

	}

	void setup_stream(
			PlaneData_Spectral &o_psi
	)
	{

		PlaneData_Physical psi_phys(o_psi.planeDataConfig);

		for (int j = 0; j < simVars.disc.space_res_physical[1]; j++)
		{
			for (int i = 0; i < simVars.disc.space_res_physical[0]; i++)
			{

				// h - lives in the center of the cell
				double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

				psi_phys.physical_set_value(j, i, stream(x, y));
			}
		}

		o_psi.loadPlaneDataPhysical(psi_phys);

	}


public:
	SWE_bench_MergeVortex(
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

		PlaneData_Spectral psi(o_h.planeDataConfig);

		/*
		 * Prepare laplace operator
		 */
		PlaneData_Spectral laplace = op.diff2_c_x + op.diff2_c_y;


		/*
		 * Setup stream function
		 */
		setup_stream(psi);
		//psi.file_physical_saveData_ascii("ouput_stream");

		//Calculate Velocities
		o_u = op.diff_c_y(psi);
		o_v = -op.diff_c_x(psi);

		//Calculate vorticity
		PlaneData_Spectral vort = op.vort(o_u, o_v);

		//Solve Poisson equation for height to get balance initial condition
		PlaneData_Spectral lap_h = (f/g)*vort;
		//o_h = lap_h.spectral_div_element_wise(laplace);
		o_h = lap_h / laplace;

	}


};


#endif
