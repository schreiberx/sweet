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
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include "PDESWEPlaneBench_BaseInterface.hpp"


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
class PDESWEPlaneBench_MergeVortex	:
		public PDESWEPlaneBench_BaseInterface
{
	double f;
	double g;
	double sx;
	double sy;

	double stream(
			double x,
			double y
			)
	{

		//double radius = shackDict.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
		double radius = 1; //shackDict.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
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
			sweet::PlaneData_Spectral &o_psi
	)
	{

		sweet::PlaneData_Physical psi_phys(o_psi.planeDataConfig);

		for (int j = 0; j < shackPlaneDataOps->space_res_physical[1]; j++)
		{
			for (int i = 0; i < shackPlaneDataOps->space_res_physical[0]; i++)
			{
				// h - lives in the center of the cell
				double x = (((double)i+0.5)/(double)shackPlaneDataOps->space_res_physical[0])*shackPlaneDataOps->plane_domain_size[0];
				double y = (((double)j+0.5)/(double)shackPlaneDataOps->space_res_physical[1])*shackPlaneDataOps->plane_domain_size[1];

				psi_phys.physical_set_value(j, i, stream(x, y));
			}
		}

		o_psi.loadPlaneDataPhysical(psi_phys);

	}


public:
	bool setupBenchmark(
			sweet::PlaneData_Spectral &o_h,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v
	)
	{
		f = shackPDESWEPlane->plane_rotating_f0;
		g = shackPDESWEPlane->gravitation;
		sx = shackPlaneDataOps->plane_domain_size[0];
		sy = shackPlaneDataOps->plane_domain_size[1];

		sweet::PlaneData_Spectral psi(o_h.planeDataConfig);

		/*
		 * Prepare laplace operator
		 */
		sweet::PlaneData_Spectral laplace = ops->diff2_c_x + ops->diff2_c_y;


		/*
		 * Setup stream function
		 */
		setup_stream(psi);
		//psi.file_physical_saveData_ascii("ouput_stream");

		//Calculate Velocities
		o_u = ops->diff_c_y(psi);
		o_v = -ops->diff_c_x(psi);

		//Calculate vorticity
		sweet::PlaneData_Spectral vort = ops->vort(o_u, o_v);

		//Solve Poisson equation for height to get balance initial condition
		sweet::PlaneData_Spectral lap_h = (f/g)*vort;
		o_h = lap_h.spectral_div_element_wise(laplace);

		return true;
	}


};


#endif
