/*
 * SWE_Plane_TS_l_irk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "SWE_Plane_TS_l_irk.hpp"

#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>




/**
 * Solve SWE with implicit Euler time stepping
 *
 * U_t = L U(0)
 *
 * (U(tau) - U(0)) / tau = L U(tau)
 *
 * <=> U(tau) - U(0) = L U(tau) tau
 *
 * <=> U(tau) - L tau U(tau) = U(0)
 *
 * <=> (1 - L tau) U(tau) = U(0)
 *
 * <=> (1/tau - L) U(tau) = U(0)/tau
 */
void SWE_Plane_TS_l_irk::runTimestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_irk: only constant time step size allowed (Please set --dt)");

#if SWEET_USE_PLANE_SPECTRAL_SPACE

	sweet::PlaneData_Spectral &eta0 = io_h;
	sweet::PlaneData_Spectral &u0 = io_u;
	sweet::PlaneData_Spectral &v0 = io_v;

	double alpha = 1.0/i_dt;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + shackPDESWEPlane->plane_rotating_f0*shackPDESWEPlane->plane_rotating_f0;

	double eta_bar = shackPDESWEPlane->h0;
	double g = shackPDESWEPlane->gravitation;

	sweet::PlaneData_Spectral rhs =
			(kappa/alpha) * eta0
			- eta_bar*(ops->diff_c_x(u0) + ops->diff_c_y(v0))
			- (shackPDESWEPlane->plane_rotating_f0*eta_bar/alpha) * (ops->diff_c_x(v0) - ops->diff_c_y(u0))
		;

	sweet::PlaneData_Spectral lhs = (-g*eta_bar*(ops->diff2_c_x + ops->diff2_c_y)).spectral_addScalarAll(kappa);
	io_h = rhs.spectral_div_element_wise(lhs);

	sweet::PlaneData_Spectral uh = u0 - g*ops->diff_c_x(io_h);
	sweet::PlaneData_Spectral vh = v0 - g*ops->diff_c_y(io_h);

	io_u = alpha/kappa * uh     + shackPDESWEPlane->plane_rotating_f0/kappa * vh;
	io_v =    -shackPDESWEPlane->plane_rotating_f0/kappa * uh + alpha/kappa * vh;

#else

	sweet::PlaneData_SpectralComplex eta0 = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(io_h);
	sweet::PlaneData_SpectralComplex u0 = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(io_u);
	sweet::PlaneData_SpectralComplex v0 = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(io_v);

	double alpha = 1.0/i_dt;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + shackPDESWEPlane->plane_rotating_f0*shackPDESWEPlane->plane_rotating_f0;

	double eta_bar = shackPDESWEPlane->h0;
	double g = shackPDESWEPlane->gravitation;

	sweet::PlaneData_SpectralComplex rhs =
			(kappa/alpha) * eta0
			- eta_bar*(opComplex.diff_c_x(u0) + opComplex.diff_c_y(v0))
			- (shackPDESWEPlane->plane_rotating_f0*eta_bar/alpha) * (opComplex.diff_c_x(v0) - opComplex.diff_c_y(u0))
		;

	sweet::PlaneData_SpectralComplex lhs = (-g*eta_bar*(opComplex.diff2_c_x + opComplex.diff2_c_y)).spectral_addScalarAll(kappa);
	sweet::PlaneData_SpectralComplex eta = rhs.spectral_div_element_wise(lhs);

	sweet::PlaneData_SpectralComplex uh = u0 - g*opComplex.diff_c_x(eta);
	sweet::PlaneData_SpectralComplex vh = v0 - g*opComplex.diff_c_y(eta);

	sweet::PlaneData_SpectralComplex u1 = alpha/kappa * uh     + shackPDESWEPlane->plane_rotating_f0/kappa * vh;
	sweet::PlaneData_SpectralComplex v1 =    -shackPDESWEPlane->plane_rotating_f0/kappa * uh + alpha/kappa * vh;

	io_h = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert(eta);
	io_u = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert(u1);
	io_v = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert(v1);
#endif
}



/*
 * Setup
 */
bool SWE_Plane_TS_l_irk::setup(
	sweet::PlaneOperators *io_ops,
	int i_order
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_irk");

	assert(i_order > 0);
	timestepping_order = i_order;

	if (timestepping_order != 1)
		SWEETError("SWE_Plane_TS_l_irk: Only 1st order IRK is supported. Please set --timestepping-order 1.");

	return true;
}




/*
 * Setup
 */
bool SWE_Plane_TS_l_irk::setup(
	sweet::PlaneOperators *io_ops
)
{
	return setup(io_ops, shackPDESWETimeDisc->timestepping_order);
}
